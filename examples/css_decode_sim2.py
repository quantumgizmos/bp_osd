import numpy as np
from tqdm import tqdm
import json
import time
import datetime
from copy import copy

from bposd import bposd_decoder
from bposd.css import css_code


class css_decode_sim2():
    def __init__(self,hx=None,hz=None,**input_dict):
    
        default_input={
            'error_rate': None,
            'xyz_error_bias': [1,1,1], 
            'channel_probs_x': [None],
            'channel_probs_z': [None],
            'max_iter': 0,
            'target_runs': 100,
            'seed':0,
            'bp_method':"minimum_sum",
            'ms_scaling_factor': 1.0,
            'osd_order': 2,
            'osd_method': "osd_cs",
            'output_file': None,
            'check_code': True,
            'tqdm_disable': False,
            'run_sim': True
            }

        for key in input_dict.keys():
            self.__dict__[key]=input_dict[key]

        for key in default_input.keys():
            if key not in input_dict:
                self.__dict__[key]=default_input[key]

        output_values={
            "K": None,
            "N":None,
            "start_date": None,
            "runtime":0,
            "runtime_readable":None,
            "run_count":0,
            "bp_converge_count":0,
            "bp_success_count":0,
            "bp_logical_error_rate":0,
            "bp_logical_error_rate_eb":0, 
            "osd0_success_count":0,
            "osd0_logical_error_rate":0,
            "osd0_logical_error_rate_eb":0,
            "osdw_success_count":0,
            "osdw_logical_error_rate":0,
            "osdw_logical_error_rate_eb":0,
            "osdw_word_error_rate":0,
            "osdw_word_error_rate_eb":0,
        }

        for key in output_values.keys():
            if key not in self.__dict__:
                self.__dict__[key]=output_values[key]

        temp=[]
        for key in self.__dict__.keys():
            temp.append(key)
        self.output_keys=temp

        self.hx=hx
        self.hz=hz

        print("Constructing CSS code from hx and hz matrices...")
        if isinstance(self.hx,np.ndarray) and isinstance(self.hz,np.ndarray):
            
            qcode=css_code(self.hx,self.hz)
            self.lx=qcode.lx
            self.lz=qcode.lz
            self.K=qcode.K
            self.N=qcode.N

            print("Checking the CSS code is valid...")
            if self.check_code and not qcode.test():
                raise Exception("Error: invalid CSS code. Check the form of your hx and hz matrices!")
        else:
            raise Exception("Invalid object type for the hx/hz matrices")


        xyz_error_bias=np.array(self.xyz_error_bias)

        if xyz_error_bias[0]==np.inf:
            self.px=self.error_rate
            self.py=0
            self.pz=0
        elif xyz_error_bias[1]==np.inf:
            self.px=0
            self.py=self.error_rate
            self.pz=0
        elif xyz_error_bias[2]==np.inf:
            self.px=0
            self.py=0
            self.pz=self.error_rate
        else:
            self.px,self.py,self.pz=self.error_rate*xyz_error_bias/np.sum(xyz_error_bias)

        self.channel_probs_x=np.ones(self.N)*(self.px+self.py)
        self.channel_probs_z=np.ones(self.N)*(self.pz+self.py)

        if self.run_sim:
            self.run_decode_sim()

    def output_dict(self):
        output_dict={}
        for key, value in self.__dict__.items():
            if key in self.output_keys:
                output_dict[key]=value
        return json.dumps(output_dict,sort_keys=True, indent=4)

    def run_decode_sim(self):
    
        if np.sum(self.channel_probs_x)==0: self.INFINITE_Z_BIAS=True
        else: self.INFINITE_Z_BIAS=False

        if np.sum(self.channel_probs_z)==0: self.INFINITE_X_BIAS=True
        else: self.INFINITE_X_BIAS=False

        print(self.INFINITE_X_BIAS)
        print(self.INFINITE_Z_BIAS)
        print(self.px,self.py,self.pz)
        print('x',self.channel_probs_x)
        print('z',self.channel_probs_z)

        #save start date
        self.start_date=datetime.datetime.fromtimestamp(time.time()).strftime("%A, %B %d, %Y %H:%M:%S")
        
        bpd_z=bposd_decoder(
            self.hx,
            channel_probs=self.channel_probs_z,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            ms_scaling_factor=self.ms_scaling_factor,
            osd_method=self.osd_method,
            osd_order=self.osd_order,
            )

        bpd_x=bposd_decoder(
            self.hz,
            channel_probs=self.channel_probs_x,
            max_iter=self.max_iter,
            bp_method=self.bp_method,
            ms_scaling_factor=self.ms_scaling_factor,
            osd_method=self.osd_method,
            osd_order=self.osd_order,
            )

        pbar=tqdm(range(self.run_count+1,self.target_runs+1),disable=self.tqdm_disable,ncols=0)

        start_time = time.time()
        save_time=start_time
        min_logical_weight=self.N

        error_x=np.zeros(self.N).astype(int)
        error_z=np.zeros(self.N).astype(int)

        for run_count in pbar:

            # print(run_count)
            current_time=time.time()
            elapsed_time=current_time-start_time
            save_loop=current_time-save_time

            #generate error
            for i in range(self.N):
                rand=np.random.random()
                if rand<self.pz:
                    error_z[i]=1
                    error_x[i]=0
                elif self.pz<=rand<(self.pz+self.px):
                    error_z[i]=0
                    error_x[i]=1
                elif (self.pz+self.px)<=rand<(self.px+self.py+self.pz):
                    error_z[i]=1
                    error_x[i]=1
                else:
                    error_z[i]=0
                    error_x[i]=0
            

            # decode z
            if not self.INFINITE_X_BIAS:
                synd_z=self.hx@error_z %2
                bpd_z.decode(synd_z)

            #decode x
            if not self.INFINITE_Z_BIAS:
                synd_x=self.hz@error_x %2
                bpd_x.decode(synd_x)

            ### osdw - checking for logical errors
            #calculate the residual error
      
            residual_x=(error_x+bpd_x.osdw_decoding)%2
            residual_z=(error_z+bpd_z.osdw_decoding)%2
            if (self.lz@residual_x %2).any() and not self.INFINITE_Z_BIAS:
                logical_weight=np.sum(residual_x)
                if logical_weight<min_logical_weight:
                    min_logical_weight=logical_weight
            elif (self.lx@residual_z %2).any() and not self.INFINITE_X_BIAS:
                logical_weight=np.sum(residual_z)
                if logical_weight<min_logical_weight:
                    min_logical_weight=logical_weight
            else: self.osdw_success_count+=1
            self.osdw_logical_error_rate=1-self.osdw_success_count/run_count
            self.osdw_logical_error_rate_eb=np.sqrt((1-self.osdw_logical_error_rate)*self.osdw_logical_error_rate/run_count)

            self.osdw_word_error_rate=1.0 - (1-self.osdw_logical_error_rate)**(1/self.K)
            self.osdw_word_error_rate_eb=self.osdw_logical_error_rate_eb*((1-self.osdw_logical_error_rate_eb)**(1/self.K -1))/self.K

            # osdw_zx_logical_error_rate=osdw_zx_logical_errors/run_count
            # osdw_zx_logical_error_rate_eb=np.sqrt((1-osdw_zx_logical_error_rate)*osdw_zx_logical_error_rate/run_count)
            # osdw_xz_logical_error_rate=osdw_xz_logical_errors/run_count
            # osdw_xz_logical_error_rate_eb=np.sqrt((1-osdw_xz_logical_error_rate)*osdw_xz_logical_error_rate/run_count)


        

            pbar.set_description(
f"d_max: {min_logical_weight}; \
Error rate:OSDW_WER: {self.osdw_word_error_rate*100:.3g}±{self.osdw_word_error_rate_eb*100:.2g}%; \
OSDW: {self.osdw_logical_error_rate*100:.3g}±{self.osdw_logical_error_rate_eb*100:.2g}%; \
OSD0: {self.osd0_logical_error_rate*100:.3g}±{self.osd0_logical_error_rate_eb*100:.2g}%;")


        #     if int(save_loop)>save_interval or run_count==target_runs:
        #         save_time=time.time()
        #         output_dict.update(dict(run_count=run_count))
        #         output_dict.update(dict(observed_min_distance=int(min_logical_weight)))
        #         output_dict.update(dict(osdw_logical_error_rate=osdw_logical_error_rate,osd0_logical_error_rate=osd0_logical_error_rate,osdw_logical_error_rate_eb=osdw_logical_error_rate_eb,osd0_logical_error_rate_eb=osd0_logical_error_rate_eb))
        #         output_dict.update(dict(elapsed_sim_time=time.strftime('%H:%M:%S', time.gmtime(elapsed_time))))
        #         output_dict.update(dict(osdw_success_count=osdw_success_count,osd0_success_count=osd0_success_count))

        #         if output_file!=None:
        #             f=open(output_file,"w+")
        #             print(json.dumps(output_dict,sort_keys=True, indent=4),file=f)
        #             f.close()

        #         if osdw_logical_error_rate_eb>0 and osdw_logical_error_rate_eb/osdw_logical_error_rate < error_bar_precision_cutoff:
        #             print("\nTarget error bar precision reached. Stopping simulation...")
        #             break

        # print(f"{100*round(1-osdw_success_count/run_count,3)}%")
        

        # return json.dumps(self.output_dict(),sort_keys=True, indent=4)

def css_decode_sim(
        hx, #x stabilisers
        hz, #z stabilisers
        error_rate=0, #error rate on each qubit
        max_iter=0, #maximum number of iterations for BP( if ==0: max_iter=block_length)
        target_runs=100, #number of runs to simulate
        seed=0, #seed for rng. Automatically selected if set to 0/
        bp_method="ms", #bp_method. Choose "minimum_sum" or "product_sum"
        ms_scaling_factor=0.625, #If applicable, choose a `minimum_sum` scaling factor
        osd_order=40, #The OSD search depth
        osd_method="osd_cs", #OSD method
        noise_type="x", #Noise type. Choose either "x" (bit) or "z" (phase)
        output_file=None, #output file name
        save_interval=3, #time interval for saving to the output file
        error_bar_precision_cutoff=1e-2,
        output_dict={},
        check_code=True,
        tqdm_disable=False
):

    #
    start_date=datetime.datetime.fromtimestamp(time.time()).strftime("%A, %B %d, %Y %H:%M:%S")
    output_dict['start_date']=start_date

    output_dict['run_type']="css_decode_sim"

    if error_rate<=0:
        raise Exception("Error: The error_rate input variable to simulate_bp_osd is not specified. It must have a value in the range 0.0<error_rate<1.0")

    print("Constructing CSS code from hx and hz matrices...")
    if isinstance(hx,np.ndarray) and isinstance(hz,np.ndarray):
        qcode=css_code(hx,hz)
        lx=qcode.lx
        lz=qcode.lz
        output_dict['K']=qcode.K
        output_dict['N']=qcode.N

        print("Checking the CSS code is valid...")
        if check_code and not qcode.test():
            raise Exception("Error: invalid CSS code. Check the form of your hx and hz matrices!")

    else:
        raise Exception("Invalid object type for the hx/hz matrices")


    if noise_type.lower()=="phase" or noise_type.lower()=='z':
        matrix_of_logicals=lx
        H=hx
    elif noise_type.lower()=="bit" or noise_type.lower()=="x":
        matrix_of_logicals=lz
        H=hz
    else:
        raise Exception("Error: noise-type must be either 'phase' or 'bit'")


    m,n=H.shape


    if max_iter==0:
        max_iter=n

    bpd=bposd_decoder(
        H,
        error_rate,
        max_iter=max_iter,
        bp_method=bp_method,
        ms_scaling_factor=ms_scaling_factor,
        osd_order=osd_order,
        osd_method=osd_method
    )


    print(f"BP Method: {bpd.bp_method}")
    print(f"OSD Method: {bpd.osd_method}")

    if bpd.bp_method=="minimum_sum":
        output_dict['ms_scaling_factor']=ms_scaling_factor

    if seed==0: seed=np.random.randint(low=1,high=2**32-1)
    np.random.seed(seed)


    print(f"RNG Seed: {seed}")

    bp_converge_count=0
    bp_success=0
    osd0_success=0
    osdw_success=0


    output_dict.update(dict(seed=seed,noise_type=noise_type,target_runs=target_runs,error_rate=error_rate,max_iter=max_iter,bp_method=bp_method,osd_order=osd_order,osd_method=osd_method))

    start_time = time.time()
    save_time=start_time
    pbar=tqdm(range(1,target_runs+1),disable=tqdm_disable,ncols=0)
    error=np.zeros(n).astype(int)
    for run_count in pbar:

        for j in range(n):
            if np.random.random()<error_rate:
                error[j]=1
            else: error[j]=0

        syndrome=H@error%2

        osdw_decoding=bpd.decode(syndrome)
        
        bp_decoding=bpd.bp_decoding
        bp_residual=(bp_decoding+error)%2
        if bpd.converge:
            bp_converge_count+=1
            if not (matrix_of_logicals@bp_residual % 2).any(): bp_success+=1

        osd0_decoding=bpd.osd0_decoding
        osd0_residual=(osd0_decoding+error)%2
        if bpd.converge or bpd.osd_order>-1:
            if not (matrix_of_logicals@osd0_residual % 2).any(): osd0_success+=1
        

        osdw_residual=(osdw_decoding+error)%2
        if bpd.converge or bpd.osd_order>-1:
            if not (matrix_of_logicals@osdw_residual % 2).any(): osdw_success+=1
        
        bp_converge_logical_error_rate=1-bp_converge_count/run_count
        bp_logical_error_rate=1-bp_success/run_count
        bp_logical_error_rate_eb=np.sqrt((1-bp_logical_error_rate)*bp_logical_error_rate/run_count)

        osd0_logical_error_rate=1-osd0_success/run_count
        osd0_logical_error_rate_eb=np.sqrt((1-osd0_logical_error_rate)*osd0_logical_error_rate/run_count)
        osdw_logical_error_rate=1-osdw_success/run_count
        osdw_logical_error_rate_eb=np.sqrt((1-osdw_logical_error_rate)*osdw_logical_error_rate/run_count)


        pbar.set_description(f"Error rate: OSDW: {osdw_logical_error_rate*100:.3g}±{osdw_logical_error_rate_eb*100:.2g}%; OSD0: {osd0_logical_error_rate*100:.3g}±{osd0_logical_error_rate_eb*100:.2g}%; BP: {bp_logical_error_rate*100:.3g}±{bp_logical_error_rate_eb*100:.2g}%")

        current_time=time.time()
        elapsed_time=current_time-start_time
        save_loop=current_time-save_time

        if int(save_loop)>save_interval or run_count==target_runs:
            save_time=time.time()
            output_dict.update(dict(run_count=run_count))
            output_dict.update(dict(bp_converge_logical_error_rate=bp_converge_logical_error_rate,osdw_logical_error_rate=osdw_logical_error_rate,osd0_logical_error_rate=osd0_logical_error_rate,bp_logical_error_rate=bp_logical_error_rate,osdw_logical_error_rate_eb=osdw_logical_error_rate_eb,osd0_logical_error_rate_eb=osd0_logical_error_rate_eb,bp_logical_error_rate_eb=bp_logical_error_rate_eb))
            output_dict.update(dict(elapsed_sim_time=time.strftime('%H:%M:%S', time.gmtime(elapsed_time))))
            output_dict.update(dict(bp_success_count=bp_success,osdw_success_count=osdw_success,osd0_success_count=osd0_success))

            if output_file!=None:
                f=open(output_file,"w+")
                print(json.dumps(output_dict,sort_keys=True, indent=4),file=f)
                f.close()

            if osdw_logical_error_rate_eb>0 and osdw_logical_error_rate_eb/osdw_logical_error_rate < error_bar_precision_cutoff:
                print("\nTarget error bar precision reached. Stopping simulation...")
                break


    return json.dumps(output_dict,sort_keys=True, indent=4)








    



    