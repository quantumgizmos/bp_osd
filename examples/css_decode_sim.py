import numpy as np
from tqdm import tqdm
import json
import sys
import time
import datetime

from bposd.bposd import bposd_decoder
from qec.css import css_code
from qec.mod2 import rank
from qec.alist import alist2numpy


def classical_decode_sim(
        pcm,
        error_rate=0,
        max_iter=0,
        target_runs=100,
        seed=0,
        bp_method="ps",
        ms_scaling_factor=0.625,
        output_file=None,
        save_interval=3,
        error_bar_precision_cutoff=1e-2,
        output_dict={}
):

    output_dict['run_type']='classical_decode_sim'

    start_date=datetime.datetime.fromtimestamp(time.time()).strftime("%A, %B %d, %Y %H:%M:%S")
    output_dict.update(dict(start_date=start_date))

    if error_rate<=0:
        raise Exception("Error: The error_rate input variable to simulate_bp_osd is not specified. It must have a value in the range 0.0<error_rate<1.0")

    m,n=pcm.shape
    if max_iter==0: max_iter=n

    output_dict['n']=n
    k=n-rank(pcm)
    output_dict['k']=k

    bpd=bposd_decoder(
        pcm,
        error_rate,
        max_iter=max_iter,
        bp_method=bp_method,
        ms_scaling_factor=ms_scaling_factor
    )


    print(f"BP Method: {bpd.bp_method}")

    if seed==0: seed=np.random.randint(low=1,high=2**32-1)
    np.random.seed(seed)
    print(f"RNG Seed: {seed}")

    if bpd.bp_method=="mininum_sum":
        output_dict['ms_scaling_factor']=ms_scaling_factor


    bp_success=0
    bp_converge_count=0

    output_dict.update(dict(seed=seed,target_runs=target_runs,error_rate=error_rate,max_iter=max_iter,bp_method=bp_method))

    start_time = time.time()
    save_time=start_time
    pbar=tqdm(range(1,target_runs+1),disable=False,ncols=0)
    error=np.zeros(n).astype(int)

    for run_count in pbar:

        for j in range(n):
            if np.random.random()<error_rate:
                error[j]=1
            else: error[j]=0

        syndrome=pcm@error%2

        bp_decoding=bpd.decode(syndrome)

        if bpd.converge:
            bp_converge_count+=1
            if np.array_equal(bp_decoding,error): bp_success+=1

        bp_frame_error_rate=1.0-bp_success/run_count
        bp_frame_error_rate_eb=np.sqrt((1-bp_frame_error_rate)*bp_frame_error_rate/run_count)
        bp_converge_failure_rate=1-bp_converge_count/run_count

        pbar.set_description(f"BP Error Rate: {bp_frame_error_rate*100:.3g}±{bp_frame_error_rate_eb*100:.2g}%")

        current_time=time.time()
        elapsed_time=current_time-start_time
        save_loop=current_time-save_time

        if int(save_loop)>save_interval or run_count==target_runs:
            save_time=time.time()
            output_dict.update(dict(bp_success_count=bp_success,bp_converge_failure_rate=bp_converge_failure_rate,run_count=run_count,bp_frame_error_rate=bp_frame_error_rate,bp_frame_error_rate_eb=bp_frame_error_rate_eb))
            output_dict['elapsed_sim_time']=time.strftime('%H:%M:%S', time.gmtime(elapsed_time))


            if output_file!=None:
                f=open(output_file,"w+")
                print(json.dumps(output_dict,sort_keys=True, indent=4),file=f)
                f.close()

            if bp_frame_error_rate_eb>0 and bp_frame_error_rate_eb/bp_frame_error_rate < error_bar_precision_cutoff:
                print("\nTarget error bar precision reached. Stopping simulation...")
                break


    return json.dumps(output_dict,sort_keys=True, indent=4)




def css_decode_sim(
        hx,
        hz,
        lx=[None],
        lz=[None],
        error_rate=0,
        max_iter=0,
        target_runs=100,
        seed=0,
        bp_method="ms",
        ms_scaling_factor=0.625,
        osd_order=40,
        osd_method="osd_cs",
        noise_type="x",
        output_file=None,
        save_interval=3,
        error_bar_precision_cutoff=1e-2,
        output_dict={},
        check_code=True,
        tqdm_disable=False
):


    start_date=datetime.datetime.fromtimestamp(time.time()).strftime("%A, %B %d, %Y %H:%M:%S")
    output_dict['start_date']=start_date

    output_dict['run_type']="css_decode_sim"




    if error_rate<=0:
        raise Exception("Error: The error_rate input variable to simulate_bp_osd is not specified. It must have a value in the range 0.0<error_rate<1.0")

    print("Constructing CSS code from hx and hz matrices...")
    if isinstance(hx,np.ndarray) and isinstance(lx,list):
        qcode=css_code(hx,hz)
        lx=qcode.lx
        lz=qcode.lz
        output_dict['K']=qcode.K
        output_dict['N']=qcode.N

        print("Checking the CSS code is valid...")
        if check_code and not qcode.test():
            raise Exception("Error: invalid CSS code. Check the form of your hx and hz matrices!")

    if isinstance(hx,np.ndarray) and isinstance(lx,np.ndarray):
        output_dict['K']=rank(lx)
        output_dict['N']=hx.shape[1]

    elif isinstance(hx,str):
        hx=alist2numpy(hx)
        hz=alist2numpy(hz)
        lx=alist2numpy(lx)
        lz=alist2numpy(lz)
        output_dict['K']=lx.shape[0]
        output_dict['N']=hx.shape[1]

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

    if bpd.bp_method=="mininum_sum":
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
        
        bp_converge_failure_rate=1-bp_converge_count/run_count
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
            output_dict.update(dict(bp_converge_failure_rate=bp_converge_failure_rate,osdw_logical_error_rate=osdw_logical_error_rate,osd0_logical_error_rate=osd0_logical_error_rate,bp_logical_error_rate=bp_logical_error_rate,osdw_logical_error_rate_eb=osdw_logical_error_rate_eb,osd0_logical_error_rate_eb=osd0_logical_error_rate_eb,bp_logical_error_rate_eb=bp_logical_error_rate_eb))
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


from qec.hgp import hgp_single
from qec.quantum_codes import surface_code, toric_code

h=np.loadtxt("examples/mkmn_16_4_6.txt").astype(int)

hgp=hgp_single(h) # construct quantum code using the symmetric hypergraph product


output_dict={"code_type": "hgp_mkmn_16_4_6"}
output_dict=css_decode_sim(
    hx=hgp.hx, # CSS hx matrix
    hz=hgp.hz, # CSS hz matrix
    error_rate=0.05, # the physical error rate
    max_iter=0, # If `0` iterations depth is set to block length n
    target_runs=1000, # the number of runs to simulate
    seed=42, # random seed
    bp_method="ms", # Options: 1) 'product_sum'; 2) 'mininum_sum'
    ms_scaling_factor=0,
    osd_method=2, # Options: 1) 'exhaustive'; 2) combination sweep
    osd_order=40, # The OSD order
    noise_type='z', # Options: 1) 'z'; 2) 'x'
    output_file="hpg_bp_osd_decode_sim_output.json", #output directory
    output_dict=output_dict, # input dictionary. Script will overwrite existing fields. Good for adding meta_data
    check_code=True
)






    



    