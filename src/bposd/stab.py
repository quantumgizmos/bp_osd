import numpy as np
from ldpc import mod2
from tqdm import tqdm

def gf2_to_gf4(bin):
    n=int(len(bin)/2)
    gf4=np.zeros(n).astype(int)
    for i in range(n):
        if bin[i]==1 and bin[i+n]==0:
            gf4[i]=1
        elif bin[i]==0 and bin[i+n]==1:
            gf4[i]=3
        elif bin[i]==1 and bin[i+n]==1:
            gf4[i]=2
        else:
            gf4[i]=0
    return gf4


class stab_code():

    def __init__(self,hx=None,hz=None,name=None):
        if name is None:
            self.name = "<Unamed stabiliser code>"
        else: self.name=name

        if hz is None or hz is None:
            self.hx=np.array([[]])
            self.hz = np.array([[]])
            self.N=np.nan
            self.K=np.nan
            self.lx=np.array([[]])
            self.lz=np.array([[]])
            self.D=np.nan

        else:
            self.hx=hx
            self.hz=hz
            self.init_code()

        self.h=np.hstack([self.hx,self.hz])
        self.l=np.hstack([self.lx,self.lz])

    def init_code(self):

        self.h=np.hstack([self.hx,self.hz])
        self.N=self.hx.shape[1]
        self.K=self.N-mod2.rank(self.h)
        self.compute_logical_operators()
        self.D=np.nan

    def compute_logical_operators(self):
        #compute logical operators
        #Kernel H

        ker_H=mod2.nullspace(np.hstack([self.hz,self.hx]))
        image_HT=mod2.row_basis(np.hstack([self.hx,self.hz]))

        log_stack=np.vstack([image_HT,ker_H])
        pivots=mod2.row_echelon(log_stack.T)[3]
        log_op_indices=[i for i in range(image_HT.shape[0],log_stack.shape[0]) if i in pivots]
        self.l=log_stack[log_op_indices]
        self.lx =self.l[:,0:self.N]
        self.lz=self.l[:,self.N:2*self.N]

        self.K=int(self.l.shape[0]/2)

    def compute_code_distance(self,return_logicals=False):

        if self.N>10:
            print("Warning: computing a code distance of codes with N>10 will take a long time.")

        re,r,_,_=mod2.row_echelon(self.h)
        stab_basis=re[0:r]
        logical_stack=np.vstack([stab_basis,self.l])
        all_logicals=mod2.row_span(logical_stack)
        # np.argmin(np.sum(all_logicals,axis=1))
        
        d_min=self.N
        min_indices=[]
        min_logicals=[]
        for i in tqdm(range(len(all_logicals))):
            logical=all_logicals[i]
            logical=gf2_to_gf4(logical)
            temp=np.count_nonzero(logical)
            if temp<d_min:
                d_min=temp
                min_indices=[i]
                min_logicals=[logical]
            elif temp==d_min:
                min_indices.append(i)
                min_logicals.append(logical)
        
        # d_min=np.min( np.sum(all_logicals,axis=1) )
        self.D=d_min

        # print(all_logicals)

        if return_logicals:
            return np.array(min_logicals)

        return d_min


    def test(self, show_tests=True):
        valid_code=True

        # if self.K==np.nan: self.compute_dimension()

        code_label=f"{self.code_params}"

        if show_tests: print(f"{self.name}, {code_label}")

        try:
            assert self.N==self.hz.shape[1]==self.lz.shape[1]==self.lx.shape[1]
            assert self.K==self.lz.shape[0]//2==self.lx.shape[0]//2
            if show_tests: print(" -Block dimensions: Pass")
        except AssertionError:
            valid_code=False
            print(" -Block dimensions incorrect")

        try:
            assert (( (self.hz@self.hx.T %2) + (self.hx@self.hz.T %2) ) %2 ).any() == 0
            if show_tests: print(" -PCMs commute hz@hx.T==0: Pass")
        except AssertionError:
            valid_code=False
            print(" -PCMs commute hz@hx.T==0: Fail")

        # if show_tests and valid_code: print("\t-PCMs commute hx@hz.T == hz@hx.T ==0: Pass")

        try:
            assert ((self.hx@self.lz.T %2 + self.hz@self.lx.T %2)%2).any()==0
        except AssertionError:
            valid_code=False
            print(" -lx \in ker{hz} AND lz \in ker{hx}: Fail")


        # if show_tests and valid_code: print("\t-lx \in ker{hz} AND lz \in ker{hx}: Pass")

        try:
            assert mod2.rank((self.lx@self.lz.T%2 + self.lz@self.lx.T%2) %2)==self.l.shape[0]
            if show_tests: print(" -lx and lz anticommute: Pass")
        except AssertionError:
            valid_code=False
            print(" -lx and lz anitcommute: Fail")

        # if show_tests and valid_code: print("\t- lx and lz anitcommute: Pass")

        if show_tests and valid_code:
            print(f"{self.name} is a valid stabiliser code w/ params {code_label}")

        return valid_code

    @property
    def code_params(self):
        return f"[[{self.N},{self.K},{self.D}]]"

    