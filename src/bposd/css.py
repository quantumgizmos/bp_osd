import numpy as np
from ldpc import mod2
from ldpc.alist import save_alist
from ldpc.code_util import compute_code_distance
from bposd import stab

class css_code():

    def __init__(self,hx=np.array([[]]),hz=np.array([[]]),code_distance=np.nan, name="<Unnamed CSS code>"):

        self.hx=hx #hx pcm
        self.hz=hz #hz pcm

        self.lx=np.array([[]]) #x logicals
        self.lz=np.array([[]]) #z logicals

        self.N=np.nan #block length
        self.K=np.nan #code dimension
        self.D=code_distance #code distance
        self.L=np.nan #max column weight
        self.Q=np.nan #max row weight

        _,nx=self.hx.shape
        _,nz=self.hz.shape
        try:
            assert nx==nz
        except AssertionError:
            raise Exception("Error: hx and hz matrices must have equal numbers of columns!")

        if nx!=0:
            self.compute_dimension()
            self.compute_ldpc_params()
            self.compute_logicals()
            if code_distance==0:
                dx=compute_code_distance(hx)
                dz=compute_code_distance(hz)
                self.D=np.min([dx,dz])
                




        self.name=name

    def compute_dimension(self):

        self.N=self.hx.shape[1]
        assert self.N == self.hz.shape[1], "Code block length (N) inconsistent!"

        self.K=self.N-mod2.rank(self.hx)-mod2.rank(self.hz)
        return self.K

    def compute_ldpc_params(self):

        #column weights
        hx_l=np.max(np.sum(self.hx,axis=0))
        hz_l=np.max(np.sum(self.hz,axis=0))
        self.L=np.max([hx_l,hz_l]).astype(int)

        #row weights
        hx_q=np.max(np.sum(self.hx,axis=1))
        hz_q=np.max(np.sum(self.hz,axis=1))
        self.Q=np.max([hx_q,hz_q]).astype(int)

    def save_sparse(self, code_name):

        self.code_name=code_name

        hx=self.hx
        hz=self.hz
        save_alist(f"{code_name}_hx.alist",hx)
        save_alist(f"{code_name}_hz.alist",hz)

        lx=self.lx
        lz=self.lz
        save_alist(f"{code_name}_lx.alist",lx)
        save_alist(f"{code_name}_lz.alist",lz)

    def to_stab_code(self):

        hx=np.vstack([np.zeros(self.hz.shape,dtype=int),self.hx])
        hz=np.vstack([self.hz,np.zeros(self.hx.shape,dtype=int)])
        return stab.stab_code(hx,hz)

    @property
    def h(self):
        hx=np.vstack([np.zeros(self.hz.shape,dtype=int),self.hx])
        hz=np.vstack([self.hz,np.zeros(self.hx.shape,dtype=int)])
        return np.hstack([hx,hz])

    @property
    def l(self):
        lx=np.vstack([np.zeros(self.lz.shape,dtype=int),self.lx])
        lz=np.vstack([self.lz,np.zeros(self.lx.shape,dtype=int)])
        return np.hstack([lx,lz])


    def compute_code_distance(self):
        temp=self.to_stab_code()
        self.D=temp.compute_code_distance()
        return self.D

    def compute_logicals(self):

        def compute_lz(hx,hz):
            #lz logical operators
            #lz\in ker{hx} AND \notin Im(Hz.T)

            ker_hx=mod2.nullspace(hx) #compute the kernel basis of hx
            im_hzT=mod2.row_basis(hz) #compute the image basis of hz.T

            #in the below we row reduce to find vectors in kx that are not in the image of hz.T.
            log_stack=np.vstack([im_hzT,ker_hx])
            pivots=mod2.row_echelon(log_stack.T)[3]
            log_op_indices=[i for i in range(im_hzT.shape[0],log_stack.shape[0]) if i in pivots]
            log_ops=log_stack[log_op_indices]
            return log_ops

        if self.K==np.nan: self.compute_dimension()
        self.lx=compute_lz(self.hz,self.hx)
        self.lz=compute_lz(self.hx,self.hz)

        return self.lx,self.lz

    def canonical_logicals(self):
        temp=mod2.inverse(self.lx@self.lz.T %2)
        self.lx=temp@self.lx %2


    @property
    def code_params(self):
        try: self.N
        except AttributeError: self.N=np.nan
        try: self.K
        except AttributeError: self.K=np.nan
        try: self.D
        except AttributeError: self.D=np.nan
        try: self.L
        except AttributeError: self.L=np.nan
        try: self.Q
        except AttributeError: self.Q=np.nan

        return f"({self.L},{self.Q})-[[{self.N},{self.K},{self.D}]]"

    def test(self, show_tests=True):
        valid_code=True

        if self.K==np.nan: self.compute_dimension()
        self.compute_ldpc_params()

        code_label=f"{self.code_params}"

        if show_tests: print(f"{self.name}, {code_label}")

        try:
            assert self.N==self.hz.shape[1]==self.lz.shape[1]==self.lx.shape[1]
            assert self.K==self.lz.shape[0]==self.lx.shape[0]
            if show_tests: print(" -Block dimensions: Pass")
        except AssertionError:
            valid_code=False
            print(" -Block dimensions incorrect")

        try:
            assert not (self.hz@self.hx.T %2).any()
            if show_tests: print(" -PCMs commute hz@hx.T==0: Pass")
        except AssertionError:
            valid_code=False
            print(" -PCMs commute hz@hx.T==0: Fail")

        try:
            assert not (self.hx@self.hz.T %2).any()
            if show_tests: print(" -PCMs commute hx@hz.T==0: Pass")
        except AssertionError:
            valid_code=False
            print(" -PCMs commute hx@hz.T==0: Fail")

        # if show_tests and valid_code: print("\t-PCMs commute hx@hz.T == hz@hx.T ==0: Pass")

        try:
            assert not (self.hz@self.lx.T %2).any()
        except AssertionError:
            valid_code=False
            print(" -lx \in ker{hz} AND lz \in ker{hx}: Fail")


        try:
            assert not (self.hx@self.lz.T %2).any()
            if show_tests: print(" -lx \in ker{hz} AND lz \in ker{hx}: Pass")
        except AssertionError:
            valid_code=False
            print(" -lx \in ker{hz} AND lz \in ker{hx}: Fail")


        # if show_tests and valid_code: print("\t-lx \in ker{hz} AND lz \in ker{hx}: Pass")

        try:
            assert mod2.rank(self.lx@self.lz.T %2)==self.K
            if show_tests: print(" -lx and lz anticommute: Pass")
        except AssertionError:
            valid_code=False
            print(" -lx and lz anitcommute: Fail")

        # if show_tests and valid_code: print("\t- lx and lz anitcommute: Pass")

        if show_tests and valid_code: print(f" -{self.name} is a valid CSS code w/ params {code_label}")

        return valid_code

