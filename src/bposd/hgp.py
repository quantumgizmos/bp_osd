import numpy as np
from ldpc.mod2 import rank
from ldpc.code_util import compute_code_distance
from ldpc.alist import save_alist
from bposd.css import css_code

class hgp(css_code):
    def __init__(self,h1,h2=None,compute_distance=False):

        super().__init__()
        #default to symmetric HGP if

        if type(h2)!=np.ndarray: h2=np.copy(h1)

        self.h1=h1
        self.h2=h2

        #setting up base codes
        self.m1,self.n1=np.shape(h1)
        i_m1=np.identity(self.m1,dtype=int)
        i_n1=np.identity(self.n1,dtype=int)
        self.r1=rank(self.h1)
        self.k1=self.n1-self.r1
        self.k1t=self.m1-self.r1

        self.m2,self.n2=np.shape(h2)
        i_m2=np.identity(self.m2,dtype=int)
        i_n2=np.identity(self.n2,dtype=int)
        self.r2=rank(self.h2)
        self.k2=self.n2-self.r2
        self.k2t=self.m2-self.r2

        #hgp code params
        self.N=self.n1*self.n2 + self.m1*self.m2
        self.K=self.k1*self.k2+self.k1t*self.k2t #number of logical qubits in hgp code
        self.D=None

        #construct hx and hz
        self.hx1=np.kron(self.h1,i_n2)
        self.hx2=np.kron(i_m1,self.h2.T)
        self.hx = np.hstack([   self.hx1,  self.hx2 ])

        self.hz1=np.kron(i_n1,self.h2)
        self.hz2=np.kron(self.h1.T,i_m2)
        self.hz = np.hstack([   self.hz1,  self.hz2 ])

        #construct the hgp logicals
        self.compute_logicals()

        ##compute code distance if the base codes are small enough for it to be tractable
        if(compute_distance==True):
            self.d1=compute_code_distance(self.h1)
            self.d1t=compute_code_distance(self.h1.T)
            self.d2=compute_code_distance(self.h2)
            self.d2t=compute_code_distance(self.h2.T)
            self.D=np.min([self.d1,self.d1t,self.d2,self.d2t]).astype(int)
        else: self.D=None

    def print_code_parameters(self):

        if self.D==None: print( f"[[{self.N},{self.K},d]]")
        else: print( f"[[{self.N},{self.K},{self.D}]]")


class hgp_single(hgp):
    def __init__(self,h1,compute_distance=False):
        super().__init__(h1,compute_distance=compute_distance)









