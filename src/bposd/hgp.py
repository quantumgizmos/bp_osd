import numpy as np
from ldpc.mod2 import rank
from ldpc.code_util import compute_exact_code_distance
from bposd.css import css_code
import scipy


class hgp(css_code):
    def __init__(self, h1, h2=None, compute_distance=False):
        super().__init__()
        # default to symmetric HGP if

        if not scipy.sparse.issparse(h1):
            h1 = scipy.sparse.csr_matrix(h1)

        if h2 is None:
            h2 = h1.copy()

        if not scipy.sparse.issparse(h2):
            h2 = scipy.sparse.csr_matrix(h2)

        self.h1 = h1
        self.h2 = h2

        # setting up base codes
        self.m1, self.n1 = h1.shape
        i_m1 = scipy.sparse.identity(self.m1, format="csr", dtype=np.uint8)
        i_n1 = scipy.sparse.identity(self.n1, format="csr", dtype=np.uint8)
        self.r1 = rank(self.h1.toarray())
        self.k1 = self.n1 - self.r1
        self.k1t = self.m1 - self.r1

        self.m2, self.n2 = h2.shape
        i_m2 = scipy.sparse.identity(self.m2, format="csr", dtype=np.uint8)
        i_n2 = scipy.sparse.identity(self.n2, format="csr", dtype=np.uint8)
        self.r2 = rank(self.h2.toarray())
        self.k2 = self.n2 - self.r2
        self.k2t = self.m2 - self.r2

        # hgp code params
        self.N = self.n1 * self.n2 + self.m1 * self.m2
        self.K = (
            self.k1 * self.k2 + self.k1t * self.k2t
        )  # number of logical qubits in hgp code
        self.D = None

        # construct hx and hz
        self.hx1 = scipy.sparse.kron(self.h1, i_n2, format="csr")
        self.hx2 = scipy.sparse.kron(i_m1, self.h2.T, format="csr")
        self.hx = scipy.sparse.hstack([self.hx1, self.hx2], format="csr")

        self.hz1 = scipy.sparse.kron(i_n1, self.h2, format="csr")
        self.hz2 = scipy.sparse.kron(self.h1.T, i_m2, format="csr")
        self.hz = scipy.sparse.hstack([self.hz1, self.hz2], format="csr")

        # construct the hgp logicals
        self.compute_logicals()

        ##compute code distance if the base codes are small enough for it to be tractable
        if compute_distance == True:
            if self.h1.shape[1] != rank(self.h1):
                self.d1 = compute_exact_code_distance(self.h1)
            else:
                self.d1 = np.inf

            if self.h2.shape[1] != rank(self.h2):
                self.d2 = compute_exact_code_distance(self.h2)
            else:
                self.d2 = np.inf

            if self.h1.T.shape[1] != rank(self.h1.T):
                self.d1t = compute_exact_code_distance(self.h1.T)
            else:
                self.d1t = np.inf

            if self.h2.T.shape[1] != rank(self.h2.T):
                self.d2t = compute_exact_code_distance(self.h2.T)
            else:
                self.d2t = np.inf

            self.D = np.min([self.d1, self.d1t, self.d2, self.d2t]).astype(int)
        else:
            self.D = None

        def print_code_parameters(self):
            if self.D == None:
                print(f"[[{self.N},{self.K},d]]")
            else:
                print(f"[[{self.N},{self.K},{self.D}]]")


class hgp_single(hgp):
    def __init__(self, h1, compute_distance=False):
        super().__init__(h1, compute_distance=compute_distance)
