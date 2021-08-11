#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

import numpy as np
from ldpc.codes import rep_code


# def test():

#     H=rep_code(3)
#     hm2=numpy2mod2sparse(H)

#     print(mod2sparse_rows(hm2))

#     cdef char* test

#     k=10
#     for i in range(k):

#         print(char2numpy(test,k))


