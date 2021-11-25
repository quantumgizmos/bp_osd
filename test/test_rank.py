from bposd.bposd import m2s_rank
from ldpc.codes import ring_code
import numpy as np
from ldpc import mod2


for i in range(100):
    mat=np.random.randint(low=0,high=2,size=(np.random.randint(3,100),np.random.randint(3,100)))
    m,n=mat.shape

    for j in range(np.random.randint(np.min([m,n]))):
        mat[np.random.randint(m),:]=np.zeros(n).astype(int)
        mat[:,np.random.randint(n)]=np.zeros(m).astype(int)

    # print(mat)

    # print(m,n,mod2.rank(mat))

    assert mod2.rank(mat)==m2s_rank(mat)