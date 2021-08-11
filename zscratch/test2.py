import numpy as np
from numpy.core.numeric import full

from ldpc.mod2 import row_echelon,reduced_row_echelon,rank

def rbm(m,n):

    id=np.identity(n).astype(int)

    mat=np.vstack([id,np.ones((m-n,n)).astype(int)])

    for row in range(m):

        for _ in range(np.random.randint(1,m)):
            rand_row=np.random.randint(m)
            if rand_row!=row:
                mat[row]=(mat[row]+mat[rand_row])%2

    return mat

mat=rbm(20,10)

# print(mat)
print(rank(mat))

print(row_echelon(mat,full=True)[0])


