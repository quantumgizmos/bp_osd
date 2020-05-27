import numpy as np

def save_alist(name, mat, j=None, k=None):

    H=np.copy(mat)
    H=H.T

    '''
    Function converts parity check matrix into the format required for the RN decoder
    '''

    if j is None:
        j=int(max(H.sum(axis=0)))


    if k is None:
        k=int(max(H.sum(axis=1)))


    m, n = H.shape # rows, cols
    f = open(name, 'w')
    print(n, m, file=f)
    print(j, k, file=f)

    for col in range(n):
        print( int(H[:, col].sum()), end=" ", file=f)
    print(file=f)
    for row in range(m):
        print( int(H[row, :].sum()), end=" ", file=f)
    print(file=f)

    for col in range(n):
        for row in range(m):
            if H[row, col]:
                print( row+1, end=" ", file=f)
        print(file=f)

    for row in range(m):
        for col in range(n):
            if H[row, col]:
                print(col+1, end=" ", file=f)
        print(file=f)
    f.close()