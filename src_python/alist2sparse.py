import numpy as np
import scipy.sparse as sps

def alist2sparse(fname):
# reads binary parity check matrix in "alist" format from file FNAME and 
# converts it to sparse matrix used in MATLAB routines. 
# This is an interface to matrices at http://wol.ra.phy.cam.ac.uk/mackay/codes/ 
# 
# Example 
#        [H] = alist2sparse('A');   % A is the ascii file in alist format 
 
 
#   Copyright (c) 1999 by Igor Kozintsev igor@ifp.uiuc.edu 
#   $Revision: 1.1 $  $Date: 2000/03/23 $ Bug fixed by Hatim Behairy 
# matlab->python  2020-05-28 14:26:11 lingr7@qq.com

    a = np.loadtxt(fname, delimiter='\n',dtype=str)
    alist_m = a.shape[0]
    list_a = []
    for i in range(alist_m-1):
        k = a[i].split()
        #print(k)
        list_a.extend(k)
    a=np.array(list_a,dtype=np.int32)
    #Read file contents as an array
    
    max_index = a.shape[0]
    mat_m = a[0]
    mat_n = a[1]
    maxincol = a[2]
    num = sum(a[4:4+mat_m])
    print(num)
    start = 4 + mat_m + mat_n
    k = 0
    position = np.zeros((mat_m,maxincol),dtype=np.int32)
    for i in range(mat_m):
        for j in range(maxincol):
            position[i,j]=a[start+k]-1
            # print(position[i,j])
            k = k+1

    ii = np.zeros(num,dtype=np.int32)#iiThe number of middle elements is the number of non-sparse elements
    jj = np.zeros(num,dtype=np.int32)#Note that the assignment of numpy adds a reference and does not create a new variable
    ss = np.zeros(num,dtype=np.int32)
    k = 0; 
    for i in range(mat_m):
        for j in range(a[4+i]):
            ii[k] = i
            jj[k] = position[i,j]
            ss[k] = 1
            k = k+1 
          
    #ii records the position of the non-zero element row, jj records the position of the non-zero element column, ss records the non-zero element value, if both are 1, it can be simplified    
    #H = sparse(ii,jj,ss,m,n) matlab sparse m*n
    H = sps.csr_matrix((ss, (ii, jj)), shape=(mat_m,mat_n),dtype= np.int32)
    # print(H.todense())

#Test    
# def main():
#     alist2sparse('./2.alist')
    
# if __name__=="__main__":
#     main()
