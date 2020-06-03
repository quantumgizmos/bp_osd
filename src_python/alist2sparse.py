import numpy as np

def alist2sparse(fname):
    a = np.loadtxt(fname, delimiter='\n',dtype=str)
    alist_n = a.shape[0]
    list_a = []
    for i in range(alist_n-1):
        k = a[i].split()
        #print(k)
        list_a.extend(k)
    a=np.array(list_a,dtype=np.int32)
    #Read file contents as an array  
    max_index = a.shape[0]
    mat_n = a[0]
    mat_m = a[1]
    maxincol = a[2]
    # print(maxincol)
    maxinrow = a[3]
    num = sum(a[4:4+mat_n])#Number of non-zero elements
    index_col_num = a[4:4+mat_n]#Number of non-zero elements column by column
    # print(num)
    start = 4 + mat_m + mat_n
    k = 0
    H = np.zeros((mat_m,mat_n),dtype=np.int32)
    for i in range(mat_n):#The index of each column must be read in to rebuild the index
        for j in range(index_col_num[i]):
            if(k==(num)):
                break
            H[a[start+k]-1,i]=1
            k = k+1
            
    #because save_alist has `.T`  ,so return `.T`        
    return H.T 
    
# def main():
#     H = alist2sparse('./hamming_d_3.alist')
#     print(H)
# if __name__=="__main__":
#     main()
