import numpy as np
from alist import save_alist

hamming_matrix=np.array([[1,0,0,1,1,0,1],
                         [0,1,0,1,0,1,1],
                         [0,0,1,1,1,1,0]])

save_alist("hamming_d_3.alist",hamming_matrix)
