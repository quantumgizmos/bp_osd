import numpy as np
from ldpc.codes import rep_code
from bposd.hgp import hgp
from css_decode_sim2 import css_decode_sim2

sc=hgp(rep_code(4))
print(sc.lx)

lk=css_decode_sim2(hx=sc.hx,hz=sc.hz,error_rate=0.05,target_runs=1000,xyz_error_bias=[np.inf,0,0])

# print(lk.output_dict())