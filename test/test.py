import numpy as np
from ldpc import mod2

from bposd import bposd_decoder

hx=np.load("test/rhombic_hx.txt").astype(int)

m,n=hx.shape

print(m,n)


mod2.rank(hx)


# exit(22)


bposd_decoder(
    hx,
    error_rate=0.05,
    max_iter=20,
    bp_method="ms",
    ms_scaling_factor=0,
    osd_method="osd_cs",
    osd_order=6
)