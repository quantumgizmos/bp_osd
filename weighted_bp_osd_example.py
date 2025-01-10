"""
Example demonstrating how to run bp+osd with a non-uniform initial error channel.
Before running, install `bposd' and `ldpc`
`pip install -U ldpc bposd`
"""

import numpy as np
from bposd.hgp import hgp
from bposd import bposd_decoder
from ldpc.codes import rep_code

"""'
Distance 2 surface code
"""
# build surface code via hypergraph product
h = rep_code(2)  # repetition code seed
surface_code = hgp(h, h)  # build surface code
hz = surface_code.hz  # z stabilisers parity check matrix
print("Distance-2 surface code hz parity check matrix")
print(hz)

"""
Example error
"""
m, n = hz.shape
error = np.zeros(n).astype(int)
error[[2]] = 1
syndrome = hz @ error % 2

print("Error", error)
print("Syndrome", syndrome)

"""
Decoding with uniform initial error channel
"""

rho = 0.1 * np.ones(n)
# rho[0]=0.0
bpd = bposd_decoder(
    hz,
    channel_probs=rho,
    bp_method="product_sum",
    osd_method="osd_e",
    osd_order=3,
    max_iter=5,
)

decoding = bpd.decode(syndrome)
print("Uniform error channel", rho)
print("Decoding (uniform initial error channel)", decoding)

"""
Decoding with non-uniform initial error channel
"""

rho[2] = (
    0.05  # here I'm setting the probability on qubit 2 to be lower than the other qubits
)
bpd = bposd_decoder(
    hz,
    channel_probs=rho,
    bp_method="product_sum",
    osd_method="osd_e",
    osd_order=3,
    max_iter=5,
)

decoding = bpd.decode(syndrome)
print("Non-uniform error channel", rho)
print("Decoding (non-uniform initial error channel)", decoding)
