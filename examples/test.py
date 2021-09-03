import numpy as np
from ldpc.codes import rep_code
from bposd.hgp import hgp
from css_decode_sim2 import css_decode_sim2


osd_options={
'error_rate': 0.05,
'target_runs': 1000,
'xyz_error_bias': [0, 0, 1],
'output_file': 'test.json',
'bp_method': "ms",
'ms_scaling_factor': 0,
'osd_method': "osd_cs",
'osd_order': 40,
'channel_update': None,
'seed': 42,
'max_iter': 0,
'output_file': "test.json"
}

qcode=hgp(rep_code(3))

h=np.loadtxt("examples/mkmn_16_4_6.txt").astype(int)
qcode=hgp(h) # construct quantum LDPC code using the symmetric hypergraph product
lk = css_decode_sim2(hx=qcode.hx, hz=qcode.hz, **osd_options)

# from bposd import bposd_decoder

# bpd=bposd_decoder(qcode.hx,error_rate=0.05)

# print(bpd.channel_probs)

# probs=np.ones(qcode.N)*0.042

# bpd.update_channel_probs(probs)

# print(bpd.channel_probs)

