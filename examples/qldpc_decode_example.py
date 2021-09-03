import numpy as np
from bposd.hgp import hgp
from css_decode_sim import css_decode_sim

h=np.loadtxt("examples/mkmn_16_4_6.txt").astype(int)
qcode=hgp(h) # construct quantum LDPC code using the symmetric hypergraph product

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

lk = css_decode_sim(hx=qcode.hx, hz=qcode.hz, **osd_options)


# import numpy as np
# from css_decode_sim import css_decode_sim as css_decode_sim
# from bposd.hgp import hgp_single

# #load the ldpc code from the file
# h=np.loadtxt("examples/mkmn_16_4_6.txt").astype(int)

# hgp=hgp_single(h) # construct quantum LDPC code using the symmetric hypergraph product

# output_dict={"code_type": "hgp_mkmn_16_4_6"}
# output_dict=css_decode_sim(
#     hx=hgp.hx, # CSS hx matrix
#     hz=hgp.hz, # CSS hz matrix
#     error_rate=0.05, # the physical error rate
#     # xyz_error_bias=[np.inf,0,0],
#     max_iter=0, # If `0` iterations depth is set to block length n
#     target_runs=1000, # the number of runs to simulate
#     seed=42, # random seed
#     bp_method="ms", # Options: 1) 'product_sum'; 2) 'minimum_sum'
#     ms_scaling_factor=0,
#     osd_method="osd_cs", # Options: 1) 'exhaustive'; 2) combination sweep
#     osd_order=40, # The OSD order
#     noise_type='z', # Options: 1) 'z'; 2) 'x'
#     output_file="hpg_bp_osd_decode_sim_output.json", #output directory
#     check_code=True
# )