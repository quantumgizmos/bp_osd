import numpy as np
from bposd.hgp import hgp

classical_seed_codes = ["mkmn_16_4_6.txt", "mkmn_20_5_8.txt", "mkmn_24_6_10.txt"]

for code in classical_seed_codes:
    seed_code = np.loadtxt(f"examples/codes/classical_seed_codes/{code}").astype(int)
    # print(seed_code)
    qcode = hgp(seed_code, compute_distance=True)

    qcode.canonical_logicals()

    qcode.test()

    # print(qcode.code_params)

    np.savetxt(f"examples/codes/hgp_codes/hgp_{qcode.code_params}_hx.txt", qcode.hx)
    np.savetxt(f"examples/codes/hgp_codes/hgp_{qcode.code_params}_hz.txt", qcode.hz)
    np.savetxt(f"examples/codes/hgp_codes/hgp_{qcode.code_params}_lx.txt", qcode.lx)
    np.savetxt(f"examples/codes/hgp_codes/hgp_{qcode.code_params}_lz.txt", qcode.lz)
