import numpy as np
from bposd.stab import stab_code


def test_stab_five_qubit_code():
    h = np.array(
        [
            [1, 0, 1, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 1, 1, 0, 1, 0, 0, 1, 1],
            [0, 1, 1, 1, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        ]
    )

    qcode = stab_code(h[:, :5], h[:, 5:])

    assert qcode.test()
    assert qcode.N == 5
    assert qcode.K == 1

    assert qcode.compute_code_distance() == 3
