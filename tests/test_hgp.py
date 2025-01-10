import numpy as np
import scipy
import pytest
import ldpc

from bposd.hgp import hgp


def test_hgp_surface():
    h = ldpc.codes.rep_code(3)

    qcode = hgp(h1=h, h2=h, compute_distance=True)

    assert qcode.test()

    assert qcode.N == 13
    assert qcode.K == 1
    assert qcode.D == 3
