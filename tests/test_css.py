import scipy
import numpy as np
import pytest
from bposd.css import css_code
import ldpc.codes


def test_css_code_steane():
    h = ldpc.codes.hamming_code(3)

    qcode = css_code(hx=h, hz=h, code_distance=3, name="Steane code")

    assert qcode.N == 7
    assert qcode.K == 1
    assert qcode.D == 3

    assert qcode.test()

    h = ldpc.codes.hamming_code(3).toarray()

    qcode = css_code(hx=h, hz=h, code_distance=3, name="Steane code")

    assert qcode.N == 7
    assert qcode.K == 1
    assert qcode.D == 3

    assert qcode.test()
