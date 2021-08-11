import numpy as np
from ldpc.codes import hamming_code,rep_code,ring_code


print(hamming_code(3))

import doctest

a=doctest.testmod(name="ldpc.hamming_code",verbose=True)
print(a)

print(ring_code(5))