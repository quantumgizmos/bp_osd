import ldpc
import ldpc.codes
import ldpc.code_util

h = ldpc.codes.rep_code(3).T

print(ldpc.code_util.compute_exact_code_distance(h))