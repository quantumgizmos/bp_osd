#cython: language_level=3, boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True
from libc.stdlib cimport malloc, calloc,free
from libc.math cimport log, tanh, pow, isnan, atanh
cimport numpy as np

from ldpc.mod2sparse cimport *
from ldpc.c_util cimport numpy2char, char2numpy, numpy2double, double2numpy
from ldpc.bp_decoder cimport bp_decoder

cdef extern from "binary_char.h":
    cdef void print_char_nonzero(char *val,int len)
    cdef int bin_char_equal(char *vec1, char *vec2, int len)
    cdef int bin_char_is_zero(char *vec1, int len)
    cdef void print_char(char *val, int len)
    cdef int bin_char_add(char *vec1, char *vec2, char *out_vec, int len)
    cdef char *decimal_to_binary_reverse(int n,int K)
    cdef int bin_char_weight(char *val,int len)
    cdef int hamming_difference(char *v1,char *v2,int len)

cdef extern from "sort.h":
    cdef void soft_decision_col_sort(double *soft_decisions,int *cols, int N)

cdef extern from "mod2sparse_extra.h":
    cdef void mod2sparse_print_terminal (mod2sparse *A)
    cdef int mod2sparse_rank(mod2sparse *A)
    
    cdef void LU_forward_backward_solve(
        mod2sparse *L,
        mod2sparse *U,
        int *rows,
        int *cols,
        char *z,
        char *x)

    cdef int mod2sparse_decomp_osd(
        mod2sparse *A,
        int R,
        mod2sparse *L,
        mod2sparse *U,
        int *rows,
        int *cols)