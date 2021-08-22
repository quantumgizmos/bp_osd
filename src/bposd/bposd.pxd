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

cdef class bposd_decoder:
    cdef int MEM_ALLOCATED
    cdef mod2sparse* H
    cdef int m, n
    cdef char* error
    cdef char* synd
    cdef char* bp_decoding_synd
    cdef char* bp_decoding
    cdef char* decoding
    cdef char* osd0_decoding
    cdef char* osdw_decoding
    cdef int iter
    cdef int converge
    cdef double* channel_probs
    cdef double* log_prob_ratios
    cdef char** osdw_encoding_inputs
    cdef double error_rate
    cdef int max_iter
    cdef int bp_method
    cdef long int encoding_input_count
    cdef int osd_order
    cdef int osd_method
    cdef double ms_scaling_factor
    cdef int rank
    cdef int k
    cdef int i, j

    cdef bp_decoder bpd

    cdef int* rows
    cdef int* cols
    cdef int* orig_cols
    cdef int* Ht_cols
    cdef char* y
    cdef char* g
    cdef char* Htx

    cdef void osd_e_setup(self)
    cdef void osd_cs_setup(self)
    cdef int osd(self)
    cpdef np.ndarray[np.int_t, ndim=1] decode(self, np.ndarray[np.int_t, ndim=1] syndrome)
    cdef char* decode_cy(self, char* syndrome)

