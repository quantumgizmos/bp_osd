#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "alloc.h"
#include "mod2sparse.h"
#include "logical_check.h"
//#include "mod2dense.h"
//#include "mod2convert.h"
#include "mod2sparse_extra.h"
#include "binary_char.h"



int check_logical_error_lx(mod2sparse* l,char *orig_error,char *decoding){

    int n = mod2sparse_cols(l);
    int k = mod2sparse_rows(l);
    char *residual=chk_alloc(n,sizeof(*residual));
    char *log_check=chk_alloc(n,sizeof(*log_check));

    int logical_error=0;
    bin_char_add(orig_error,decoding,residual,n);
    mod2sparse_mulvec(l,residual,log_check);
    if(!bin_char_is_zero(log_check,k)) logical_error=1;

    free(log_check);
    free(residual);

    return logical_error;

}

int check_logical_error_hz(mod2sparse* h,char *orig_error,char *decoding){

    int logical_error = 0;

    int n = mod2sparse_cols(h);
    int c = mod2sparse_rows(h);

    char *residual=chk_alloc(n,sizeof(*residual));

    bin_char_add(orig_error,decoding,residual,n);

    mod2sparse *ht = mod2sparse_allocate(n,c);
    mod2sparse_transpose(h,ht);

    mod2sparse *htf = mod2sparse_allocate(n,n);

    mod2sparse_merge_vec(ht, residual, n, htf);

    int rank_ht = mod2sparse_rank(h); //we have rank(A^T) = rank(A)

    int rank_htf = mod2sparse_rank(htf);

    if (rank_ht < rank_htf) {
        logical_error = 1;
    }

    free(residual);

    mod2sparse_free(ht);
    mod2sparse_free(htf);

    return logical_error;

}