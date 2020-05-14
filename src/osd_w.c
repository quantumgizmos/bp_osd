/* Ordered statistics decoding */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "alloc.h"
#include "mod2sparse.h"
#include "mod2sparse_extra.h"
#include "binary_char.h"
#include "syndrome.h"
#include "sort.h"
#include "osd_0.h"
#include "osd_w.h"

int osd_w(
        mod2sparse *A, //parity check matrix
        char *synd, //syndrome
        char *osd0_decoding, //char string for OSD0 decoding
        char *osdw_decoding, //char string for higher order decoding
        double *log_prob_ratios, //the soft-decisions (from BP)
        int osd_order, //OSD order
        int A_rank, //rank of parity check matrix
        char **encoding_operator_inputs, //the list of higher order terms to search over
        int encoding_input_count //the number of higher order terms to search over
        ){

    int M = mod2sparse_rows(A);
    int N = mod2sparse_cols(A);

    int k= N - A_rank;
    mod2sparse *L=mod2sparse_allocate(M, A_rank);
    mod2sparse *U=mod2sparse_allocate(A_rank, N);
    mod2sparse *Ht=mod2sparse_allocate(M,k);
    int *Ht_cols=chk_alloc(k,sizeof(*Ht_cols));
    int *cols=chk_alloc(N,sizeof(*cols));
    int *rows=chk_alloc(M,sizeof(*rows));

    char *x;
    char *y=chk_alloc(N,sizeof(*y));
    char *g=chk_alloc(M,sizeof(*g));
    char *Htx=chk_alloc(M,sizeof(*Htx));

    //first we find the OSD-0 solutions
    osd_0_solve(
            A,
            L,
            U,
            cols,
            rows,
            synd,
            osd0_decoding,
            log_prob_ratios,
            A_rank,
            1);

    //create the H_T matrix
    for(int col_no=0;col_no<k;col_no++){
        Ht_cols[col_no]=cols[col_no + A_rank];
    }
    mod2sparse_copycols(A,Ht,Ht_cols);


    /*Calculate the weight of the OSD_0 solution. Higher order solutions are compared against this*/
    int osd_0_weight=bin_char_weight(osd0_decoding,N);

    /*Set OSD_0 as the benchmark to beat*/
    int osd_min_weight=osd_0_weight;
    for(int bit_no=0;bit_no<N;bit_no++) osdw_decoding[bit_no]=osd0_decoding[bit_no]; //osdw_decoding is set to osd0_decoding initially

    /*Search through the encoding strings*/
    int solution_weight;
    for(long unsigned int j=0; j<encoding_input_count;j++){

        x=encoding_operator_inputs[j];
        mod2sparse_mulvec(Ht,x,Htx);
        bin_char_add(synd,Htx,g,M);
        LU_forward_backward_solve(L, U, rows, cols, g, y);

        for(int col_no=0;col_no<k;col_no++){
                y[Ht_cols[col_no]]=x[col_no];
            }

        solution_weight=bin_char_weight(y,N);

        if(solution_weight<osd_min_weight){
            osd_min_weight=solution_weight;
            for(int bit_no=0;bit_no<N;bit_no++) osdw_decoding[bit_no]=y[bit_no];
        }

    }


    free(Htx);
    free(Ht_cols);
    free(y);
    free(g);
    mod2sparse_free(Ht);
    free(rows);
    free(cols);
    mod2sparse_free(L);
    mod2sparse_free(U);

    return 0;

}




