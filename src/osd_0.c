/* Ordered statistics decoding */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "alloc.h"
#include "mod2sparse.h"
#include "osd.h"
#include "mod2sparse_extra.h"
#include "binary_char.h"
#include "syndrome.h"
#include "sort.h"
#include "osd_0.h"

//This function is a wrapper for osd_0_solve (see below)
int osd_0(
        mod2sparse *A,
        char *synd,
        char *osd0_decoding,
        double *log_prob_ratios,
        int A_rank){


    int M = mod2sparse_rows(A);
    int N = mod2sparse_cols(A);

    int k= N - A_rank;
    mod2sparse *L=mod2sparse_allocate(M, A_rank);
    mod2sparse *U=mod2sparse_allocate(A_rank, N);
    int *cols=chk_alloc(N,sizeof(*cols));
    int *rows=chk_alloc(M,sizeof(*rows));


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
            0);


    mod2sparse_free(L);
    mod2sparse_free(U);
    free(rows);
    free(cols);


  }


  int osd_0_solve(
          mod2sparse *A, //parity check matrix
          mod2sparse *L,
          mod2sparse *U,
          int *cols,
          int *rows,
          char *synd,
          char *osd0_decoding,
          double *log_prob_ratios,
          int A_rank,
          int post_processing) {

      int M = mod2sparse_rows(A);
      int N = mod2sparse_cols(A);
      int k = N - A_rank;

      //sort columns on the basis of soft decisions
      soft_decision_col_sort(log_prob_ratios, cols, N);

      //save the original sorted column list
      int *orig_cols = chk_alloc(N, sizeof(*rows));
      for (int i = 0; i < N; i++) {
          orig_cols[i] = cols[i];
      }

      mod2sparse_decomp_osd
              (A,    /* Input matrix, M by N */
               A_rank,        /* Size of sub-matrix to find LU decomposition of */
               L,    /* Matrix in which L is stored, M by A_rank */
               U,    /* Matrix in which U is stored, A_rank by N */
               rows,        /* Array where row indexes are stored, M long */
               cols        /* Array where column indexes are stored, N long */
              );

      //this next step is necessary to ensure the bits outside of the info set are ordered according to the soft decisions. The decomp function scrambles them!
      if(post_processing==1) {
          //sort the bits outside the info set
          int check;
          int counter = 0;
          int in_pivot;
          for (int i = 0; i < N; i++) {
              check = orig_cols[i];
              in_pivot = 0;
              for (int i = 0; i < A_rank; i++) {
                  if (cols[i] == check) {
                      in_pivot = 1;
                      break;
                  }
              }

              if (in_pivot == 0) {
                  cols[counter + A_rank] = check;
                  counter = counter + 1;
              }

          }

      }

      //solve the syndrome equation using forwards/backwards substitution
      LU_forward_backward_solve(
              L,
              U,
              rows,
              cols,
              synd,
              osd0_decoding);


      //cleanup
      free(orig_cols);

  }

