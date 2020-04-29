#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "alloc.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "mod2sparse_extra.h"


void mod2sparse_print_terminal(mod2sparse *A){
  int M,N;
  mod2dense *Adense;
  M=mod2sparse_rows(A);
  N=mod2sparse_cols(A);

  Adense=mod2dense_allocate(M,N);

  mod2sparse_to_dense(A,Adense);
for (int i = 0; i<mod2dense_rows(Adense); i++)
  { for (int j = 0; j<mod2dense_cols(Adense); j++)
    { printf("%d",mod2dense_get(Adense,i,j));
    }
    printf("\n");
  }

  mod2dense_free(Adense);

}

int mod2sparse_rank(mod2sparse *A){
    int M,N;

    M=mod2sparse_rows(A);
    N=mod2sparse_cols(A);

    int nnf,rank;
    mod2sparse *L;
    mod2sparse *U;
    int *rows;
    int *cols;

    cols=chk_alloc(N,sizeof(*rows));
    rows=chk_alloc(M,sizeof(*cols));

    L=mod2sparse_allocate(M,M);
    U=mod2sparse_allocate(M,N);

    int abandon_number=0;  	/* Number of columns to abandon at some point *//* When to abandon these columns */
    int abandon_when=0;
    mod2sparse_strategy strategy =Mod2sparse_first;/* Strategy to follow in picking rows/columns */

    nnf=mod2sparse_decomp
            (A,	/* Input matrix, M by N */
             M,		/* Size of sub-matrix to find LU decomposition of */
             L,	/* Matrix in which L is stored, M by R */
             U,	/* Matrix in which U is stored, R by N */
             rows,		/* Array where row indexes are stored, M long */
             cols,		/* Array where column indexes are stored, N long */
             strategy, /* Strategy to follow in picking rows/columns */
             abandon_number,	/* Number of columns to abandon at some point */
             abandon_when	/* When to abandon these columns */
            );


    free(rows);
    free(cols);
    free(L);
    free(U);

    rank=M-nnf;

    return rank;


}