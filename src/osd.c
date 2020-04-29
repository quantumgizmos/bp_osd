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


// Factorial function:
// http://www.rosettacode.org/wiki/Factorial#C


// long int factorial(int n)
// 	{
// 		if(n==0) return 1;

// 	long i;
//     long int fac=1;
//     for(i=1;i<=n;i++)
//       {
//         fac=fac*i;
//       }
//     return fac;

// 	}


// long int factorial(int n) {
//     return n == 0 ? 1 : n * factorial(n - 1);
// }


unsigned long long
ncr(int n, int k) {
    if (k > n) {
        return 0;
    }
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}




//https://stackoverflow.com/a/36714204
//Modified slighty to sort indices in ascending order.

struct str
{
    double value;
    int index;
};

int cmp(const void *a, const void *b)
{
    struct str *a1 = (struct str *)a;
    struct str *a2 = (struct str *)b;
    if ((*a1).value > (*a2).value)
        return 1;
    else if ((*a1).value < (*a2).value)
        return -1;
    else
        return 0;
}


void soft_decision_col_sort(double *soft_decisions,int *cols, int N){
  struct str *objects;
    objects=chk_alloc (N, sizeof *objects);
    for (int i = 0; i < N; i++)
    {
        objects[i].value = soft_decisions[i];
        objects[i].index = i;
    }
    //sort objects array according to value using qsort
    qsort(objects, N, sizeof(objects[0]), cmp);
    for (int i = 0; i < N; i++)
        cols[i]=objects[i].index;
    
    free(objects);

}


//https://stackoverflow.com/a/36714204
//Integer sort version Modified slighty to sort indices in ascending order.

struct str_int
{
    int value;
    int index;
};

int cmp_int(const void *a, const void *b)
{
    struct str_int *a1 = (struct str_int *)a;
    struct str_int *a2 = (struct str_int *)b;
    if ((*a1).value > (*a2).value)
        return 1;
    else if ((*a1).value < (*a2).value)
        return -1;
    else
        return 0;
}


void col_sort_int(int *integers,int *cols, int N){
  struct str_int *objects;
    objects=chk_alloc (N, sizeof *objects);
    for (int i = 0; i < N; i++)
    {
        objects[i].value = integers[i];
        objects[i].index = i;
    }
    //sort objects array according to value using qsort
    qsort(objects, N, sizeof(objects[0]), cmp_int);
    for (int i = 0; i < N; i++)
        cols[i]=objects[i].index;
    
    free(objects);

}




void reorder_char(char *bit_char, char *reordered_bit_char, int *indices){
  int length=sizeof(bit_char);
  char *temp;
  temp=chk_alloc (sizeof(bit_char), sizeof *temp);
  for(int j=0; j<length; j++){
    temp[j]=bit_char[j];
  }

  for(int j=0; j<length; j++){
    reordered_bit_char[j]=temp[indices[j]];
  }

  free(temp);
}






osd_struct *create_osd_cs_struct(mod2sparse *H, int w)
    {
      int M,N;

      N=mod2sparse_cols(H);
      
      int osd_methods_count=2;
      
      
      //setup OSD-0
      osd_struct *osd_data;
      osd_data=malloc(osd_methods_count*sizeof(*osd_data));
      osd_data[0].osd_methods_count=osd_methods_count;
      osd_data[0].order=0;
      sprintf(osd_data[0].name,"OSD-%i",0);
      for(int i=0;i<osd_methods_count;i++) 
        {
          osd_data[i].decoding=chk_alloc(N,sizeof(char*));
          osd_data[i].success=0;
        }
      if(w==0) return osd_data;

      sprintf(osd_data[1].name,"OSD_CS_%i",w);
      osd_data[1].order=w;
      


      return osd_data;

    } 


osd_struct *create_osd_e_struct(mod2sparse *H, int w, int k)
    {
      int M,N;

      N=mod2sparse_cols(H);
      
      int osd_methods_count=2;
      
      
      //setup OSD-0
      osd_struct *osd_data;
      osd_data=malloc(osd_methods_count*sizeof(*osd_data));
      osd_data[0].osd_methods_count=osd_methods_count;
      osd_data[0].order=0;
      sprintf(osd_data[0].name,"OSD-%i",0);
      for(int i=0;i<osd_methods_count;i++) 
        {
          osd_data[i].decoding=chk_alloc(N,sizeof(char*));
          osd_data[i].success=0;
        }
      if(w==0) return osd_data;

      //Setup OSD-W strategy
      osd_search_strategy *osd_strategy;
      osd_strategy=malloc(pow(2,w)*sizeof(osd_search_strategy*));
      if (osd_strategy == NULL)
        {
          printf("Error allocating data for search strategy");
          exit(EXIT_FAILURE);
        }
      //Fill strategy with inputs for the encoding operator
      for(int i=0;i<pow(2,w);i++)
        {
          osd_strategy[i].encoding_op_input=decimal_to_binary_reverse(i,k);
          // print_char_nonzero(osd_strategy[i].encoding_op_input,k);
          // printf("\n");
        }
      sprintf(osd_data[1].name,"OSD-%i",w);
      osd_data[1].strategy=osd_strategy;
      osd_data[1].order=w;
      osd_data[1].strategy_size=pow(2,w);


      return osd_data;

    } 


// int in_list(int *list, int check, int len)
//   {

//     for(i=0;i<len;i++)
//       {
//         if(list[i]==check){return 1;}
//       }

//   return 0;

//   }



  int osd_e
( mod2sparse *A,	/* Input matrix, M by N */
  char *synd,
  double *log_prob_ratios,
  int R,		/* Size of sub-matrix to find LU decomposition of */
  osd_struct *osd_data
)
  {
    int M,N,k;

    M = mod2sparse_rows(A);
    N = mod2sparse_cols(A);
    k=N-R;

    int nnf;
    mod2sparse *L;
    mod2sparse *U;
    int *rows;
    int *cols;
    int *orig_cols;

    cols=chk_alloc(N,sizeof(*rows));
    orig_cols=chk_alloc(N,sizeof(*rows));
    rows=chk_alloc(M,sizeof(*cols));




    L=mod2sparse_allocate(M,R);
    U=mod2sparse_allocate(R,N);

    soft_decision_col_sort(log_prob_ratios,cols,N);


    for(int i=0;i<N;i++){
      orig_cols[i]=cols[i];
    }

    // for(int i=0;i<M;i++){
    //   printf("%i ",cols[i]);
    // }
    // printf("\n");


    int abandon_number=0;  	/* Number of columns to abandon at some point *//* When to abandon these columns */
    int abandon_when=0;
    mod2sparse_strategy strategy =Mod2sparse_first;/* Strategy to follow in picking rows/columns */

      nnf=mod2sparse_decomp_osd
    (A,	/* Input matrix, M by N */
    R,		/* Size of sub-matrix to find LU decomposition of */
    L,	/* Matrix in which L is stored, M by R */
    U,	/* Matrix in which U is stored, R by N */
    rows,		/* Array where row indexes are stored, M long */
    cols,		/* Array where column indexes are stored, N long */
    strategy, /* Strategy to follow in picking rows/columns */
    abandon_number,	/* Number of columns to abandon at some point */
    abandon_when	/* When to abandon these columns */
    );

  
    // for(int i=0;i<M;i++){
    //   printf("%i ",cols[i]);
    // }
    // printf("\n");

    int check;
    int counter=0;
    int in_pivot=1;
    for(int i=0;i<N;i++)
    {
      check=orig_cols[i];
      in_pivot=1;
      for(int i=0;i<M;i++)
      {
        if(cols[i]==check)
        {
          in_pivot=0;
          break;
        }
      }

        if(in_pivot==1){
          cols[counter+M]=check;
          counter=counter+1;
        }

    }

    // for(int i=M;i<N;i++){
    //   printf("%i ",cols[i]);
    // }
    // printf("\n.........\n");

    // exit(1);

    LU_solve_osd
        (L,
        U,
        rows,
        cols,
        synd,
        osd_data[0].decoding);

    // printf("%s",osd_data[0].decoding);
    if(osd_data[0].osd_methods_count==1)
      {    

        free(rows);
        free(cols);
        free(orig_cols);
        mod2sparse_free(L);
        mod2sparse_free(U);
        return 0;
      }



    mod2sparse *Ht;
    Ht=mod2sparse_allocate(M,k);
    int *Ht_cols;
    Ht_cols=chk_alloc(k,sizeof(*Ht_cols));
    for(int col_no=0;col_no<k;col_no++)
      {
        Ht_cols[col_no]=cols[col_no+R];
      }
    mod2sparse_copycols(A,Ht,Ht_cols);
    
    char *x;
    char *Htx;
    char *g;
    char *y;
    int encoding_op_count;
    int osd_min_weight;
    int osd_0_weight;
    int solution_weight;

    /*Calculate the weight of the OSD_0 solution. Higher order solutions are compared against this*/
    osd_0_weight=bin_char_weight(osd_data[0].decoding,N);

    y=chk_alloc(N,sizeof(*y));
    g=chk_alloc(M,sizeof(*g));
    Htx=chk_alloc(M,sizeof(*Htx));

    for(int method_no=1; method_no<osd_data[0].osd_methods_count;method_no++)
      {

        /*Set OSD_0 as the benchmark to beat*/
        osd_min_weight=osd_0_weight;
        for(int bit_no=0;bit_no<N;bit_no++) osd_data[method_no].decoding[bit_no]=osd_data[0].decoding[bit_no];

        /*Search through the encoding strings*/
        encoding_op_count=osd_data[method_no].strategy_size;
        for(int j=0; j<encoding_op_count;j++)

          {
            


            x=osd_data[method_no].strategy[j].encoding_op_input;

            // print_char_nonzero(osd_data[method_no].strategy[j].encoding_op_input,k);
            // // printf("Encoding op: ");
            // // print_char(x,k);
            // printf("\n");

           

            

            // char *test;
            // test=chk_alloc(k,sizeof(char*));
            // // mod2sparse_print_terminal(Ht);

            // print_char(test,k);

            mod2sparse_mulvec(Ht,x,Htx);
            // printf("\nHtx: ");
            // print_char(Htx,M);
            // printf("\n");
            bin_char_add(synd,Htx,g,M);

            LU_solve_osd(L,U,rows,cols,g,y);


            for(int col_no=0;col_no<k;col_no++)
                  {
                    y[Ht_cols[col_no]]=x[col_no];
                  }

            solution_weight=bin_char_weight(y,N);
            // print_char_nonzero(y,N);
            // printf("\n");
            // syndrome(A,y,g);
            // print_char_nonzero(g,M);
            // printf("\n");
            // printf("Solution weight: %i, OSD min weight: %i\n",solution_weight,osd_min_weight);

            if(solution_weight<osd_min_weight)
              {
                // printf("hello\n");
                osd_min_weight=solution_weight;
                for(int bit_no=0;bit_no<N;bit_no++) osd_data[method_no].decoding[bit_no]=y[bit_no];
                // print_char_nonzero(osd_data[method_no].decoding,N);
                // print_char_nonzero(x,k);
                // printf("\n");
              }
            // else if(solution_weight<osd_0_weight)
            // {
            //     print_char_nonzero(x,k);
            //     printf("\n");
            // }

            // printf("%i\n",solution_weight);

          }

          //  exit(1);






      }

    


    // mod2sparse_print_terminal(Ht);

    // printf("\nOSD-W Syndrome");
    // syndrome(A,y,g);
    // print_char_nonzero(g,M);
    // printf("\n");
    // printf("\nOSD-W Decode");
    // // syndrome(A,y,g);
    // print_char_nonzero(y,M);
    // printf("\n");


    free(Htx);
    free(Ht_cols);
    free(y);
    // free(x);
    free(g);
    mod2sparse_free(Ht);
    free(rows);
    free(cols);
    free(orig_cols);
    mod2sparse_free(L);
    mod2sparse_free(U);

    return 0;

  }


int osd_cs
( mod2sparse *A,	/* Input matrix, M by N */
  char *synd,
  double *log_prob_ratios,
  int R,		/* Size of sub-matrix to find LU decomposition of */
  osd_struct *osd_data
)
  {

    int M,N,k;
    int nnf;
    mod2sparse *L;
    mod2sparse *U;
    int *rows;
    int *orig_cols;
    int *cols;

    M = mod2sparse_rows(A);
    N = mod2sparse_cols(A);
    k=N-R;


    cols=chk_alloc(N,sizeof(*rows));
    rows=chk_alloc(M,sizeof(*cols));
    orig_cols=chk_alloc(N,sizeof(*rows));


    

    L=mod2sparse_allocate(M,R);
    U=mod2sparse_allocate(R,N);

    soft_decision_col_sort(log_prob_ratios,cols,N);

    // printf("\n");
    // for(int i=0;i<10;i++) printf("%lf ",log_prob_ratios[cols[i]]);


    for(int i=0;i<N;i++){
      orig_cols[i]=cols[i];
    }


    int abandon_number=0;  	/* Number of columns to abandon at some point *//* When to abandon these columns */
    int abandon_when=0;
    mod2sparse_strategy strategy =Mod2sparse_first;/* Strategy to follow in picking rows/columns */


    nnf=mod2sparse_decomp_osd
    (A,	/* Input matrix, M by N */
    R,		/* Size of sub-matrix to find LU decomposition of */
    L,	/* Matrix in which L is stored, M by R */
    U,	/* Matrix in which U is stored, R by N */
    rows,		/* Array where row indexes are stored, M long */
    cols,		/* Array where column indexes are stored, N long */
    strategy, /* Strategy to follow in picking rows/columns */
    abandon_number,	/* Number of columns to abandon at some point */
    abandon_when	/* When to abandon these columns */
    );

    int check;
    int counter=0;
    int in_pivot=0;
    int dep_col_count=0;
    for(int i=0;i<N;i++)
    {
      check=orig_cols[i];
      in_pivot=0;
      for(int j=0;j<M;j++)
      {
        if(cols[j]==check)
        {
          in_pivot=1;
          break;
        }
      }

        if(in_pivot==0){
          if(i<M){dep_col_count+=1;}
          cols[counter+M]=check;
          counter=counter+1;
        }

    } 



    LU_solve_osd
        (L,
        U,
        rows,
        cols,
        synd,
        osd_data[0].decoding);

    // printf("%s",osd_data[0].decoding);
    if(osd_data[0].osd_methods_count==1)
      {    

        free(rows);
        free(cols);
        free(orig_cols);
        mod2sparse_free(L);
        mod2sparse_free(U);
        return 0;
      }

    // exit(1);

    mod2sparse *Ht;
    Ht=mod2sparse_allocate(M,k);
    int *Ht_cols;
    Ht_cols=chk_alloc(k,sizeof(*Ht_cols));
    for(int col_no=0;col_no<k;col_no++)
      {
        Ht_cols[col_no]=cols[col_no+R];
      }
    mod2sparse_copycols(A,Ht,Ht_cols);
    
    char *test_x;
    char *Htx;
    char *g;
    char *y;
    int encoding_op_count;
    int osd_min_weight;
    int osd_0_weight;
    int solution_weight;
    int search_depth;
    int *osd_w1_weights;
    int *osd_w1_cols;

    int abandon_search;


    int *flipped_bits;
    // int *solution_weight_level;
    int level;
    int level_min_weight;


    /*Calculate the weight of the OSD_0 solution. Higher order solutions are compared against this*/
    osd_0_weight=bin_char_weight(osd_data[0].decoding,N);

    y=chk_alloc(N,sizeof(*y));
    g=chk_alloc(M,sizeof(*g));

    test_x=chk_alloc(k,sizeof(*test_x));
    osd_w1_weights=chk_alloc(k,sizeof(*osd_w1_weights));
    osd_w1_cols=chk_alloc(k,sizeof(*osd_w1_weights));
   

    Htx=chk_alloc(M,sizeof(*Htx));


    /*Set OSD_0 as the benchmark to beat*/
    osd_min_weight=osd_0_weight;
    for(int bit_no=0;bit_no<N;bit_no++) osd_data[1].decoding[bit_no]=osd_data[0].decoding[bit_no];

    /*Search through the encoding strings*/
    // search_depth=osd_data[1].order;
    search_depth=osd_0_weight;

    flipped_bits=chk_alloc(search_depth,sizeof(int*));
    

    for(level=0; level<search_depth;level++) flipped_bits[level]=-1;
    // solution_weight_level=chk_alloc(search_depth,sizeof(int*));

    // printf("Dep col count: %i\n", dep_col_count);

    for(int j=0;j<k;j++)
    {
      test_x[j]=1;

      mod2sparse_mulvec(Ht,test_x,Htx);
      bin_char_add(synd,Htx,g,M);
      LU_solve_osd(L,U,rows,cols,g,y);

                for(int col_no=0;col_no<k;col_no++)
                      {
                        y[Ht_cols[col_no]]=test_x[col_no];
                      }

                solution_weight=bin_char_weight(y,N);

            osd_w1_weights[j]=solution_weight;

            if(solution_weight<osd_min_weight)
                  {
                    osd_min_weight=solution_weight;
                    for(int bit_no=0;bit_no<N;bit_no++) osd_data[1].decoding[bit_no]=y[bit_no];
                // print_char_nonzero(test_x,k);
                // printf("\n");
                  }
      test_x[j]=0;
    }



    col_sort_int(osd_w1_weights,osd_w1_cols,k);

    // for(int bit_no=0;bit_no<k;bit_no++){
    //   printf("%i ",osd_w1_cols[bit_no]);
    //   }
    // printf("\n");

    // exit(0);



    for(int i=0;i<osd_data[1].order;i++)
    {
      for(int j=0;j<osd_data[1].order;j++){
      
      if(i==j){continue;}

      // test_x[i]=1;
      // test_x[j]=1;

      test_x[osd_w1_cols[i]]=1;
      test_x[osd_w1_cols[j]]=1;


      mod2sparse_mulvec(Ht,test_x,Htx);
      bin_char_add(synd,Htx,g,M);
      LU_solve_osd(L,U,rows,cols,g,y);

                for(int col_no=0;col_no<k;col_no++)
                      {
                        y[Ht_cols[col_no]]=test_x[col_no];
                      }

                solution_weight=bin_char_weight(y,N);

            if(solution_weight<osd_min_weight)
                  {
                    osd_min_weight=solution_weight;
                    for(int bit_no=0;bit_no<N;bit_no++) osd_data[1].decoding[bit_no]=y[bit_no];
                // print_char_nonzero(test_x,k);
                // printf("\n");
                  }
      test_x[osd_w1_cols[i]]=0;
      test_x[osd_w1_cols[j]]=0;

    // test_x[i]=0;
    // test_x[j]=0;

    }
    }


    free(flipped_bits);
    free(test_x);
    free(Htx);
    free(Ht_cols);
    free(y);
    free(g);
    mod2sparse_free(Ht);
    free(rows);
    free(cols);
    free(orig_cols);
    mod2sparse_free(L);
    mod2sparse_free(U);
    free(osd_w1_weights);
    free(osd_w1_cols);
   

    return 0;

  }


 


void LU_solve_osd
(mod2sparse *L,
mod2sparse *U,
int *rows,
int *cols,
char *z,
char *x)
{
    int N,R;
    char *forward_b;
    N=mod2sparse_cols(U);
    R=mod2sparse_cols(L);
    forward_b=chk_alloc(R,sizeof(*forward_b));

    for(int bit_no=0;bit_no<N;bit_no++) x[bit_no]=0;

    mod2sparse_forward_sub
    ( L,	/* Matrix that is lower triangular after reordering */
      rows,		/* Array of indexes (from 0) of rows for new order */
      z,		/* Vector on right of equation, also reordered */
      forward_b		/* Place to store solution */
    );
    mod2sparse_backward_sub
    ( U,	/* Matrix that is lower triangular after reordering */
      cols,		/* Array of indexes (from 0) of cols for new order */
      forward_b,		/* Vector on right of equation, also reordered */
      x		/* Place to store solution */
    );
    free(forward_b);
}





/* FIND AN LU DECOMPOSITION OF A SPARSE MATRIX. */

int mod2sparse_decomp_osd
( mod2sparse *A,	/* Input matrix, M by N */
  int R,		/* Size of sub-matrix to find LU decomposition of */
  mod2sparse *L,	/* Matrix in which L is stored, M by R */
  mod2sparse *U,	/* Matrix in which U is stored, R by N */
  int *rows,		/* Array where row indexes are stored, M long */
  int *cols,		/* Array where column indexes are stored, N long */
  mod2sparse_strategy strategy, /* Strategy to follow in picking rows/columns */
  int abandon_number,	/* Number of columns to abandon at some point */
  int abandon_when	/* When to abandon these columns */
)
{  
  int *rinv, *cinv, *acnt, *rcnt;
  mod2sparse *B;
  int M, N;

  mod2entry *e, *f, *fn, *e2;
  int i, j, k, cc, cc2, cc3, cr2, pr;
  int found, nnf;

  M = mod2sparse_rows(A);
  N = mod2sparse_cols(A);

  if (mod2sparse_cols(L)!=R || mod2sparse_rows(L)!=M
   || mod2sparse_cols(U)!=N || mod2sparse_rows(U)!=R)
  { fprintf (stderr,
      "mod2sparse_decomp: Matrices have incompatible dimensions\n");
    exit(1);
  }

  if (abandon_number>N-R)
  { fprintf(stderr,"Trying to abandon more columns than allowed\n");
    exit(1);
  }

  rinv = chk_alloc (M, sizeof *rinv);
  cinv = chk_alloc (N, sizeof *cinv);

  if (abandon_number>0)
  { acnt = chk_alloc (M+1, sizeof *acnt);
  }

  if (strategy==Mod2sparse_minprod)
  { rcnt = chk_alloc (M, sizeof *rcnt);
  }


  //these are the problematic functions!
  mod2sparse_clear(L);
  mod2sparse_clear(U);

  /* Copy A to B.  B will be modified, then discarded. */

  B = mod2sparse_allocate(M,N);
  mod2sparse_copy(A,B);

  /* Count 1s in rows of B, if using minprod strategy. */

  if (strategy==Mod2sparse_minprod)
  { for (i = 0; i<M; i++) 
    { rcnt[i] = mod2sparse_count_row(B,i);
    }
  }

  /* Set up initial row and column choices. */

  for (i = 0; i<M; i++) rows[i] = rinv[i] = i;
  for (j = 0; j<N; j++) {
    
    cinv[cols[j]]=j;
 
  }
  /* Find L and U one column at a time. */

  nnf = 0;

  for (i = 0; i<R; i++)
  { 
    /* Choose the next row and column of B. */

    switch (strategy)
    {
      case Mod2sparse_first: 
      { 
        found = 0;

        for (k = i; k<N; k++)
        { e = mod2sparse_first_in_col(B,cols[k]);
          while (!mod2sparse_at_end(e))
          { if (rinv[mod2sparse_row(e)]>=i)
            { found = 1;
              goto out_first;
            }
            e = mod2sparse_next_in_col(e);
          }
        }

      out_first:
        break;
      }

      case Mod2sparse_mincol:
      { 
        found = 0;

        for (j = i; j<N; j++)
        { cc2 = mod2sparse_count_col(B,cols[j]);
          if (!found || cc2<cc)
          { e2 = mod2sparse_first_in_col(B,cols[j]);
            while (!mod2sparse_at_end(e2))
            { if (rinv[mod2sparse_row(e2)]>=i)
              { found = 1;
                cc = cc2;
                e = e2;
                k = j;
                break;
              }
              e2 = mod2sparse_next_in_col(e2);
            }
          }
        }

        break;
      }

      case Mod2sparse_minprod:
      { 
        found = 0;

        for (j = i; j<N; j++)
        { cc2 = mod2sparse_count_col(B,cols[j]);
          e2 = mod2sparse_first_in_col(B,cols[j]);
          while (!mod2sparse_at_end(e2))
          { if (rinv[mod2sparse_row(e2)]>=i)
            { cr2 = rcnt[mod2sparse_row(e2)];
              if (!found || cc2==1 || (cc2-1)*(cr2-1)<pr)
              { found = 1;
                pr = cc2==1 ? 0 : (cc2-1)*(cr2-1);
                e = e2;
                k = j;
              }
            }
            e2 = mod2sparse_next_in_col(e2);
          }
        }

        break;
      }

      default:
      { fprintf(stderr,"mod2sparse_decomp: Unknown stategy\n");
        exit(1);
      }
    }

    if (!found) 
    { nnf += 1;
    }

 

    /* Update 'rows' and 'cols'.  Looks at 'k' and 'e' found above. */

    if (found)
    { 
      // if (cinv[mod2sparse_col(e)]!=k) abort();
      if (cinv[mod2sparse_col(e)]!=k){
        printf("\n e: %i, k: %i",mod2sparse_col(e),k);
        printf("\nError. Exiting.");
        exit(1);
      }


      cols[k] = cols[i];
      cols[i] = mod2sparse_col(e);

      cinv[cols[k]] = k;
      cinv[cols[i]] = i;

      k = rinv[mod2sparse_row(e)];

      if (k<i) abort();

      rows[k] = rows[i];
      rows[i] = mod2sparse_row(e);

      rinv[rows[k]] = k;
      rinv[rows[i]] = i;
    }



    /* Update L, U, and B. */

    f = mod2sparse_first_in_col(B,cols[i]); 

    while (!mod2sparse_at_end(f))
    { 
      fn = mod2sparse_next_in_col(f);
      k = mod2sparse_row(f);
      if (rinv[k]>i)
      { mod2sparse_add_row(B,k,B,mod2sparse_row(e));
        if (strategy==Mod2sparse_minprod) 
        { rcnt[k] = mod2sparse_count_row(B,k);
        }
        mod2sparse_insert(L,k,i);
      }
      else if (rinv[k]<i)
      { mod2sparse_insert(U,rinv[k],cols[i]);
      }
      else
      { mod2sparse_insert(L,k,i);
        mod2sparse_insert(U,i,cols[i]);
      }

      f = fn;
    }


    /* Get rid of all entries in the current column of B, just to save space. */

    for (;;)
    { f = mod2sparse_first_in_col(B,cols[i]);
      if (mod2sparse_at_end(f)) break;
      mod2sparse_delete(B,f);
    }

    /* Abandon columns of B with lots of entries if it's time for that. */

    if (abandon_number>0 && i==abandon_when)
    { 
      for (k = 0; k<M+1; k++) 
      { acnt[k] = 0;
      }
      for (j = 0; j<N; j++) 
      { k = mod2sparse_count_col(B,j);
        acnt[k] += 1;
      }

      cc = abandon_number;
      k = M;
      while (acnt[k]<cc)
      { cc -= acnt[k];
        k -= 1;
        if (k<0) abort();
      }

      cc2 = 0;
      for (j = 0; j<N; j++)
      { cc3 = mod2sparse_count_col(B,j);
        if (cc3>k || cc3==k && cc>0)
        { if (cc3==k) cc -= 1;
          for (;;)
          { f = mod2sparse_first_in_col(B,j);
            if (mod2sparse_at_end(f)) break;
            mod2sparse_delete(B,f);
          }
          cc2 += 1;
        }
      }

      if (cc2!=abandon_number) abort();

      if (strategy==Mod2sparse_minprod)
      { for (j = 0; j<M; j++) 
        { rcnt[j] = mod2sparse_count_row(B,j);
        }
      }
    }
  }

  /* Get rid of all entries in the rows of L past row R, after reordering. */

  for (i = R; i<M; i++)
  { for (;;)
    { f = mod2sparse_first_in_row(L,rows[i]);
      if (mod2sparse_at_end(f)) break;
      mod2sparse_delete(L,f);
    }
  }

  mod2sparse_free(B);
  free(rinv);
  free(cinv);
  if (strategy==Mod2sparse_minprod) free(rcnt);
  if (abandon_number>0) free(acnt);

  return nnf;
}
