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
  int A_rank,		/* Size of sub-matrix to find LU decomposition of */
  osd_struct *osd_data
)
  {
    int M,N,k;

    M = mod2sparse_rows(A);
    N = mod2sparse_cols(A);
    k= N - A_rank;

    int nnf;
    mod2sparse *L;
    mod2sparse *U;
    int *rows;
    int *cols;
    int *orig_cols;

    cols=chk_alloc(N,sizeof(*rows));
    orig_cols=chk_alloc(N,sizeof(*rows));
    rows=chk_alloc(M,sizeof(*cols));


    L=mod2sparse_allocate(M, A_rank);
    U=mod2sparse_allocate(A_rank, N);

    soft_decision_col_sort(log_prob_ratios,cols,N);


    for(int i=0;i<N;i++){
      orig_cols[i]=cols[i];
    }


    int abandon_number=0;  	/* Number of columns to abandon at some point *//* When to abandon these columns */
    int abandon_when=0;
    mod2sparse_strategy strategy =Mod2sparse_first;/* Strategy to follow in picking rows/columns */

    nnf=mod2sparse_decomp_osd
    (A,	/* Input matrix, M by N */
    A_rank,		/* Size of sub-matrix to find LU decomposition of */
    L,	/* Matrix in which L is stored, M by A_rank */
    U,	/* Matrix in which U is stored, A_rank by N */
    rows,		/* Array where row indexes are stored, M long */
    cols		/* Array where column indexes are stored, N long */
    );


    int check;
    int counter=0;
    int in_pivot;
    for(int i=0;i<N;i++)
    {
      check=orig_cols[i];
      in_pivot=0;
      for(int i=0;i<A_rank;i++)
      {
        if(cols[i]==check)
        {
          in_pivot=1;
          break;
        }
      }

        if(in_pivot==0){
          cols[counter+A_rank]=check;
          counter=counter+1;
        }

    }


      lu_forward_backward_solve
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
        Ht_cols[col_no]=cols[col_no + A_rank];
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

              lu_forward_backward_solve(L, U, rows, cols, g, y);


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
    cols
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


      lu_forward_backward_solve
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
        lu_forward_backward_solve(L, U, rows, cols, g, y);

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
          lu_forward_backward_solve(L, U, rows, cols, g, y);

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


 



