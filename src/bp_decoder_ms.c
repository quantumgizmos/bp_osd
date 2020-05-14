/* Belief propagation decoder */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "alloc.h"
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"
#include "rand.h"
#include "rcode.h"
#include "check.h"

#include "syndrome.h"
#include "bp_decoder_ms.h"
#include "binary_char.h"


double log_prob_ratio_initial;

void bp_setup_ms(mod2sparse *H,
    double error_prob,
    double *log_prob_ratios,
    char *decoding)
    {

  mod2entry *e;
  int N;
  int j;
  

  log_prob_ratio_initial=log((1-error_prob)/error_prob);
  N = mod2sparse_cols(H);

  for (j = 0; j<N; j++)
    { 
      for (e = mod2sparse_first_in_col(H,j);
            !mod2sparse_at_end(e);
            e = mod2sparse_next_in_col(e))
        { 
          e->pr = log_prob_ratio_initial;
          e->lr = 1;
        }
    // log_prob_ratios[j]=log_prob_ratio_initial;
    decoding[j]=0;
    }

}

void bp_update_ms(mod2sparse *H,
    char *synd,
    double *log_prob_ratios,
    int iter,
    char *decoding)
    {

  double pr, dR,nX,dX;
  mod2entry *e;
  int N, M;
  int i, j;
  double sgn;
  int mod2row_weight;
  int synd_sgn;
  double min1, min2;
  double alpha;


  M = mod2sparse_rows(H);
  N = mod2sparse_cols(H);

  char *temp_synd;
  char *temp_dec;

  //scaling factor
  iter++;
  alpha=1.0 - pow(2.0,  (-1*(double) iter)/1  );


  //Recompute likelihood ratios.
  for (i = 0; i<M; i++)
  { 
    
    mod2row_weight=0;
    if(synd[i]==0)
      { 
        mod2row_weight=0;
      }
    else if(synd[i]==1)
      { 
        mod2row_weight=1;
      }

    min1=1e308;
    min2=1e308;

    for (e = mod2sparse_first_in_row(H,i);
         !mod2sparse_at_end(e);
         e = mod2sparse_next_in_row(e))
      { 
        if (fabs(e->pr) < min1)
            {
                min2=min1;
                min1=fabs(e->pr);
            }
        else if(fabs(e->pr) < min2)
            {
                min2=fabs(e->pr);
            }
        
        if(e->pr <= 0) mod2row_weight++;
      }

    for (e = mod2sparse_last_in_row(H,i);
         !mod2sparse_at_end(e);
         e = mod2sparse_prev_in_row(e))
      {
        if(e->pr <= 0) sgn=sgn+ mod2row_weight;
        else sgn=mod2row_weight;

        sgn=pow(-1,sgn);
        if (fabs(e->pr)==min1) e->lr=alpha*sgn*min2;
        else e->lr=alpha*sgn*min1;

      }

  }

  //Recompute log-probability-ratios for the bits
  for (j = 0; j<N; j++)
    { 
      pr = log_prob_ratio_initial;
      for (e = mod2sparse_first_in_col(H,j);
          !mod2sparse_at_end(e);
          e = mod2sparse_next_in_col(e))
        {
          pr+=e->lr;
        }
      log_prob_ratios[j]=pr;

      decoding[j] = pr<=0;

      for (e = mod2sparse_last_in_col(H,j);
          !mod2sparse_at_end(e);
          e = mod2sparse_prev_in_col(e))
        {
          e->pr = pr-(e->lr);

        }
    }

}

int bp_decode_ms(
mod2sparse *H,
char *synd,
double error_prob,
int max_iter,
int *converge,
int *iter,
char *decoding,
double *log_prob_ratios)
  {

    int M,N;
    int it;
    M = mod2sparse_rows(H);
    N = mod2sparse_cols(H);
    int hamming_weight;
    char *candidate_synd;

    candidate_synd=chk_alloc(M,sizeof(*candidate_synd));

    bp_setup_ms(H,error_prob,log_prob_ratios,decoding);

    for(it=0; it<max_iter;it++)
      {
        bp_update_ms(H,synd, log_prob_ratios, it, decoding);
        syndrome(H,decoding,candidate_synd);

        //check hamming weight between candidate syndrome and target syndrome
        hamming_weight=0;
        for(int check_no=0;check_no<M;check_no++)
          {
            hamming_weight+=synd[check_no]^candidate_synd[check_no];
          }

        //exit if candidate syndrome has converged to the target solution
        if(hamming_weight==0)
          {
            free(candidate_synd);
            *iter=it;
            *converge=1;
            return 1;
          }


      }
 
    free(candidate_synd);
    *iter=it;
    *converge=0;
    return 0;
  }



int bp_decode_ms_min_synd(
        mod2sparse *H,
        char *synd,
        double error_prob,
        int max_iter,
        int *converge,
        int *iter,
        char *decoding,
        double *log_prob_ratios)
{

    int M,N;
    int it;
    M = mod2sparse_rows(H);
    N = mod2sparse_cols(H);
    int hamming_weight;
    char *candidate_synd;
    double *log_prob_ratios_save;

    int decreased=0;

    int synd_weight=bin_char_weight(synd,M);
    int hamming_weight_min=synd_weight;

    log_prob_ratios_save=chk_alloc(N,sizeof(*log_prob_ratios_save));
    candidate_synd=chk_alloc(M,sizeof(*candidate_synd));

    bp_setup_ms(H,error_prob,log_prob_ratios,decoding);

    for(it=0; it<max_iter;it++)
    {
        bp_update_ms(H,synd, log_prob_ratios, it, decoding);
        syndrome(H,decoding,candidate_synd);

        hamming_weight=0;
        for(int check_no=0;check_no<M;check_no++)
        {
            hamming_weight+=synd[check_no]^candidate_synd[check_no];
        }

        if(hamming_weight==0)
        {
            free(candidate_synd);
            free(log_prob_ratios_save);
            *iter=it;
            *converge=1;
            return 1;
        }


        if(hamming_weight<hamming_weight_min){
            hamming_weight_min=hamming_weight;
            memcpy(log_prob_ratios_save,log_prob_ratios,N*sizeof(*log_prob_ratios));
            decreased=1;
        }

    }

    if(decreased) memcpy(log_prob_ratios,log_prob_ratios_save,N*sizeof(*log_prob_ratios));


    free(candidate_synd);
    free(log_prob_ratios_save);
    *iter=it;
    *converge=0;
    return 0;
}