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
#include "bp_decoder_ps.h"
#include "binary_char.h"


double prob_ratio_intial;

void bp_setup_ps(mod2sparse *H,
    double error_prob,
    double *log_prob_ratios,
    char *decoding)
    {

  mod2entry *e;
  int N;
  int j;
  

  prob_ratio_intial=error_prob/(1-error_prob);
  N = mod2sparse_cols(H);

  for (j = 0; j<N; j++)
    { 
      for (e = mod2sparse_first_in_col(H,j);
            !mod2sparse_at_end(e);
            e = mod2sparse_next_in_col(e))
        { 
          e->pr = prob_ratio_intial;
          e->lr = 1;
        }
    decoding[j]=0;
    }

}

void bp_update_ps(mod2sparse *H,
    char *synd,
    double *log_prob_ratios,
    int iter,
    char *decoding)
    {
        double pr, dl, t;
        mod2entry *e;
        int N, M;
        int i, j;

        M = mod2sparse_rows(H);
        N = mod2sparse_cols(H);

        /* Recompute likelihood ratios. */

        for (i = 0; i<M; i++)
        {
            if(synd[i]==0) dl = 1;
            if(synd[i]==1) dl=-1;
            for (e = mod2sparse_first_in_row(H,i);
                 !mod2sparse_at_end(e);
                 e = mod2sparse_next_in_row(e))
            { e->lr = dl;
                dl *= 2/(1+e->pr) - 1;
            }
            dl = 1;
            for (e = mod2sparse_last_in_row(H,i);
                 !mod2sparse_at_end(e);
                 e = mod2sparse_prev_in_row(e))
            { t = e->lr * dl;
                e->lr = (1-t)/(1+t);
                dl *= 2/(1+e->pr) - 1;
            }
        }

        /* Recompute probability ratios.  Also find the next guess based on the
           individually most likely values. */

        for (j = 0; j<N; j++)
        { pr = prob_ratio_intial;
            for (e = mod2sparse_first_in_col(H,j);
                 !mod2sparse_at_end(e);
                 e = mod2sparse_next_in_col(e))
            { e->pr = pr;
                pr *= e->lr;
            }
            if (isnan(pr))
            { pr = 1;
            }
            log_prob_ratios[j]=log(1/pr);
            // if (bprb) bprb[j] = 1 - 1/(1+pr);
            decoding[j] = pr>=1;
            pr = 1;
            for (e = mod2sparse_last_in_col(H,j);
                 !mod2sparse_at_end(e);
                 e = mod2sparse_prev_in_col(e))
            { e->pr *= pr;
                if (isnan(e->pr))
                { e->pr = 1;
                }
                pr *= e->lr;
            }
        }

}

int bp_decode_ps(
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

    bp_setup_ps(H,error_prob,log_prob_ratios,decoding);

    for(it=0; it<max_iter;it++)
      {
        bp_update_ps(H,synd, log_prob_ratios, it, decoding);
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



int bp_decode_ps_min_synd(
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

    bp_setup_ps(H,error_prob,log_prob_ratios,decoding);

    for(it=0; it<max_iter;it++)
    {
        bp_update_ps(H,synd, log_prob_ratios, it, decoding);
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