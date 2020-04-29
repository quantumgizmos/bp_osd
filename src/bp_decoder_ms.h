/* DEC.H - Interface to decoding procedures. */

/* Copyright (c) 1995-2012 by Radford M. Neal.
 *
 * Permission is granted for anyone to copy, use, modify, and distribute
 * these programs and accompanying documents for any purpose, provided
 * this copyright notice is retained and prominently displayed, and note
 * is made of any changes made to these programs.  These programs and
 * documents are distributed without any warranty, express or implied.
 * As the programs were written for research purposes only, they have not
 * been tested to the degree that would be advisable in any important
 * application.  All use of these programs is entirely at the user's own
 * risk.
 */


/* DECODING METHOD, ITS PARAMETERS, AND OTHER VARIABLES.  The global variables 
   declared here are located in dec.c. */

typedef struct
{ 
   int converge; 
   int iter;
} bp_out;

//double log_prob_ratio_initial;
void bp_setup_ms(mod2sparse *,double, double *,char *);
void bp_update_ms(mod2sparse*, char *, double *, int , char *);

int bp_decode_ms(
        mod2sparse *H,
        char *synd,
        double error_prob,
        int max_iter,
        int *converge,
        int *iter,
        char *decoding,
        double *log_prob_ratios);

int bp_decode_ms_min_synd(
        mod2sparse *H,
        char *synd,
        double error_prob,
        int max_iter,
        int *converge,
        int *iter,
        char *decoding,
        double *log_prob_ratios);












