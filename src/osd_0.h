int osd_0(mod2sparse *A,char *synd, char *osd0_decoding,double *log_prob_ratios,int A_rank);

int osd_0_solve(
        mod2sparse *A,
        mod2sparse *L,
        mod2sparse *U,
        int *cols,
        int *rows,
        char *synd,
        char *osd0_decoding,
        double *log_prob_ratios,
        int A_rank,
        int post_processing);
