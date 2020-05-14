void bp_setup_ms(
        mod2sparse *,
        double,
        double *,
        char *);
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












