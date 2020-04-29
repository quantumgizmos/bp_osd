//
// Created by joschka on 09/02/2020.
//

#ifndef MATCHINGBP_BP_OSD_H
#define MATCHINGBP_BP_OSD_H


class bp_osd {

    public:
        mod2sparse *H;

        int M,N,rank;

        double *log_prob_ratios;
        double channel_prob;

        int max_iter;
        int *iter;
        int *converge;

        int osd_order;
        int osd_method;

//        char *synd;
        char *bp_decoding;
        char *osd0_decoding;
        char *osdw_decoding;

        osd_struct *osd_data;



    bp_osd(mod2sparse *H,double channel_prob, int max_iter=0, double osd_order=10, int osd_method=0); //constructor

        char *bp_decode(char *synd);
        char *bp_osd_decode(char *synd);
        char *osd_post_process(double *, char *);
        void test();


};


#endif //MATCHINGBP_BP_OSD_H
