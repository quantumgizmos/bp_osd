//
// Created by joschka on 09/02/2020.
//



#include <stdio.h>
#include <vector>
#include<iostream>
#include <cstring>

//C libraries
extern "C" {
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"

#include "load_alist.h"
#include "syndrome.h"
#include "binary_char.h"
#include "bp_decoder_ms.h"
#include "osd.h"
#include "mod2sparse_extra.h"
}

#include "bp_osd.h"

using namespace std;


bp_osd::bp_osd(mod2sparse *H, double channel_prob, int max_iter, double osd_order,int osd_method){

    this->H=H;
    this->N=mod2sparse_cols(H);
    this->M=mod2sparse_rows(H);

    this->rank=mod2sparse_rank(H);

    this->channel_prob=channel_prob;

    if(max_iter==0){this->max_iter=N;}
    else if(max_iter>N){this->max_iter=N;}
    else{this->max_iter=max_iter;}

    this->osd_order=osd_order;
    this->osd_method=osd_method;
    this->bp_decoding=new char[N]();
    this->osd0_decoding=new char[N]();
    this->osdw_decoding=new char[N]();
    this->log_prob_ratios=new double[N]();
    this->converge=new int[1]();
    this->iter=new int[1]();

if(osd_method==0) {
        this->osd_data=create_osd_e_struct(H,osd_order,N-rank);
    }
    else if(osd_method==1){
        this->osd_data=create_osd_cs_struct(H,osd_order);
    }
    else{
        cout<<"ERROR. Function <bp_osd::constructor>. OSD Method not valid"<<endl;
    }

};

char *bp_osd::bp_decode(char *synd) {

    bp_decode_ms_min_synd(
            H,
            synd,
            channel_prob,
            max_iter,
            converge,
            iter,
            bp_decoding,
            log_prob_ratios);





    return bp_decoding;



}


char *bp_osd::osd_post_process(double *soft_decisions, char *synd){

    if(osd_method==0) osd_e(H,synd,soft_decisions,rank,osd_data);
    else if(osd_method==1) osd_cs(H,synd,soft_decisions,rank,osd_data);
    else {
        cout << "ERROR. <Function bp_osd::osd_post_process>. Invalid OSD method!" << endl;
        exit(22);
    }
    osd0_decoding=osd_data[0].decoding;
    osdw_decoding=osd_data[1].decoding;

    return osdw_decoding;

}


char *bp_osd::bp_osd_decode(char *synd) {

    bp_decoding=bp_decode(synd);


    if(*converge){
        osd0_decoding=bp_decoding;
        osdw_decoding=bp_decoding;
        return bp_decoding;
    }

    if(osd_method==0) osd_e(H,synd,log_prob_ratios,rank,osd_data);
    else if(osd_method==1) osd_cs(H,synd,log_prob_ratios,rank,osd_data);
    else {
        cout << "ERROR. <Function bp_osd::bp_osd_decode>. Invalid OSD method!" << endl;
        exit(22);
    }

    osd0_decoding=osd_data[0].decoding;
    osdw_decoding=osd_data[1].decoding;

    return osdw_decoding;

}

void bp_osd::test(){
    printf("\nTest");
    for(int i=0; i<N;i++){printf("%i",bp_decoding[i]);}

}

