//
// Created by joschka on 09/02/2020.
//



#include <stdio.h>
#include <vector>
#include<iostream>
#include <cstring>
#include <assert.h>
#include <math.h>

//C include
extern "C" {
#include "mod2sparse.h"
#include "mod2dense.h"
#include "mod2convert.h"

#include "load_alist.h"
#include "syndrome.h"
#include "binary_char.h"
#include "bp_decoder_ms.h"
#include "mod2sparse_extra.h"
#include "osd_0.h"
#include "osd_w.h"
#include "osd_g.h"
#include "sort.h"
}

#include "bp_osd.h"

using namespace std;


bp_osd::bp_osd(mod2sparse *H, double channel_prob, int max_iter, double osd_order,int osd_method){

    this->H=H;
    this->N=mod2sparse_cols(H);
    this->M=mod2sparse_rows(H);

    this->rank=mod2sparse_rank(H);
    this->K=N-rank;

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
        assert(osd_order<=N-rank);
        this->encoding_input_count=pow(2,osd_order);
        this->osd_w_encoding_inputs = new char*[encoding_input_count];
        for(int i = 0; i < encoding_input_count; ++i)
            osd_w_encoding_inputs[i] = decimal_to_binary_reverse(i, N - rank);
        }
    else if(osd_method==1){
        assert(osd_order<=K);
        int total_count=0;
        int w2_count=ncr(osd_order,2);
        this->encoding_input_count = K+w2_count;

        this->osd_w_encoding_inputs = new char*[encoding_input_count];
        for(int i=0;i<K;i++){
            osd_w_encoding_inputs[total_count]=new char[K]();
            osd_w_encoding_inputs[total_count][i]=1;
            total_count++;
        }
        for(int i=0;i<osd_order;i++){
            for(int j=0;j<osd_order;j++){
                if(i<j){
                    osd_w_encoding_inputs[total_count]=new char[K]();
                    osd_w_encoding_inputs[total_count][i]=1;
                    osd_w_encoding_inputs[total_count][j]=1;
                    total_count++;
                }
            }
        }
        assert(total_count==this->encoding_input_count);

    }
    else if(osd_method==2) 1;
    else if(osd_method==3){
        assert(osd_order<=K);
        int w2_count=ncr(osd_order,2);
        this->encoding_input_count = K+w2_count;
    }

    else{
        cout<<"ERROR. Function <bp_osd::constructor>. OSD Method "<<osd_method<< " not valid"<<endl;
        exit(22);
    }

};

char *bp_osd::bp_decode(char *synd) {

    bp_decode_ms(
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

    if((osd_method==0)||(osd_method==1))
        osd_w(
                H,
                synd,
                osd0_decoding,
                osdw_decoding,
                soft_decisions,
                osd_order, rank,
                osd_w_encoding_inputs,
                encoding_input_count
        );

    else if(osd_method==2){
        osd_0(
                H,
                synd,
                osd0_decoding,
                soft_decisions,
                rank
        );
        osdw_decoding=osd0_decoding;
    }

    else if(osd_method==3){
        osd_g(
                H,
                synd,
                osd0_decoding,
                osdw_decoding,
                soft_decisions,
                osd_order,
                rank
        );
    }

    else {
        cout << "ERROR. <Function bp_osd::osd_post_process>. Invalid OSD method! Method: "<<osd_method << endl;
        exit(22);
    }

    return osdw_decoding;

}


char *bp_osd::bp_osd_decode(char *synd) {

    bp_decoding=bp_decode(synd);

    if(*converge){
        for(int i=0;i<N;i++) {
            osd0_decoding[i] = bp_decoding[i];
            osdw_decoding[i] = bp_decoding[i];
        }
        return bp_decoding;
    }

    if((osd_method==0)||(osd_method==1))
        osd_w(
                H,
                synd,
                osd0_decoding,
                osdw_decoding,
                log_prob_ratios,
                osd_order, rank,
                osd_w_encoding_inputs,
                encoding_input_count
        );

    else if(osd_method==2){
        osd_0(
                H,
                synd,
                osd0_decoding,
                log_prob_ratios,
                rank
                );
        osdw_decoding=osd0_decoding;
    }

    else if(osd_method==3){
        osd_g(
                H,
                synd,
                osd0_decoding,
                osdw_decoding,
                log_prob_ratios,
                osd_order,
                rank
        );
    }

    else {
        cout << "ERROR. <Function bp_osd::bp_osd_decode>. Invalid OSD method! Method: "<<osd_method << endl;
        exit(22);
    }

//    test_correct_synd(osd0_decoding,synd);
//    test_correct_synd(osdw_decoding,synd);

//    osd0_decoding=osd_data[0].decoding;
//    osdw_decoding=osd_data[1].decoding;

    return osdw_decoding;

}

void bp_osd::test_correct_synd(char *decoding, char *synd){
    char *test_synd=new char[M]();
//    cout<<"testing"<<endl;
    syndrome(H,decoding,test_synd);
    assert(bin_char_equal(synd,test_synd,M)==1);
    delete[] test_synd;
}

