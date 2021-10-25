//
// Created by joschka on 28/04/2020.
//

#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>
#include <ctime>

#include "bp_osd_c.hpp"
using namespace std;




int main(int argc, char *argv[])
{

    if(argv[1]==nullptr){
        cout<<"ERROR: Please provide a location for the input directory"<<endl;
        exit(22);
    };

    if(argv[2]==nullptr){
        cout<<"ERROR: Please provide a location for the output directory"<<endl;
        exit(22);
    };

//    if(argv[3]==nullptr){
//        cout<<"ERROR: Please select a method for verifying the recovered error"<<endl;
//        exit(22);
//    } else if ( strcmp(argv[3],"-lx") != 0 && strcmp(argv[3],"-hz") != 0  ) {
//        cout<<"ERROR: Please select a valid method for verifying the recovered error"<<endl;
//        exit(22);
//    };

    //timing functions setup
    timing simtime;
    cout<<"Start time: "<<simtime.start_time_string<<endl;
    cout<<"Input file: "<<argv[1]<<endl;
    int elapsed_seconds;
    elapsed_seconds=simtime.elapsed_time_seconds();

    //read in json input file
    ifstream json_input_file(argv[1]);
    json json_input;
    json_input_file>>json_input;
    cout<<json_input.dump(1)<<endl;
    json output=json_input;

    //get input file from command line arguments
    output["input_file"]=argv[1];

    //setup random number generator. If seed is <=0 then use current time for seed.
    int random_seed = json_read_safe(output,"input_seed");
    if (random_seed <= 0) {
      time_t now = time(0);
      random_seed = (int)now;
    }
    MTRand r=setup_mtwister_rng(random_seed);
    output["seed"]=random_seed;

    //save start time
    output["start_time"] = simtime.start_time_string;

    //generate UI from start time + random seed
    unsigned long long int ui =simtime.start_time_seconds+random_seed;
    output["ui"]=ui;

    //read in sim parameters
    double bit_error_rate=json_read_safe(json_input,"bit_error_rate");
    int osd_order = json_read_safe(json_input,"osd_order");
    int max_iter=json_read_safe(json_input,"max_iter");
    long long unsigned int target_runs=json_read_safe(json_input,"target_runs");

    //read in Logical check method
    string logical_check_method = json_read_safe(json_input,"logical_check_method","lx");
    if((logical_check_method=="lx")||(logical_check_method=="hz")) 1;
    else{
        cout<<"ERROR: Please select a valid method for verifying the recovered error. `lx' or `hz'."<<endl;
        exit(22);
    }

    //read in OSD method and check that it is valid
    string osd_method = json_read_safe(json_input,"osd_method");
    int osd_method_i=-1;
    if((osd_method=="exhaustive")||(osd_method=="osd_e")) osd_method_i=0;
    else if((osd_method=="combination_sweep")||(osd_method=="osd_cs")) osd_method_i=1;
    else if(osd_method=="osd_0"){
        osd_method_i=2;
        output["osd_order"]=0;
    }
    else if(osd_method=="osd_g") osd_method_i=3;
    else{
        cout<<"ERROR: Invalid OSD method: "<<osd_method<<endl;
        exit(22);
    }
    if(osd_order==0){
        osd_method="osd_0";
        output["osd_method"]=osd_method;
        osd_method_i=2;
    }

    //read in BP method and check that it is valid
    string bp_method = json_read_safe(json_input,"bp_method","min_sum");
    int bp_method_i=-1;
    if(bp_method=="min_sum") bp_method_i=0;
    else if(bp_method=="min_sum_min_synd") bp_method_i=1;
    else if(bp_method=="product_sum") bp_method_i=2;
    else if(bp_method=="product_sum_min_synd") bp_method_i=3;
    else{
        cout<<"ERROR: Invalid BP method: "<<bp_method<<endl;
        exit(22);
    }
    output["bp_method"]=bp_method;

    

    //set up output file
    string p_label =to_string(bit_error_rate);
    replace( p_label.begin(), p_label.end(), '.', '_');
    string label =output["code_label"];
    stringstream output_filename;
    output_filename<<argv[2]<<"/"<<label<<";p:"<<p_label<<";ui:"<<ui<<".json";
    cout<<"\nOutput filename: "<< output_filename.str() <<endl;


    //LOAD ALIST FILES
    cout<<"\nLoading alist files"<<endl;

    //load hx parity check matrix
    mod2sparse *hx=load_alist_cpp(json_read_safe(json_input,"hx_filename"));
    int hx_m=mod2sparse_rows(hx);
    int hx_n=mod2sparse_cols(hx);
    output["hx_n"]=hx_n;
    output["hx_m"]=hx_m;

    //load logical check matrix, lx
    mod2sparse *lx;
    if(logical_check_method=="lx") {
        lx = load_alist_cpp(json_read_safe(json_input, "lx_filename"));
        int lx_k = mod2sparse_rows(lx);
        int lx_n = mod2sparse_cols(lx);
        output["lx_n"] = lx_n;
        output["lx_k"] = lx_k;
        assert(lx_n == hx_n);
    }

    mod2sparse *hz;
    if(logical_check_method=="hz") {
        hz=load_alist_cpp(json_read_safe(json_input,"hz_filename"));
        int hz_m=mod2sparse_rows(hz);
        int hz_n=mod2sparse_cols(hz);
        output["hz_n"]=hz_n;
        output["hz_m"]=hz_m;
        assert(hz_n==hx_n);
    }


    //Setup BP+OSD decoder for hx
    bp_osd hx_bp_osd(hx,bit_error_rate,max_iter,osd_order,osd_method_i,bp_method_i);
    output["max_iter"]=hx_bp_osd.max_iter;
    output["osdw_encoding_operator_count"]=hx_bp_osd.encoding_input_count;


    //MEMORY ALLOCATION (we can do this now that we know the size of the matrices we are dealing with!)
    double *soft_decisions = new double[hx_n]();
    char *bit_error = new char[hx_n]();
    char *bp_decoding = new char[hx_n]();
    char *osd_decoding = new char[hx_n]();
    char *osdw_decoding = new char[hx_n]();
    char *bit_syndrome = new char[hx_m]();

    //declarations for various simulation variables
    long long unsigned int bp_converge_count=0;
    long long unsigned int bp_success_count=0;
    long long unsigned int osd0_success_count=0;
    long long unsigned int osdw_success_count=0;
    bool logical_error;
    bool zero_error;
    double osdw_ler;
    double osd0_ler;
    double bp_ler;
    double average_iterations_before_converge=-1;
    int min_iterations_before_converge=hx_n;
    int max_iterations_before_converge=-1;

    //MAIN SIM LOOP

    cout<<endl;
    cout<<"Simulating "<<target_runs<<" error correction cycles..."<<endl;
    for(long long unsigned int run_count=1; run_count <= target_runs; run_count++) {

        //generate error
        zero_error=true;
        for (int bit_no = 0; bit_no < hx_n; bit_no++) {
            bit_error[bit_no] = genRand(&r) < bit_error_rate;
            if(bit_error[bit_no]==1) zero_error=false;
        }

        //Skip decoding if generated error is a zero vector.
        if(zero_error) {
            bp_success_count++;
            osd0_success_count++;
            osdw_success_count++;
        }
        else {
            //calculate syndrome
            syndrome(hx, bit_error, bit_syndrome);

            //decode
            osdw_decoding = hx_bp_osd.bp_osd_decode(bit_syndrome);


            //check for logical errors
            if (logical_check_method=="lx") {
                logical_error = check_logical_error_lx(lx, bit_error, osdw_decoding);
            } else {
                logical_error = check_logical_error_hz(hz, bit_error, osdw_decoding);
            }

            if (!logical_error) osdw_success_count++;
            
//            osdw_decoding = hx_bp_osd.osd0_decoding;
            if (logical_check_method=="lx") {
                logical_error = check_logical_error_lx(lx, bit_error, hx_bp_osd.osd0_decoding);
            } else {
                logical_error = check_logical_error_hz(hz, bit_error, hx_bp_osd.osd0_decoding);
            }
            
            if (!logical_error) osd0_success_count++;

            if (*hx_bp_osd.converge == 1) {
                bp_converge_count++;

		if (*hx_bp_osd.iter>max_iterations_before_converge){
                    max_iterations_before_converge=*hx_bp_osd.iter+1;
                }

                if (*hx_bp_osd.iter<min_iterations_before_converge){
                    min_iterations_before_converge=*hx_bp_osd.iter+1;
                }

                if (average_iterations_before_converge==-1){
                    average_iterations_before_converge=*hx_bp_osd.iter+1;
                }
                else{

                    average_iterations_before_converge=((bp_converge_count-1)*average_iterations_before_converge+*hx_bp_osd.iter+1)/bp_converge_count;
                }

                if (logical_check_method=="lx") {
                    logical_error = check_logical_error_lx(lx, bit_error, hx_bp_osd.bp_decoding);
                } else {
                    logical_error = check_logical_error_hz(hz, bit_error, hx_bp_osd.bp_decoding);
                }

                if (!logical_error) bp_success_count++;
            }
        }


        //write to output every 3 seconds. Change this as you see fit
        if((simtime.elapsed_time_seconds()-elapsed_seconds>3) | (run_count==target_runs) ){

            //calculate logical error rates
            osdw_ler= 1.0 - osdw_success_count / (double(run_count));
            osd0_ler= 1.0 - osd0_success_count / (double(run_count));
            bp_ler= 1.0 - bp_success_count / (double(run_count));

            elapsed_seconds=simtime.elapsed_time_seconds();

            output["run_count"] = run_count;
            output["osdw_success_count"]=osdw_success_count;
            output["osd0_success_count"]=osd0_success_count;
            output["bp_success_count"]=bp_success_count;
            output["bp_converge_count"]=bp_converge_count;
            output["osdw_ler"] = osdw_ler;
            output["osd0_ler"] = osd0_ler;
            output["bp_ler"] = bp_ler;
            output["runtime_seconds"] = simtime.elapsed_time_seconds();
            output["runtime_readable"] = simtime.elapsed_time_readable();
            output["max_iterations_before_converge"]=max_iterations_before_converge;
            output["min_iterations_before_converge"]=min_iterations_before_converge;
            output["average_iterations_before_converge"]=average_iterations_before_converge;
            cout << "Runs: " << run_count << "; OSDW_LER: " << osdw_ler << "; OSD0_LER: " << osd0_ler << "; BP_LER: " << bp_ler << "; Iter. (min/avg/max)" << min_iterations_before_converge << "/"
		 << average_iterations_before_converge << "/" << max_iterations_before_converge <<
	      "; Runtime: " << output["runtime_readable"] << endl;

            //command line output
            ofstream output_file(output_filename.str(),ofstream::trunc);
            output_file<<output.dump(1);
            output_file.close();

        }





    }



}
