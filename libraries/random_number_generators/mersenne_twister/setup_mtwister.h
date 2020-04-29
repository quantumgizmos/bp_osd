//
// Created by joschka on 26/02/2020.
//

#ifndef JOSCHKA_TEST3_SETUP_MTWISTER_H
#define JOSCHKA_TEST3_SETUP_MTWISTER_H

#include <iostream>
#include <stdio.h>
#include <memory>
#include <vector>
#include<time.h>

// basic file operations
#include <iostream>
#include <fstream>

extern "C" {
    #include "mtwister.h"
}


using namespace std;


unsigned long long int get_random_seed(){
    int iterator;
    unsigned long int random_seed;
    ifstream infile;
    ofstream outfile;

    infile.open("random_number_generators/mersenne_twister/random_org_seeds_iterator.txt");
    infile >> iterator;

//    cout<<"hello world"<<endl;
//    cout<<iterator<<endl;
    infile.close();

    iterator++;
    outfile.open("random_number_generators/mersenne_twister/random_org_seeds_iterator.txt");
    outfile << iterator << "\n";
    outfile.close();

    infile.open("random_number_generators/mersenne_twister/random_org_seeds.txt");

    for (int j = 0; j <= iterator; j++) infile >> random_seed;

    cout << random_seed << endl;

    infile.close();

    return random_seed;


}

MTRand setup_mtwister_rng(unsigned long long int random_seed){
    MTRand r;
    if(random_seed==0) random_seed=get_random_seed(); //gets random seed from file
    r=seedRand(random_seed);
    cout<<"Mtwister seed: "<<random_seed<<endl;
    return r;
}




#endif //JOSCHKA_TEST3_SETUP_MTWISTER_H
