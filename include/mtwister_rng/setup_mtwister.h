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

MTRand setup_mtwister_rng(unsigned long long int random_seed){
    MTRand r;
    r=seedRand(random_seed);
    cout<<"Mtwister seed: "<<random_seed<<endl;
    return r;
}




#endif //JOSCHKA_TEST3_SETUP_MTWISTER_H
