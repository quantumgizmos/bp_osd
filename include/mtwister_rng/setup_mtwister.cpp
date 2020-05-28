#include "setup_mtwister.h"

MTRand setup_mtwister_rng(unsigned long long int random_seed){
    MTRand r;
    r=seedRand(random_seed);
    cout<<"Mtwister seed: "<<random_seed<<endl;
    return r;
}