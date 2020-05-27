//
// Created by joschka on 24/02/2020.
//

#ifndef SIM_UTIL_H
#define SIM_UTIL_H

extern "C" {
#include "mtwister.h"
}

int check_logical_error(mod2sparse *l,char *orginal_error,char *decoding);
int gen_error(char *bit_error, int error_len, double bit_error_rate, MTRand *r);
mod2sparse *load_alist_cpp(std::string filename);

#endif //SIM_UTIL_H
