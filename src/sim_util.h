//
// Created by joschka on 24/02/2020.
//

#ifndef JOSCHKA_TEST3_SIM_UTIL_H
#define JOSCHKA_TEST3_SIM_UTIL_H

int check_logical_error(mod2sparse *l,char *orginal_error,char *decoding);
int gen_error(char *bit_error, int error_len, double bit_error_rate, MTRand *r);
mod2sparse *load_alist_cpp(std::string filename);

#endif //JOSCHKA_TEST3_SIM_UTIL_H
