//
// Created by joschka on 17/04/2020.
//

#ifndef JOSCHKA_TEST_PERIODIC_JSON_UTIL_H
#define JOSCHKA_TEST_PERIODIC_JSON_UTIL_H


#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

#include <json.hpp>
using json = nlohmann::json;

auto json_read_safe(json j, string key){

    try {
        auto data_item = j[key];

        if(data_item==nullptr) {
            throw 22;
        }

        return data_item;

    }
    catch (int e)
    {
        // output exception information
        cout<<"Error. Key '"<<key<<"' does not exist!"<<endl;
        exit(23);
    }

}





#endif //JOSCHKA_TEST_PERIODIC_JSON_UTIL_H
