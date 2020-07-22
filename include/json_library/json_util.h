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

json json_read_safe(json j, string key,string default_to="0");

#endif //JOSCHKA_TEST_PERIODIC_JSON_UTIL_H
