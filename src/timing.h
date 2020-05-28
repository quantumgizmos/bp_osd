//
// Created by joschka on 16/03/2020.
//

#ifndef TIMING_H
#define TIMING_H

#include <iostream>
#include <stdio.h>
#include <memory>
#include <vector>
#include<time.h>
#include <ctime>
#include<chrono>
#include <string.h>

using namespace std;

class timing {

    public:
        chrono::time_point<std::chrono::system_clock> start, current;
        time_t start_time,current_time;
        double diff_t;
        char runtime[100];
        char *start_time_string;
        int start_time_seconds;


        timing();
        double elapsed_time_seconds();
        char *elapsed_time_readable();

};




#endif //TIMING_H
