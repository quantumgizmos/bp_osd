//
// Created by joschka on 16/03/2020.
//

#ifndef JOSCHKA_SINGLE_SHOT_TEST_TIMING_H
#define JOSCHKA_SINGLE_SHOT_TEST_TIMING_H

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

    timing(){

        start = chrono::system_clock::now();
        start_time = chrono::system_clock::to_time_t(start);
        start_time_string=strtok(ctime(&start_time),"\n");
        cout<<start_time_string<<endl;
        start_time_seconds=int(start_time);
        cout<<start_time_seconds<<endl;

    }


    double elapsed_time_seconds(){

        current = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = current - start;
        return elapsed_seconds.count();

    }

    char *elapsed_time_readable(){

        double diff_time=elapsed_time_seconds();
        int delta_time=int(diff_time);

        int runtime_hours=delta_time/3600;
        int runtime_minutes=(delta_time % 3600)/60;
        int runtime_seconds=(delta_time % 3600)%60;
        sprintf(runtime,"%ih:%im:%is       ",runtime_hours,runtime_minutes,runtime_seconds);

//        cout<<"Runtime: " <<diff_time<<endl;

        return runtime;

    }





};


#endif //JOSCHKA_SINGLE_SHOT_TEST_TIMING_H
