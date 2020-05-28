#include "timing.h"

timing::timing(){
    start = chrono::system_clock::now();
    start_time = chrono::system_clock::to_time_t(start);
    start_time_string=strtok(ctime(&start_time),"\n");
    start_time_seconds=int(start_time);
}


double timing::elapsed_time_seconds(){

    current = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = current - start;
    return elapsed_seconds.count();

}

char *timing::elapsed_time_readable(){

    double diff_time=elapsed_time_seconds();
    int delta_time=int(diff_time);

    int runtime_hours=delta_time/3600;
    int runtime_minutes=(delta_time % 3600)/60;
    int runtime_seconds=(delta_time % 3600)%60;
    sprintf(runtime,"%ih:%im:%is",runtime_hours,runtime_minutes,runtime_seconds);

    return runtime;

}
