#ifndef _TIMING_HPP_
#define _TIMING_HPP_

#include <chrono>

using std::chrono::high_resolution_clock;

#define TIMING_INIT std::chrono::time_point<std::chrono::high_resolution_clock> a_time1, a_time2; std::chrono::duration<double> a_time12; double a_time;

#define TIMING_START a_time1 = high_resolution_clock::now();

#define TIMING_END \
{ \
  a_time2 = high_resolution_clock::now(); \
  a_time12 = a_time2 - a_time1; \
  a_time = a_time12.count(); \
};

#define TIMING_FINISH \
{ \
  a_time2 = high_resolution_clock::now(); \
  a_time12 = a_time2 - a_time1; \
  a_time = a_time12.count(); \
};

#define TIMING_STOP \
{ \
  a_time2 = high_resolution_clock::now(); \
  a_time12 = a_time2 - a_time1; \
  a_time = a_time12.count(); \
};

//#define TIMING_START  {a_time = omp_get_wtime();};
//#define TIMING_FINISH {a_time = omp_get_wtime() - a_time;};
//#define TIMING_STOP   {a_time = omp_get_wtime() - a_time;};

#endif
