#ifndef Stopwatch_header
#define Stopwatch_header

#include <chrono>

// Measures elapsed wall-clock time. clock() would add up the CPU time of
// all OpenMP threads.
class Stopwatch
{
  private:
    std::chrono::steady_clock::time_point start, end;
  public:
    Stopwatch() {start = std::chrono::steady_clock::now(); end = start;}
    void tic() {start = std::chrono::steady_clock::now();}
    void toc() {end = std::chrono::steady_clock::now();}
    double takeTime() {
        return std::chrono::duration<double>(end - start).count();
    }
};

#endif

/*-----------------------------------------------------------------------
  Usage:
  Declare a class as:
    Stopwatch sw;
  Then the time a piece of code takes can be recorded as:
    sw.tic();
    **** code ****
    sw.toc();
  And the result can be outputted using:
    cout << sw.takeTime() << endl;
-----------------------------------------------------------------------*/
