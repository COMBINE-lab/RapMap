#ifndef __SCOPED_TIMER_HPP__
#define __SCOPED_TIMER_HPP__
// from https://gist.github.com/justgord/4482447
#include <chrono>
#include <iostream>

struct ScopedTimer
{
    std::chrono::high_resolution_clock::time_point t0;

    ScopedTimer()
        : t0(std::chrono::high_resolution_clock::now())
    { }
    ~ScopedTimer(void)
    {
        auto  t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedSec =  t1 - t0;
        std::cerr << "Elapsed time: " << elapsedSec.count() << "s\n";
    }
};

#endif //__SCOPED_TIMER_HPP__
