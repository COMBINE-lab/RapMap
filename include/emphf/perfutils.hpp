#pragma once

#include <map>
#include <sys/time.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cmath>

namespace emphf {

    inline double get_time_usecs() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return double(tv.tv_sec) * 1000000 + double(tv.tv_usec);
    }

    // stolen from folly
    template <class T>
    inline void do_not_optimize_away(T&& datum) {
        asm volatile("" : "+r" (datum));
    }

    struct stats_accumulator {
        stats_accumulator()
            : m_n(0)
            , m_mean(0)
            , m_m2(0)
        {}

        void add(double x)
        {
            m_n += 1;
            auto delta = x - m_mean;
            m_mean += delta / m_n;
            m_m2 += delta * (x - m_mean);
        }

        double mean() const
        {
            return m_mean;
        }

        double variance() const
        {
            return m_m2 / (m_n - 1);
        }

        double relative_stddev() const
        {
            return std::sqrt(variance()) / mean() * 100;
        }

    private:
        double m_n;
        double m_mean;
        double m_m2;
    };

}
