#ifndef __PAIR_ALIGNMENT_FORMATTER_HPP__
#define __PAIR_ALIGNMENT_FORMATTER_HPP__

#include "RapMapUtils.hpp"

template <typename IndexPtrT>
struct PairAlignmentFormatter {
    PairAlignmentFormatter(IndexPtrT indexIn) : index(indexIn),
    read1Temp(1000, 'A'),
    qual1Temp(1000, '~'),
    read2Temp(1000, 'A'),
    qual2Temp(1000, '~'),
    cigarStr1(buff1, 1000),
    cigarStr2(buff2, 1000) {
    }

    // Data members
    IndexPtrT index;
    std::string read1Temp;
    std::string qual1Temp;
    std::string read2Temp;
    std::string qual2Temp;
    char buff1[1000];
    char buff2[1000];
    rapmap::utils::FixedWriter cigarStr1;
    rapmap::utils::FixedWriter cigarStr2;
};

#endif //__PAIR_ALIGNMENT_FORMATTER_HPP__
