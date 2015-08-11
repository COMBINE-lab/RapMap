#ifndef __SINGLE_ALIGNMENT_FORMATTER_HPP__
#define __SINGLE_ALIGNMENT_FORMATTER_HPP__

#include "RapMapUtils.hpp"

template <typename IndexPtrT>
struct SingleAlignmentFormatter {
    SingleAlignmentFormatter(IndexPtrT indexIn) : index(indexIn),
    readTemp(1000, 'A'),
    qualTemp(1000, '~'),
    cigarStr(buff, 1000){
    }

    // Data members
    IndexPtrT index;
    std::string readTemp;
    std::string qualTemp;
    char buff[1000];
    rapmap::utils::FixedWriter cigarStr;
};

#endif //__PAIR_ALIGNMENT_FORMATTER_HPP__
