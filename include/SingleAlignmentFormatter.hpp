//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

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
