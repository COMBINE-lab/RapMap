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

#ifndef __RAPMAP_SA_INDEX_HPP__
#define __RAPMAP_SA_INDEX_HPP__

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"

#include "bit_array.h"
//#include "bitmap.h"
//#include "shared.h"
#include "rank9b.h"

#include <cstdio>
#include <vector>
#include <memory>

#include <fstream>
#include "RapMapUtils.hpp"

template <typename IndexT, typename HashT>
class RapMapSAIndex {
    public:
    using IndexType = IndexT;
    using HashType = HashT;

      struct BitArrayDeleter {
        void operator()(BIT_ARRAY* b) {
          if(b != nullptr) {
            bit_array_free(b);
          }
        }
      };

	  using BitArrayPointer = std::unique_ptr<BIT_ARRAY, BitArrayDeleter>;

    RapMapSAIndex();

  	// Given a position, p, in the concatenated text,
  	// return the corresponding transcript
  	IndexT transcriptAtPosition(IndexT p);
    bool isDecoy(IndexT p);
    uint64_t getNumDecoys();
    bool load(const std::string& indDir);

    std::vector<IndexT> SA;

    BitArrayPointer bitArray{nullptr};
    uint64_t numDecoys{0};
    uint64_t firstDecoyIndex{std::numeric_limits<uint64_t>::max()};
  //BitArrayPointer decoyArray{nullptr};
    std::unique_ptr<rank9b> rankDict{nullptr};

    std::string seq;
    std::vector<std::string> txpNames;
    std::vector<IndexT> txpOffsets;
    std::vector<IndexT> txpLens;
    std::vector<IndexT> positionIDs;
    std::vector<uint32_t> txpCompleteLens;
    std::vector<rapmap::utils::SAIntervalWithKey<IndexT>> kintervals;
    HashT khash;
};

#endif //__RAPMAP_SA_INDEX_HPP__
