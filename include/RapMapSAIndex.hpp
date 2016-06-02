#ifndef __RAPMAP_SA_INDEX_HPP__
#define __RAPMAP_SA_INDEX_HPP__

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"

#include "google/dense_hash_map"
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

    bool load(const std::string& indDir);

    std::vector<IndexT> SA;

    BitArrayPointer bitArray{nullptr};
    std::unique_ptr<rank9b> rankDict{nullptr};

    std::string seq;
    std::vector<std::string> txpNames;
    std::vector<IndexT> txpOffsets;
    std::vector<IndexT> txpLens;
    std::vector<IndexT> positionIDs;
    std::vector<rapmap::utils::SAIntervalWithKey<IndexT>> kintervals;
    HashT khash;
};

#endif //__RAPMAP_SA_INDEX_HPP__
