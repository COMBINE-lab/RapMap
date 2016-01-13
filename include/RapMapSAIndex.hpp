#ifndef __RAPMAP_SA_INDEX_HPP__
#define __RAPMAP_SA_INDEX_HPP__

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "jellyfish/file_header.hpp"
#include "jellyfish/binary_dumper.hpp"
#include "jellyfish/hash_counter.hpp"
#include "jellyfish/mapped_file.hpp"
#include "JFRaw.hpp"

#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"
#include "RSDic.hpp"
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

template <typename IndexT>
class RapMapSAIndex {
    using FileMerArray = jellyfish::large_hash::array_raw<rapmap::utils::my_mer>;

    public:
    using IndexType = IndexT;

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

    // Given a key (k-mer), return a pointer to the SA
    // interval to which it corresponds.
    // NOTE: Returns nullptr if the key is not found.
    inline rapmap::utils::SAInterval<IndexT>*
        intervalForKmer(rapmap::utils::my_mer& key) {
            uint64_t val;
            size_t kID;
            bool found = khash->get_val_for_key(
                    key,
                    &val,
                    searchBuffer_,
                    &kID);
            return found ? &saIntervals[val] : nullptr;
        }

    bool load(const std::string& indDir);

    std::vector<IndexT> SA;
    //rsdic::RSDic rankDictSafe;

    BitArrayPointer bitArray{nullptr};
    std::unique_ptr<rank9b> rankDict{nullptr};

    std::string seq;
    std::vector<std::string> txpNames;
    std::vector<IndexT> txpOffsets;
    std::vector<IndexT> txpLens;
    std::vector<IndexT> positionIDs;

    std::unique_ptr<char> rawHashMem{nullptr};
    std::unique_ptr<FileMerArray> khash{nullptr};
    std::vector<rapmap::utils::SAInterval<IndexT>> saIntervals;

    private:
    rapmap::utils::my_mer searchBuffer_;
    /*
    google::dense_hash_map<uint64_t,
                        rapmap::utils::SAInterval<IndexT>,
                        rapmap::utils::KmerKeyHasher> khash;
                        */
                       // ::NopointerSerializer(), &hashStream);
        /*
    std::unordered_map<uint64_t,
                       rapmap::utils::SAInterval,
                       rapmap::utils::KmerKeyHasher> khash;
                       */
};

#endif //__RAPMAP_SA_INDEX_HPP__
