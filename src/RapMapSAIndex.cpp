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

#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
#include "RapMapSAIndex.hpp"
#include "IndexHeader.hpp"
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>


#include <future>
#include <thread>

/*
void set_empty_key(spp::sparse_hash_map<uint64_t,
                       rapmap::utils::SAInterval<IndexT>,
		   rapmap::utils::KmerKeyHasher>& khash,
		   uint64_t k) {
}

void set_empty_key(spp::sparse_hash_map<uint64_t,
                       rapmap::utils::SAInterval<IndexT>,
		   rapmap::utils::KmerKeyHasher>& khash,
		   uint64_t k) {
}
*/

    // Set the SA and text pointer if this is a perfect hash
template <typename IndexT>
void setPerfectHashPointers(RegHashT<uint64_t,
                            rapmap::utils::SAInterval<IndexT>,
                            rapmap::utils::KmerKeyHasher>& khash, std::vector<IndexT>& SA, std::string& seq) {
    // do nothing
}

template <typename IndexT>
void setPerfectHashPointers(PerfectHashT<uint64_t,
                            rapmap::utils::SAInterval<IndexT>>& khash, std::vector<IndexT>& SA, std::string& seq) {
    khash.setSAPtr(&SA);
    khash.setTextPtr(seq.c_str(), seq.length());
}

// These are **free** functions that are used for loading the
// appropriate type of hash.
template <typename IndexT>
bool loadHashFromIndex(const std::string& indexDir,
                       RegHashT<uint64_t,
                       rapmap::utils::SAInterval<IndexT>,
                       rapmap::utils::KmerKeyHasher>& khash) {
      std::ifstream hashStream(indexDir + "hash.bin", std::ios::binary);
      khash.unserialize(typename spp_utils::pod_hash_serializer<uint64_t, rapmap::utils::SAInterval<IndexT>>(),
			&hashStream);
      return true;
}

template <typename IndexT>
bool loadHashFromIndex(const std::string& indexDir,
		       PerfectHashT<uint64_t, rapmap::utils::SAInterval<IndexT>> & h) {
    std::string hashBase = indexDir + "hash_info";
    h.load(hashBase);
    return true;
}

template <typename IndexT, typename HashT>
RapMapSAIndex<IndexT, HashT>::RapMapSAIndex() {}

// Given a position, p, in the concatenated text,
// return the corresponding transcript
template <typename IndexT, typename HashT>
IndexT RapMapSAIndex<IndexT, HashT>::transcriptAtPosition(IndexT p) {
    return rankDict->rank(p);
}

// Return true if the corresponding bit in the decoyArray is set
template <typename IndexT, typename HashT>
bool RapMapSAIndex<IndexT, HashT>::isDecoy(IndexT p) {
  return (numDecoys > 0 and p >= static_cast<int64_t>(firstDecoyIndex));
  //return (decoyArray) ? bit_array_get_bit(decoyArray.get(), p) : false;
}

// Return true if the corresponding bit in the decoyArray is set
template <typename IndexT, typename HashT>
uint64_t RapMapSAIndex<IndexT, HashT>::getNumDecoys() {
  return numDecoys;
}

template <typename IndexT, typename HashT>
bool RapMapSAIndex<IndexT, HashT>::load(const std::string& indDir) {

    auto logger = spdlog::get("stderrLog");

    IndexHeader h;
    std::ifstream indexStream(indDir + "header.json");
    {
      cereal::JSONInputArchive ar(indexStream);
      ar(h);
    }
    indexStream.close();
    uint32_t idxK = h.kmerLen();
    rapmap::utils::my_mer::k(idxK);

    // This part takes the longest, so do it in it's own asynchronous task
    std::future<bool> loadingHash = std::async(std::launch::async, [this, logger, indDir]() -> bool {
	return loadHashFromIndex(indDir, khash);
    });

    std::ifstream saStream(indDir + "sa.bin");
    {
        logger->info("Loading Suffix Array ");
        cereal::BinaryInputArchive saArchive(saStream);
        saArchive(SA);
        //saArchive(LCP);
    }
    saStream.close();

    std::ifstream seqStream(indDir + "txpInfo.bin");
    {
        logger->info("Loading Transcript Info ");
        cereal::BinaryInputArchive seqArchive(seqStream);
        seqArchive(txpNames);
        seqArchive(numDecoys);
        seqArchive(firstDecoyIndex);
        seqArchive(txpOffsets);
        //seqArchive(positionIDs);
        seqArchive(seq);
        seqArchive(txpCompleteLens);
    }
    seqStream.close();

    std::string rsFileName = indDir + "rsd.bin";
    FILE* rsFile = fopen(rsFileName.c_str(), "r");
    {
        logger->info("Loading Rank-Select Bit Array");
        bitArray.reset(bit_array_create(0));
        if (!bit_array_load(bitArray.get(), rsFile)) {
            logger->error("Couldn't load bit array from {}!", rsFileName);
            std::exit(1);
        }
        logger->info("There were {:n} set bits in the bit array", bit_array_num_bits_set(bitArray.get()));
        rankDict.reset(new rank9b(bitArray->words, bitArray->num_of_bits));
    }
    fclose(rsFile);

    {
        logger->info("Computing transcript lengths");
        txpLens.resize(txpOffsets.size());
        if (txpOffsets.size() > 1) {
            for(size_t i = 0; i < txpOffsets.size() - 1; ++i) {
                auto nextOffset = txpOffsets[i+1];
                auto currentOffset = txpOffsets[i];
                txpLens[i] = (nextOffset - 1) - currentOffset;
            }
        }
        // The last length is just the length of the suffix array - the last offset
        txpLens[txpOffsets.size()-1] = (SA.size() - 1) - txpOffsets[txpOffsets.size() - 1];
    }

    /*
    std::string decoyBVFileName = indDir + "isdecoy.bin";
    FILE* decoyBVFile = fopen(decoyBVFileName.c_str(), "r");
    {
      logger->info("Loading decoy bitvector.");
      decoyArray.reset(bit_array_create(0));
      if (!bit_array_load(decoyArray.get(), decoyBVFile)) {
        logger->error("Couldn't load decoy vector from {}!", decoyBVFileName);
        std::exit(1);
      }
      logger->info("Finished loading decoy bitvector.");
    }
    */

    logger->info("Waiting to finish loading hash");
    loadingHash.wait();
    auto hashLoadRes = loadingHash.get();
    if (!hashLoadRes) {
        logger->error("Failed to load hash!");
        std::exit(1);
    }
    // Set the SA and text pointer if this is a perfect hash
    setPerfectHashPointers(khash, SA, seq); 
    logger->info("Done loading index");
    return true;
}

template class RapMapSAIndex<int32_t,  RegHashT<uint64_t,
                      rapmap::utils::SAInterval<int32_t>,
                      rapmap::utils::KmerKeyHasher>>;
template class RapMapSAIndex<int64_t,  RegHashT<uint64_t,
                      rapmap::utils::SAInterval<int64_t>,
                      rapmap::utils::KmerKeyHasher>>;
template class RapMapSAIndex<int32_t, PerfectHashT<uint64_t, rapmap::utils::SAInterval<int32_t>>>;
template class RapMapSAIndex<int64_t, PerfectHashT<uint64_t, rapmap::utils::SAInterval<int64_t>>>;
