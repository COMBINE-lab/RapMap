//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016, 2017 Rob Patro, Avi Srivastava, Hirak Sarkar
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

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <random>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <set>

#include "tclap/CmdLine.h"

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
#include "xxhash.h"

#include "spdlog/spdlog.h"

#include "FastxParser.hpp"
// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"

#include "divsufsort.h"
#include "divsufsort64.h"

#include "RapMapFileSystem.hpp"
#include "RapMapUtils.hpp"
#include "ScopedTimer.hpp"
#include "bit_array.h"

#include "JFRaw.hpp"
#include "jellyfish/binary_dumper.hpp"
#include "jellyfish/file_header.hpp"
#include "jellyfish/hash_counter.hpp"
#include "jellyfish/mer_iterator.hpp"
#include "jellyfish/mer_overlap_sequence_parser.hpp"
#include "jellyfish/thread_exec.hpp"
#include "rank9b.h"

#include "IndexHeader.hpp"

// sha functionality
#include "picosha2.h"

#include <chrono>

using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using MerMapT = jellyfish::cooperative::hash_counter<rapmap::utils::my_mer>;

bool buildSA(const std::string& outputDir, std::string& concatText, size_t tlen,
             std::vector<int64_t>& SA) {
  // IndexT is the signed index type
  // UIndexT is the unsigned index type
  using IndexT = int64_t;
  using UIndexT = uint64_t;
  bool success{false};

  std::ofstream saStream(outputDir + "sa.bin", std::ios::binary);
  {
    ScopedTimer timer;
    SA.resize(tlen, 0);
    IndexT textLen = static_cast<IndexT>(tlen);
    std::cerr << "Building suffix array . . . ";
    auto ret = divsufsort64(
        reinterpret_cast<unsigned char*>(const_cast<char*>(concatText.data())),
        SA.data(), tlen);

    success = (ret == 0);
    if (success) {
      std::cerr << "success\n";
      {
        ScopedTimer timer2;
        std::cerr << "saving to disk . . . ";
        cereal::BinaryOutputArchive saArchive(saStream);
        saArchive(SA);
        std::cerr << "done\n";
      }
    } else {
      std::cerr << "FAILURE: return code from libdivsufsort64() was " << ret
                << "\n";
      saStream.close();
      std::exit(1);
    }
    std::cerr << "done\n";
  }
  saStream.close();
  return success;
}

// IndexT is the index type.
// int32_t for "small" suffix arrays
// int64_t for "large" ones
template <typename IndexT>
bool buildPerfectHash(const std::string& outputDir, std::string& concatText,
                      size_t tlen, uint32_t k, std::vector<IndexT>& SA,
                      uint32_t numHashThreads,
                      std::shared_ptr<spdlog::logger> log) {
  //BooMap<uint64_t, rapmap::utils::SAInterval<IndexT>> intervals;
  PerfectHashT<uint64_t, rapmap::utils::SAInterval<IndexT>> intervals;
  intervals.setSAPtr(&SA);
  intervals.setTextPtr(concatText.data(), concatText.length());

  // The start and stop of the current interval
  IndexT start = 0, stop = 0;
  // An iterator to the beginning of the text
  auto textB = concatText.begin();
  auto textE = concatText.end();
  // The current k-mer as a string
  rapmap::utils::my_mer mer;
  bool currentValid{false};
  std::string currentKmer;
  std::string nextKmer;
  while (stop < tlen) {
    // Check if the string starting at the
    // current position is valid (i.e. doesn't contain $)
    // and is <= k bases from the end of the string
    nextKmer = concatText.substr(SA[stop], k);
    if (nextKmer.length() == k and
        nextKmer.find_first_of('$') == std::string::npos) {
      // If this is a new k-mer, then hash the current k-mer
      if (nextKmer != currentKmer) {
        if (currentKmer.length() == k and
            currentKmer.find_first_of('$') == std::string::npos) {
          mer = rapmap::utils::my_mer(currentKmer);
          auto bits = mer.get_bits(0, 2 * k);
          intervals.add(std::move(bits), {start, stop});
          // push_back(std::make_pair<uint64_t,
          // rapmap::utils::SAInterval<IndexT>>(std::move(bits), {start,
          // stop}));
        }
        currentKmer = nextKmer;
        start = stop;
      }
    } else {
      // If this isn't a valid suffix (contains a $)
      // If the previous interval was valid, put it
      // in the hash.
      if (currentKmer.length() == k and
          currentKmer.find_first_of('$') == std::string::npos) {
        mer = rapmap::utils::my_mer(currentKmer);
        auto bits = mer.get_bits(0, 2 * k);
        // intervals.push_back(std::make_pair<uint64_t,
        // rapmap::utils::SAInterval<IndexT>>(std::move(bits), {start, stop}));
        intervals.add(std::move(bits), {start, stop});
      }
      // The current interval is invalid and empty
      currentKmer = nextKmer;
      start = stop;
    }
    if (stop % 1000000 == 0) {
      std::cerr << "\r\rprocessed " << stop << " positions";
    }
    // We always update the end position
    ++stop;
  }
  if (start < tlen) {
    if (currentKmer.length() == k and
        currentKmer.find_first_of('$') != std::string::npos) {
      mer = rapmap::utils::my_mer(currentKmer);
      auto bits = mer.get_bits(0, 2 * k);
      // intervals.push_back(std::make_pair<uint64_t,
      // rapmap::utils::SAInterval<IndexT>>(std::move(bits), {start, stop}));
      intervals.add(std::move(bits), {start, stop});
    }
  }

  // std::cerr << "\nthere are " << intervals.size() << " intervals of the
  // selected depth\n";

  std::cout << "building perfect hash function\n";
  intervals.build(numHashThreads);
  std::cout << "\ndone.\n";
  std::string outputPrefix = outputDir + "hash_info";
  std::cout << "saving the perfect hash and SA intervals to disk ... ";
  intervals.save(outputPrefix);
  std::cout << "done.\n";

  return true;
}

bool buildSA(const std::string& outputDir, std::string& concatText, size_t tlen,
             std::vector<int32_t>& SA) {
  // IndexT is the signed index type
  // UIndexT is the unsigned index type
  using IndexT = int32_t;
  using UIndexT = uint32_t;
  bool success{false};

  std::ofstream saStream(outputDir + "sa.bin", std::ios::binary);
  {
    ScopedTimer timer;
    SA.resize(tlen, 0);
    IndexT textLen = static_cast<IndexT>(tlen);
    std::cerr << "Building suffix array . . . ";
    auto ret = divsufsort(
        reinterpret_cast<unsigned char*>(const_cast<char*>(concatText.data())),
        SA.data(), tlen);

    success = (ret == 0);
    if (success) {
      std::cerr << "success\n";
      {
        ScopedTimer timer2;
        std::cerr << "saving to disk . . . ";
        cereal::BinaryOutputArchive saArchive(saStream);
        saArchive(SA);
        std::cerr << "done\n";
      }
    } else {
      std::cerr << "FAILURE: return code from libdivsufsort() was " << ret
                << "\n";
      saStream.close();
      std::exit(1);
    }
    std::cerr << "done\n";
  }
  saStream.close();
  return success;
}


template <typename IndexT>
uint32_t findLCPLength(std::string &concatText,size_t tlen, std::vector<IndexT>& SA,
		IndexT start,IndexT stop){
		uint32_t lcplen = 0 ;
		auto o1 = SA[start] + lcplen ;
		auto o2 = SA[stop] + lcplen ;
		auto maxIndex = std::max(o1,o2) ;
		while(maxIndex + lcplen < tlen and concatText[o1+lcplen] == concatText[o2+lcplen]){
			if(concatText[o1+lcplen] == '$') { break; }
			if(lcplen >= 255) { return 255; }
			++lcplen ;
		}
		return (lcplen) ;
}



template <typename IndexT, typename HashT>
bool updateSafe(std::string& concatText,
		  size_t tlen,
		  std::vector<IndexT>& SA,
		  rapmap::utils::kmerVal<IndexT>& val,
		  HashT& khash,
		  std::unique_ptr<rank9b>& rankDir,
                uint32_t k,
                std::vector<rapmap::utils::TGroup>& eqClasses
                ){
	uint32_t lcpLength = val.lcpLength ;
  uint32_t safeLCP = k;
  uint32_t maxSafeLCP = std::numeric_limits<decltype(val.safeLength)>::max();
  //assume readLen = 100
  uint32_t readLen = 100;
  maxSafeLCP = readLen ;

	if(lcpLength > k){
		auto start = val.interval.begin();
		auto stop = val.interval.end();
		rapmap::utils::my_mer mer,rcMer ;
		auto startIndex = SA[start];
    auto endIndex = startIndex + k - 1;
    uint32_t shift = 1;
    auto groundEqId = val.eqId;
    auto& groundEqClass = eqClasses[groundEqId];
    /*
		std::set<uint32_t> groundTidSet ;
		//make the set of transcripts
		for(auto i=start; i < stop ; ++i){
			groundTidSet.insert(rankDir->rank(SA[i])) ;
		}
    auto notFoundIt = groundTidSet.end();
    */
		//sort the transcript ids
		//std::sort(groundTidSet.begin(), groundTidSet.end()) ;

    mer = concatText.substr(startIndex, k);
    rcMer = mer.get_reverse_complement() ;

    auto mer_c = mer;
   		while((shift < lcpLength-k) and (startIndex + shift + k) < tlen){
      // for now, we don't care about longer safe LCPs
            auto hashItc = khash.find(mer.word(0)) ;

                        if (safeLCP >= maxSafeLCP) { safeLCP = maxSafeLCP; break; }


                        mer.shift_left(concatText[startIndex+shift+k-1]);
                        rcMer = mer.get_reverse_complement() ;


			//nextKmer = concatText.substr(startIndex+shift,k);
			//mer = nextKmer ;
			auto bits = mer.word(0);
			auto hashIt = khash.find(bits);
                        auto rcHashIt = khash.find(rcMer.word(0)) ;

			if(hashIt != khash.end()){
			    //find out the transcript set
                            auto& thisVal = hashIt->second ;
                            if (thisVal.eqId != groundEqId and !(groundEqClass.contains(eqClasses[thisVal.eqId]))) {
                                    val.safeLength = safeLCP + 1;


                                   return false;
                            }

                            if(rcHashIt != khash.end()){
                                /*
                                 auto& rcVal = rcHashIt->second ;
                                 if(rcVal.eqId != groundEqId and !(groundEqClass.contains(eqClasses[rcVal.eqId]))){
                                    val.safeLength = safeLCP;
                                    return false;
                                 }*/
                                val.safeLength = safeLCP + 1;
                                return false;

                            }
                        }else{
				//this case should never happen I mean come on
				//get real
                            //std::cerr << "\n" << mer << "\n" ;
                            std::cerr << "Impossible is the opposite of possible (https://www.youtube.com/watch?v=nAV0sxwx9rY)" << std::endl;
                            return false;
			}
			++shift;
                        ++safeLCP;
	            }

	}else{
		val.safeLength = safeLCP;
		return true;
	}

	val.safeLength = safeLCP + 1;
	return true;

}

template<typename IndexT>
bool build9kmerHash(const std::string& outputDir, std::string& concatText,
                    size_t tlen, uint32_t k, std::vector<IndexT>& SA) {
    using WordT = rapmap::utils::my_mer9::base_type ;
    std::vector<rapmap::utils::SAInterval<IndexT>> all9mers ;
    RegHashT<WordT, rapmap::utils::SAInterval<IndexT>,
                         rapmap::utils::KmerKeyHasher> khash;

    all9mers.resize(std::pow(4,9)) ;


    IndexT start = 0, stop = 0;
  // An iterator to the beginning of the text
  auto textB = concatText.begin();
  auto textE = concatText.end();
  // The current k-mer as a string
  rapmap::utils::my_mer9 mer;
  bool currentValid{false};
  std::string currentKmer;
  std::string nextKmer;
  while (stop < tlen) {
    // Check if the string starting at the
    // current position is valid (i.e. doesn't contain $)
    // and is <= k bases from the end of the string
    nextKmer = concatText.substr(SA[stop], k);
    /*
    if(SA[stop] + k < tlen){
        nextKmer = concatText.substr(SA[stop], k);
    }else{
        nextKmer = concatText.substr(SA[stop],tlen-SA[stop]);
    }*/
    if (nextKmer.length() == k and
        nextKmer.find_first_of('$') == std::string::npos) {
      // If this is a new k-mer, then hash the current k-mer
      if (nextKmer != currentKmer) {
        if (currentKmer.length() == k and
            currentKmer.find_first_of('$') == std::string::npos) {
          mer = currentKmer;
          auto bits = mer.word(0);
          auto hashIt = khash.find(bits);
          if (hashIt == khash.end()) {
            if (start > 1) {
              if (concatText.substr(SA[start - 1], k) ==
                  concatText.substr(SA[start], k)) {
                std::cerr << "T[SA[" << start - 1 << "]:" << k
                          << "] = " << concatText.substr(SA[start - 1], k)
                          << " = T[SA[" << start << "]:" << k << "]\n";
                std::cerr << "start = " << start << ", stop = " << stop << "\n";
                std::cerr << "[fatal (1)] THIS SHOULD NOT HAPPEN\n";
                std::exit(1);
              }
            }
            if (start == stop) {
              std::cerr << "[fatal (2)] Interval is empty! (start = " << start
                        << ") = (stop =  " << stop << ")\n";
            }
            if (start == stop) {
              std::cerr << "[fatal (3)] Interval is empty! (start = " << start
                        << ") = (stop =  " << stop << ")\n";
            }
            khash[bits] = {start,stop};
            all9mers[bits] = {start, stop};
            /*
            IndexT len = stop - start;
            bool overflow = (len >= std::numeric_limits<uint8_t>::max());
            uint8_t blen = overflow ? std::numeric_limits<uint8_t>::max() :
                static_cast<uint8_t>(len);
            khash[bits] = {start, blen};
            if (overflow) {
                overflowhash[bits] = len;
            }
            */
          } else {
            std::cerr << "\nERROR (1): trying to add same suffix "
                      << currentKmer << " (len = " << currentKmer.length()
                      << ") multiple times!\n";
            auto prevInt = hashIt->second;
            std::cerr << "existing interval is [" << prevInt.begin() << ", "
                      << prevInt.end() << ")\n";
            for (auto x = prevInt.begin(); x < prevInt.end(); ++x) {
              auto suff = concatText.substr(SA[x], k);
              for (auto c : suff) {
                std::cerr << "*" << c << "*";
              }
              std::cerr << " (len = " << suff.length() << ")\n";
            }
            std::cerr << "new interval is [" << start << ", " << stop << ")\n";
            for (auto x = start; x < stop; ++x) {
              auto suff = concatText.substr(SA[x], k);
              for (auto c : suff) {
                std::cerr << "*" << c << "*";
              }
              std::cerr << "\n";
            }
          }
        }
        currentKmer = nextKmer;
        start = stop;
      }
    } else {
      // If this isn't a valid suffix (contains a $)

      // If the previous interval was valid, put it
      // in the hash.
      if (currentKmer.length() == k and
          currentKmer.find_first_of('$') == std::string::npos) {
        mer = currentKmer.c_str();
        auto bits = mer.word(0);
        auto hashIt = khash.find(bits);
        if (hashIt == khash.end()) {
          if (start > 2) {
            if (concatText.substr(SA[start - 1], k) ==
                concatText.substr(SA[start], k)) {
              std::cerr << "T[SA[" << start - 1 << "]:" << k
                        << "] = " << concatText.substr(SA[start - 1], k)
                        << " = T[SA[" << start << "]:" << k << "]\n";
              std::cerr << "start = " << start << ", stop = " << stop << "\n";
              std::cerr << "[fatal (4)] THIS SHOULD NOT HAPPEN\n";
              std::exit(1);
            }
          }
          //auto lcpLength = findLCPLength(concatText,tlen,SA,start,stop);
          khash[bits] = {start, stop};
          all9mers[bits] = {start,stop};

          //khash[bits] = {start, stop};
          /*
          IndexT len = stop - start;
          bool overflow = (len >= std::numeric_limits<uint8_t>::max());
          uint8_t blen = overflow ? std::numeric_limits<uint8_t>::max() :
              static_cast<uint8_t>(len);
          khash[bits] = {start, blen};
          if (overflow) {
              overflowhash[bits] = len;
          }
          */
        } else {
          std::cerr << "\nERROR (2): trying to add same suffix " << currentKmer
                    << "multiple times!\n";
          auto prevInt = hashIt->second;
          std::cerr << "existing interval is [" << prevInt.begin() << ", "
                    << prevInt.end() << ")\n";
          for (auto x = prevInt.begin(); x < prevInt.end(); ++x) {
            std::cerr << concatText.substr(SA[x], k) << "\n";
          }
          std::cerr << "new interval is [" << start << ", " << stop << ")\n";
          for (auto x = start; x < stop; ++x) {
            std::cerr << concatText.substr(SA[x], k) << "\n";
          }
        }
      }
      // The current interval is invalid and empty
      currentKmer = nextKmer;
      start = stop;
    }
    if (stop % 1000000 == 0) {
      std::cerr << "\r\rprocessed " << stop << " positions";
    }
    // We always update the end position
    ++stop;
  }
  if (start < tlen) {
    if (currentKmer.length() == k and
        currentKmer.find_first_of('$') != std::string::npos) {
      mer = currentKmer.c_str();
      //auto lcpLength = findLCPLength(concatText,tlen,SA,start,stop);
      khash[mer.word(0)] = {start, stop};
      all9mers[mer.word(0)] = {start,stop};
      /*
      IndexT len = stop - start;
      bool overflow = (len >= std::numeric_limits<uint8_t>::max());
      uint8_t blen = overflow ? std::numeric_limits<uint8_t>::max() :
          static_cast<uint8_t>(len);
      khash[mer.get_bits(0, 2 * k)] = {start, blen};
      if (overflow) {
          overflowhash[mer.get_bits(0, 2 * k)] = len;
      }
      */
    }
  }
  std::cerr << "\nNumber of 9 mers " << all9mers.size() << "\n";
  std::ofstream all9mStream(outputDir + "all9mers.bin", std::ios::binary);
  {
    ScopedTimer timer;
    std::cerr << "saving 9mers to disk . . . ";
    //khash.serialize(typename spp_utils::pod_hash_serializer<WordT, rapmap::utils::SAInterval<IndexT>>(),
    //                &hashStream);
    cereal::BinaryOutputArchive all9mArchive(all9mStream);
    all9mArchive(all9mers);
    std::cerr << "done\n";
  }
  all9mStream.close();
  return true;

}

// IndexT is the index type.
// int32_t for "small" suffix arrays
// int64_t for "large" ones


template <typename IndexT, typename WordT, typename MerT>
bool buildHash(const std::string& outputDir, std::string& concatText,
               size_t tlen, uint32_t k, std::vector<IndexT>& SA,
               std::shared_ptr<spdlog::logger> log) {
  // Now, build the k-mer lookup table
    // The base type should always be uint64_t
    //using WordT = (k == 9)?(rapmap::utils::my_mer9::base_type):(rapmap::utils::my_mer::base_type);

  //I am done with building hash here
  //now I will test the safe
  //quick hack to make rankDir alive
  //very hacky need to change
  struct BitArrayDeleter {
    void operator()(BIT_ARRAY* b) {
      if(b != nullptr) {
        bit_array_free(b);
      }
    }
  };
  using BitArrayPointer = std::unique_ptr<BIT_ARRAY, BitArrayDeleter>;
  BitArrayPointer bitArray{nullptr};

  std::unique_ptr<rank9b> rankDict{nullptr} ;
  std::string rsFileName = outputDir + "rsd.bin";
  FILE* rsFile = fopen(rsFileName.c_str(), "r");
  {
    //logger->info("Loading Rank-Select Bit Array");
    bitArray.reset(bit_array_create(0));
    if (!bit_array_load(bitArray.get(), rsFile)) {
      //std::ce("Couldn't load bit array from {}!", rsFileName);
      std::exit(1);
    }
    //logger->info("There were {} set bits in the bit array", bit_array_num_bits_set(bitArray.get()));
    rankDict.reset(new rank9b(bitArray->words, bitArray->num_of_bits));
  }
  fclose(rsFile);
  //ends here


  spp::sparse_hash_map<rapmap::utils::TGroup, uint32_t, rapmap::utils::TGHasher> eqTable;

  auto getEqId = [&eqTable, &SA, &rankDict](IndexT start, IndexT stop) -> uint32_t {
    rapmap::utils::TGroup kset;
    for (size_t ival = start; ival < stop; ++ival) {
      kset.push_back(rankDict->rank(SA[ival]));
    }
    kset.prepare();
    uint32_t eqId{0};
    auto eqIt = eqTable.find(kset);
    if ( eqIt == eqTable.end()){
      eqId = eqTable.size();
      eqTable[kset] = eqId;//eqTable.size();
    } else {
      eqId = eqIt->second;
    }
    return eqId;
  };

  //using WordT = MerT::base_type ;
  MerT mer;
  RegHashT<WordT, rapmap::utils::kmerVal<IndexT>,
                         rapmap::utils::KmerKeyHasher> khash;
    //RegHashT<uint64_t, IndexT, rapmap::utils::KmerKeyHasher> overflowhash;
    //std::cerr << "sizeof(SAInterval<IndexT>) = " << sizeof(rapmap::utils::SAInterval<IndexT>) << '\n';
    //khash.set_empty_key(std::numeric_limits<uint64_t>::max());
  // The start and stop of the current interval
  IndexT start = 0, stop = 0;
  // An iterator to the beginning of the text
  auto textB = concatText.begin();
  auto textE = concatText.end();
  // The current k-mer as a string
  bool currentValid{false};
  std::string currentKmer;
  std::string nextKmer;
  while (stop < tlen) {
    // Check if the string starting at the
    // current position is valid (i.e. doesn't contain $)
    // and is <= k bases from the end of the string
    nextKmer = concatText.substr(SA[stop], k);
    if (nextKmer.length() == k and
        nextKmer.find_first_of('$') == std::string::npos) {
      // If this is a new k-mer, then hash the current k-mer
      if (nextKmer != currentKmer) {
        if (currentKmer.length() == k and
            currentKmer.find_first_of('$') == std::string::npos) {
          mer = currentKmer;
          auto bits = mer.word(0);
          auto hashIt = khash.find(bits);
          if (hashIt == khash.end()) {
            if (start > 1) {
              if (concatText.substr(SA[start - 1], k) ==
                  concatText.substr(SA[start], k)) {
                std::cerr << "T[SA[" << start - 1 << "]:" << k
                          << "] = " << concatText.substr(SA[start - 1], k)
                          << " = T[SA[" << start << "]:" << k << "]\n";
                std::cerr << "start = " << start << ", stop = " << stop << "\n";
                std::cerr << "[fatal (1)] THIS SHOULD NOT HAPPEN\n";
                std::exit(1);
              }
            }
            if (start == stop) {
              std::cerr << "[fatal (2)] Interval is empty! (start = " << start
                        << ") = (stop =  " << stop << ")\n";
            }
            if (start == stop) {
              std::cerr << "[fatal (3)] Interval is empty! (start = " << start
                        << ") = (stop =  " << stop << ")\n";
            }

            auto lcpLength = findLCPLength(concatText,tlen,SA,start,stop-1);
            auto eqId = getEqId(start, stop);
            khash[bits] = {start, stop, lcpLength, k, eqId};
            //khash[bits] = {start, stop, lcpLength};
            /*
            IndexT len = stop - start;
            bool overflow = (len >= std::numeric_limits<uint8_t>::max());
            uint8_t blen = overflow ? std::numeric_limits<uint8_t>::max() :
                static_cast<uint8_t>(len);
            khash[bits] = {start, blen};
            if (overflow) {
                overflowhash[bits] = len;
            }
            */
          } else {
            std::cerr << "\nERROR (1): trying to add same suffix "
                      << currentKmer << " (len = " << currentKmer.length()
                      << ") multiple times!\n";
            auto prevInt = hashIt->second.interval;
            std::cerr << "existing interval is [" << prevInt.begin() << ", "
                      << prevInt.end() << ")\n";
            for (auto x = prevInt.begin(); x < prevInt.end(); ++x) {
              auto suff = concatText.substr(SA[x], k);
              for (auto c : suff) {
                std::cerr << "*" << c << "*";
              }
              std::cerr << " (len = " << suff.length() << ")\n";
            }
            std::cerr << "new interval is [" << start << ", " << stop << ")\n";
            for (auto x = start; x < stop; ++x) {
              auto suff = concatText.substr(SA[x], k);
              for (auto c : suff) {
                std::cerr << "*" << c << "*";
              }
              std::cerr << "\n";
            }
          }
        }
        currentKmer = nextKmer;
        start = stop;
      }
    } else {
      // If this isn't a valid suffix (contains a $)

      // If the previous interval was valid, put it
      // in the hash.
      if (currentKmer.length() == k and
          currentKmer.find_first_of('$') == std::string::npos) {
        mer = currentKmer.c_str();
        auto bits = mer.word(0);
        auto hashIt = khash.find(bits);
        if (hashIt == khash.end()) {
          if (start > 2) {
            if (concatText.substr(SA[start - 1], k) ==
                concatText.substr(SA[start], k)) {
              std::cerr << "T[SA[" << start - 1 << "]:" << k
                        << "] = " << concatText.substr(SA[start - 1], k)
                        << " = T[SA[" << start << "]:" << k << "]\n";
              std::cerr << "start = " << start << ", stop = " << stop << "\n";
              std::cerr << "[fatal (4)] THIS SHOULD NOT HAPPEN\n";
              std::exit(1);
            }
          }
          auto lcpLength = findLCPLength(concatText,tlen,SA,start,stop-1);

          auto eqId = getEqId(start, stop);
          khash[bits] = {start, stop, lcpLength, k, eqId};
          //khash[bits] = {start, stop};
          /*
          IndexT len = stop - start;
          bool overflow = (len >= std::numeric_limits<uint8_t>::max());
          uint8_t blen = overflow ? std::numeric_limits<uint8_t>::max() :
              static_cast<uint8_t>(len);
          khash[bits] = {start, blen};
          if (overflow) {
              overflowhash[bits] = len;
          }
          */
        } else {
          std::cerr << "\nERROR (2): trying to add same suffix " << currentKmer
                    << "multiple times!\n";
          auto prevInt = hashIt->second.interval;
          std::cerr << "existing interval is [" << prevInt.begin() << ", "
                    << prevInt.end() << ")\n";
          for (auto x = prevInt.begin(); x < prevInt.end(); ++x) {
            std::cerr << concatText.substr(SA[x], k) << "\n";
          }
          std::cerr << "new interval is [" << start << ", " << stop << ")\n";
          for (auto x = start; x < stop; ++x) {
            std::cerr << concatText.substr(SA[x], k) << "\n";
          }
        }
      }
      // The current interval is invalid and empty
      currentKmer = nextKmer;

      start = stop;
    }
    if (stop % 1000000 == 0) {
      std::cerr << "\r\rprocessed " << stop << " positions";
    }
    // We always update the end position
    ++stop;
  }
  if (start < tlen) {
    if (currentKmer.length() == k and
        currentKmer.find_first_of('$') == std::string::npos) {

      mer = currentKmer.c_str();
      auto lcpLength = findLCPLength(concatText,tlen,SA,start,stop-1);
      auto eqId = getEqId(start, stop);
      khash[mer.word(0)] = {start, stop, lcpLength, k, eqId};
      /*
      IndexT len = stop - start;
      bool overflow = (len >= std::numeric_limits<uint8_t>::max());
      uint8_t blen = overflow ? std::numeric_limits<uint8_t>::max() :
          static_cast<uint8_t>(len);
      khash[mer.get_bits(0, 2 * k)] = {start, blen};
      if (overflow) {
          overflowhash[mer.get_bits(0, 2 * k)] = len;
      }
      */
    }
  }

  std::vector<rapmap::utils::TGroup> eqClasses(eqTable.size());

  std::cerr << "\n \n  EqTable size : " <<eqTable.size() << "\n\n" ;

  if(k != 9){
    for (auto it = eqTable.begin(); it != eqTable.end(); ++it) {
        //std::cerr << "\n WTF " << it->second << "\n" ;
        eqClasses[it->second] = it->first;
      }
  }

 if(k != 9){
   log->info("computing longest k-safe common prefixes");
   std::string lcpLogFilename = outputDir + "lcpLog.tsv";
   std::ofstream lcpLog(lcpLogFilename);
   lcpLog << "LCPLength\tSafeLCPLength\n";
   size_t numProc{0};
	 for(auto& it : khash){
		  auto& val = it.second ;
		  auto& interval = val.interval ;
		  updateSafe(concatText,tlen,SA,val,khash,rankDict,k, eqClasses);
      lcpLog << static_cast<uint32_t>(val.lcpLength) << '\t' << static_cast<uint32_t>(val.safeLength) << '\n';
      if (numProc % 100000 == 0) {
        std::cerr << "\r\rprocessed " << numProc << " of " << khash.size() << "intervals";
      }
      ++numProc;
	  }
   lcpLog.close();
 }

  std::cerr << "\nkhash had " << khash.size() << " keys\n";
  std::string fileName = (k == 9)?"hash9.bin":"hash.bin";
  std::ofstream hashStream(outputDir + fileName, std::ios::binary);
  {
    ScopedTimer timer;
    std::cerr << "saving hash to disk . . . ";
    //khash.serialize(typename spp_utils::pod_hash_serializer<WordT, rapmap::utils::SAInterval<IndexT>>(),
    //                &hashStream);
    khash.serialize(typename spp_utils::pod_hash_serializer<WordT, rapmap::utils::kmerVal<IndexT>>(),
                    &hashStream);
    std::cerr << "done\n";
  }
  hashStream.close();
  return true;
}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT> //, typename CoverageCalculator>
void indexTranscriptsSA(ParserT* parser,
                        std::string& outputDir,
                        bool noClipPolyA, bool usePerfectHash,
                        uint32_t numHashThreads,
                        std::string& sepStr,
                        std::mutex& iomutex,
                        std::shared_ptr<spdlog::logger> log) {
  // Create a random uniform distribution
  std::default_random_engine eng(271828);

  std::uniform_int_distribution<> dis(0, 3);

  // Hashers for getting txome signature
  picosha2::hash256_one_by_one seqHasher; seqHasher.init();
  picosha2::hash256_one_by_one nameHasher; nameHasher.init();

  uint32_t n{0};
  uint32_t k = rapmap::utils::my_mer::k();
  std::vector<std::string> transcriptNames;
  std::vector<int64_t> transcriptStarts;
  // std::vector<uint32_t> positionIDs;
  constexpr char bases[] = {'A', 'C', 'G', 'T'};
  uint32_t polyAClipLength{10};
  uint32_t numPolyAsClipped{0};
  uint32_t numNucleotidesReplaced{0};
  std::string polyA(polyAClipLength, 'A');

  using TranscriptList = std::vector<uint32_t>;
  using eager_iterator = MerMapT::array::eager_iterator;
  using KmerBinT = uint64_t;

  bool clipPolyA = !noClipPolyA;

  // http://biology.stackexchange.com/questions/21329/whats-the-longest-transcript-known
  // longest human transcript is Titin (108861), so this gives us a *lot* of
  // leeway before
  // we issue any warning.
  size_t tooLong = 200000;
  size_t numDistinctKmers{0};
  size_t numKmers{0};
  size_t currIndex{0};
  std::cerr << "\n[Step 1 of 4] : counting k-mers\n";

  // rsdic::RSDicBuilder rsdb;
  std::vector<uint64_t>
      onePos; // Positions in the bit array where we should write a '1'
  // remember the initial lengths (e.g., before clipping etc., of all transcripts)
  std::vector<uint32_t> completeLengths;
  // the stream of transcript sequence
  fmt::MemoryWriter txpSeqStream;

  std::string nonpolyAFile = outputDir + "/tra.fa" ;
  std::ofstream nonpolyAFileHandle(nonpolyAFile, std::ofstream::out | std::ofstream::app) ;


  {
    ScopedTimer timer;
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();

    while (parser->refill(rg)) {
      for (auto& read : rg) { // for each sequence
        std::string& readStr = read.seq;
        readStr.erase(
            std::remove_if(readStr.begin(), readStr.end(),
                           [](const char a) -> bool { return !(isprint(a)); }),
            readStr.end());

        seqHasher.process(readStr.begin(), readStr.end());

        uint32_t readLen = readStr.size();
        uint32_t completeLen = readLen;
        // First, replace non ATCG nucleotides
        for (size_t b = 0; b < readLen; ++b) {
          readStr[b] = ::toupper(readStr[b]);
          int c = jellyfish::mer_dna::code(readStr[b]);
          // Replace non-ACGT bases with pseudo-random bases
          if (jellyfish::mer_dna::not_dna(c)) {
            char rbase = bases[dis(eng)];
            c = jellyfish::mer_dna::code(rbase);
            readStr[b] = rbase;
            ++numNucleotidesReplaced;
          }
        }

        // Now, do Kallisto-esque clipping of polyA tails
        if (clipPolyA) {
          if (readStr.size() > polyAClipLength and
              readStr.substr(readStr.length() - polyAClipLength) == polyA) {

            auto newEndPos = readStr.find_last_not_of("Aa");
            // If it was all As
            if (newEndPos == std::string::npos) {
              log->warn("Entry with header [{}] appeared to be all A's; it "
                        "will be removed from the index!",
                        read.name);
              readStr.resize(0);
            } else {
              readStr.resize(newEndPos + 1);
            }
            ++numPolyAsClipped;
          }
        }

        readLen = readStr.size();
        // If the transcript was completely removed during clipping, don't
        // include it in the index.
        if (readStr.size() > 0) {
          // If we're suspicious the user has fed in a *genome* rather
          // than a transcriptome, say so here.
          if (readStr.size() >= tooLong) {
            log->warn("Entry with header [{}] was longer than {} nucleotides.  "
                      "Are you certain that "
                      "we are indexing a transcriptome and not a genome?",
                      read.name, tooLong);
          } else if (readStr.size() < k) {
            log->warn("Entry with header [{}], had length less than "
                      "the k-mer length of {} (perhaps after poly-A clipping)",
                      read.name, k);
          }

          uint32_t txpIndex = n++;

          // The name of the current transcript
          auto& recHeader = read.name;
          auto processedName = recHeader.substr(0, recHeader.find_first_of(sepStr));
          transcriptNames.emplace_back(processedName);
          nameHasher.process(processedName.begin(), processedName.end());

          // The position at which this transcript starts
          transcriptStarts.push_back(currIndex);
          // The un-molested length of this transcript
          completeLengths.push_back(completeLen);

          txpSeqStream << readStr;
          txpSeqStream << '$';
          currIndex += readLen + 1;
          onePos.push_back(currIndex - 1);

          //write transcript after polyA clipping
          nonpolyAFileHandle << ">" << processedName << "\n";
          nonpolyAFileHandle << readStr << "\n" ;

        } else {
            log->warn("Discarding entry with header [{}], since it had length 0 "
                      "(perhaps after poly-A clipping)",
                      read.name);
        }
      }
      if (n % 10000 == 0) {
        std::cerr << "\r\rcounted k-mers for " << n << " transcripts";
      }
    }
  }

  nonpolyAFileHandle.close() ;
  std::cerr << "\n";

  std::cerr << "Replaced " << numNucleotidesReplaced
            << " non-ATCG nucleotides\n";
  std::cerr << "Clipped poly-A tails from " << numPolyAsClipped
            << " transcripts\n";

  // Put the concatenated text in a string
  std::string concatText = txpSeqStream.str();
  // And clear the stream
  txpSeqStream.clear();

  // Build the suffix array
  size_t tlen = concatText.length();
  size_t maxInt = std::numeric_limits<int32_t>::max();
  bool largeIndex = (tlen + 1 > maxInt);

  // Make our dense bit arrray
  BIT_ARRAY* bitArray = bit_array_create(concatText.length());
  for (auto p : onePos) {
    bit_array_set_bit(bitArray, p);
  }

  onePos.clear();
  onePos.shrink_to_fit();

  std::string rsFileName = outputDir + "rsd.bin";
  FILE* rsFile = fopen(rsFileName.c_str(), "w");
  {
    ScopedTimer timer;
    std::cerr << "Building rank-select dictionary and saving to disk ";
    bit_array_save(bitArray, rsFile);
    std::cerr << "done\n";
  }
  fclose(rsFile);
  bit_array_free(bitArray);

  std::ofstream seqStream(outputDir + "txpInfo.bin", std::ios::binary);
  {
    ScopedTimer timer;
    std::cerr << "Writing sequence data to file . . . ";
    cereal::BinaryOutputArchive seqArchive(seqStream);
    seqArchive(transcriptNames);
    if (largeIndex) {
      seqArchive(transcriptStarts);
    } else {
      std::vector<int32_t> txpStarts(transcriptStarts.size(), 0);
      size_t numTranscriptStarts = transcriptStarts.size();
      for (size_t i = 0; i < numTranscriptStarts; ++i) {
        txpStarts[i] = static_cast<int32_t>(transcriptStarts[i]);
      }
      transcriptStarts.clear();
      transcriptStarts.shrink_to_fit();
      { seqArchive(txpStarts); }
    }
    // seqArchive(positionIDs);
    seqArchive(concatText);
    seqArchive(completeLengths);
    std::cerr << "done\n";
  }
  seqStream.close();


  // clear stuff we no longer need
  // positionIDs.clear();
  // positionIDs.shrink_to_fit();
  transcriptStarts.clear();
  transcriptStarts.shrink_to_fit();
  transcriptNames.clear();
  transcriptNames.shrink_to_fit();
  // done clearing


  if (largeIndex) {
    largeIndex = true;
    std::cerr << "[info] Building 64-bit suffix array "
                 "(length of generalized text is "
              << tlen << " )\n";
    using IndexT = int64_t;
    std::vector<IndexT> SA;
    bool success = buildSA(outputDir, concatText, tlen, SA);
    if (!success) {
      std::cerr << "[fatal] Could not build the suffix array!\n";
      std::exit(1);
    }


    if (usePerfectHash) {
      success = buildPerfectHash<IndexT>(outputDir, concatText, tlen, k, SA,
                                         numHashThreads, log);
    } else {
      success = buildHash<IndexT, rapmap::utils::my_mer::base_type, rapmap::utils::my_mer>(outputDir, concatText, tlen, k, SA, log);
    }
    if (!success) {
      std::cerr << "[fatal] Could not build the suffix interval hash!\n";
      std::exit(1);
    }
    std::cerr << "[info] building all 9 mer hash on the"
               "suffix array"
               <<"\n" ;
    bool success2 = buildHash<IndexT, rapmap::utils::my_mer9::base_type, rapmap::utils::my_mer9>(outputDir, concatText, tlen, 9, SA, log);
    if(!success2) {
        std::cerr << "[fatal] Could not build the 9 mer vector!\n";
        std::exit(1);
    }
  } else {
    std::cerr << "[info] Building 32-bit suffix array "
                 "(length of generalized text is "
              << tlen << ")\n";
    using IndexT = int32_t;
    std::vector<IndexT> SA;
    bool success = buildSA(outputDir, concatText, tlen, SA);
    if (!success) {
      std::cerr << "[fatal] Could not build the suffix array!\n";
      std::exit(1);
    }

    if (usePerfectHash) {
      success = buildPerfectHash<IndexT>(outputDir, concatText, tlen, k, SA,
                                         numHashThreads, log);
    } else {
      success = buildHash<IndexT, rapmap::utils::my_mer::base_type, rapmap::utils::my_mer>(outputDir, concatText, tlen, k, SA, log);
    }
    if (!success) {
      std::cerr << "[fatal] Could not build the suffix interval hash!\n";
      std::exit(1);
    }
    std::cerr << "[info] building all 9 mer hash on the"
               "suffix array"
               <<"\n" ;
    bool success2 = buildHash<IndexT, rapmap::utils::my_mer9::base_type, rapmap::utils::my_mer9>(outputDir, concatText, tlen, 9, SA, log);
    if(!success2) {
        std::cerr << "[fatal] Could not build the 9 mer vector!\n";
        std::exit(1);
    }
  }


  seqHasher.finish();
  nameHasher.finish();

  std::string indexVersion = "q5";
  IndexHeader header(IndexType::QUASI,
                     indexVersion,
                     true, k, largeIndex,
                     usePerfectHash);
  // Set the hash info
  std::string seqHash;
  std::string nameHash;
  picosha2::get_hash_hex_string(seqHasher, seqHash);
  picosha2::get_hash_hex_string(nameHasher, nameHash);
  header.setSeqHash(seqHash);
  header.setNameHash(nameHash);

  // Finally (since everything presumably succeeded) write the header
  std::ofstream headerStream(outputDir + "header.json");
  {
    cereal::JSONOutputArchive archive(headerStream);
    archive(header);
  }
  headerStream.close();
}

int rapMapSAIndex(int argc, char* argv[]) {
  std::cerr << "RapMap Indexer\n";

  TCLAP::CmdLine cmd("RapMap Indexer");
  TCLAP::ValueArg<std::string> transcripts("t", "transcripts",
                                           "The transcript file to be indexed",
                                           true, "", "path");
  TCLAP::ValueArg<std::string> index(
      "i", "index", "The location where the index should be written", true, "",
      "path");
  TCLAP::ValueArg<uint32_t> kval("k", "klen", "The length of k-mer to index",
                                 false, 31, "positive integer less than 32");

  TCLAP::ValueArg<std::string> customSeps("s", "headerSep", "Instead of a space or tab, break the header at the first "
                                          "occurrence of this string, and name the transcript as the token before "
                                          "the first separator", false, " \t", "string");
  TCLAP::SwitchArg noClip(
      "n", "noClip",
      "Don't clip poly-A tails from the ends of target sequences", false);
  TCLAP::SwitchArg perfectHash(
      "p", "perfectHash", "Use a perfect hash instead of sparse hash --- "
                          "somewhat slows construction, but uses less memory",
      false);
  /*
  TCLAP::SwitchArg perfectHash(
      "f", "frugalPerfectHash", "Use a frugal variant of the perfect hash --- "
                          "this will considerably slow construction, and somewhat slow lookup, but "
                          "hash construction and the subsequent mapping will require the least memory."
      false);
  */
  TCLAP::ValueArg<uint32_t> numHashThreads(
      "x", "numThreads",
      "Use this many threads to build the perfect hash function", false, 4,
      "positive integer <= # cores");
  cmd.add(transcripts);
  cmd.add(index);
  cmd.add(kval);
  cmd.add(noClip);
  cmd.add(perfectHash);
  cmd.add(customSeps);
  cmd.add(numHashThreads);
  cmd.parse(argc, argv);

  // stupid parsing for now
  std::string transcriptFile(transcripts.getValue());
  std::vector<std::string> transcriptFiles({transcriptFile});

  std::string sepStr = customSeps.getValue();

  uint32_t k = kval.getValue();
  if (k % 2 == 0) {
    std::cerr << "Error: k must be an odd value, you chose " << k << '\n';
    std::exit(1);
  } else if (k > 31) {
    std::cerr << "Error: k must not be larger than 31, you chose " << k << '\n';
    std::exit(1);
  }
  rapmap::utils::my_mer::k(k);
  rapmap::utils::my_mer9::k(9);

  std::string indexDir = index.getValue();
  if (indexDir.back() != '/') {
    indexDir += '/';
  }
  bool dirExists = rapmap::fs::DirExists(indexDir.c_str());
  bool dirIsFile = rapmap::fs::FileExists(indexDir.c_str());
  if (dirIsFile) {
    std::cerr << "The requested index directory already exists as a file.";
    std::exit(1);
  }
  if (!dirExists) {
    rapmap::fs::MakeDir(indexDir.c_str());
  }

  std::string logPath = indexDir + "quasi_index.log";
  auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_st>(logPath);
  auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_st>();
  auto consoleLog = spdlog::create("stderrLog", {consoleSink});
  auto fileLog = spdlog::create("fileLog", {fileSink});
  auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

  size_t numThreads{1};

  std::unique_ptr<single_parser> transcriptParserPtr{nullptr};

  size_t numProd = 1;
  transcriptParserPtr.reset(
			    new single_parser(transcriptFiles, numThreads, numProd));
  transcriptParserPtr->start();
  bool noClipPolyA = noClip.getValue();
  bool usePerfectHash = perfectHash.getValue();
  uint32_t numPerfectHashThreads = numHashThreads.getValue();
  std::mutex iomutex;
  indexTranscriptsSA(transcriptParserPtr.get(), indexDir, noClipPolyA,
                     usePerfectHash, numPerfectHashThreads, sepStr, iomutex, jointLog);

  // Output info about the reference
  std::ofstream refInfoStream(indexDir + "refInfo.json");
  {
    cereal::JSONOutputArchive archive(refInfoStream);
    archive(cereal::make_nvp("ReferenceFiles", transcriptFiles));
  }
  refInfoStream.close();

  return 0;
}
