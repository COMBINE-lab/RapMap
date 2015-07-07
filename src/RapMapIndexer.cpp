#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <fstream>

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "xxhash.h"
#include "btree/btree_map.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include "PairSequenceParser.hpp"
#include "RapMapUtils.hpp"

using paired_parser = pair_sequence_parser<char**>;
using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;

uint64_t encode(uint32_t tid, uint32_t pos) {
    uint64_t res = tid;
    res = res << 32;
    res |= pos;
    return res;
}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT>//, typename CoverageCalculator>
void processTranscripts(ParserT* parser,
                        std::mutex& iomutex) {
    // Seed with a real random value, if available
    std::random_device rd;

    // Create a random uniform distribution
    std::default_random_engine eng(rd());

    std::uniform_int_distribution<> dis(0, 3);

    uint32_t n{0};
    uint32_t k = rapmap::utils::my_mer::k();
    std::vector<std::string> transcriptNames;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    using TranscriptList = std::vector<uint32_t>;
    btree::btree_map<uint64_t, uint32_t> cmap;
    std::vector<std::string> transcriptSeqs;
    size_t numKmers{0};
  while(true) {
    typename ParserT::job j(*parser); 
    if(j.is_empty()) break;           

    for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
        std::string& readStr = j->data[i].seq;
        uint32_t readLen  = readStr.size();
        uint32_t txpIndex = n++;
        transcriptNames.push_back(j->data[i].header);

        rapmap::utils::my_mer mer;
        mer.polyT();
        for (size_t b = 0; b < readLen; ++b) {
            int c = jellyfish::mer_dna::code(readStr[b]);
            if (jellyfish::mer_dna::not_dna(c)) {
                char rbase = bases[dis(eng)];
                c = jellyfish::mer_dna::code(rbase);
                readStr[b] = rbase; 
            }
            mer.shift_left(c);
            if (b >= k) {
                auto canonicalMer = mer.get_canonical();
                auto key = canonicalMer.get_bits(0, 2*k);
                auto it = cmap.find(key);
                // If we found the k-mer, increment the count
                if (it != cmap.end()) {
                    it->second++;
                } else { // Otherwise, add it 
                    cmap[key]++;
                }
                // No matter what, our k-mer count increased
                numKmers++;
            }
        }
        transcriptSeqs.push_back(j->data[i].seq);
        if (n % 10000 == 0) {
            std::cerr << "counted k-mers for " << n << " transcripts\n";
        }
    }
  }
  
  auto numDistinctKmers = cmap.size();
  std::cerr << "parsed " << transcriptNames.size() << " transcripts\n";
  std::cerr << "There were " << numDistinctKmers << " distinct k-mers (canonicalized)\n";

  // now, we can prepare the vector for our "efficient map"
  uint32_t invalid = std::numeric_limits<uint32_t>::max();
  std::unordered_map<uint64_t, rapmap::utils::KmerInterval, rapmap::utils::KmerKeyHasher> hmap;
  std::vector<uint32_t> tidVec(numKmers, invalid);
  std::vector<uint32_t> posVec(numKmers, invalid);

  size_t tid{0};
  uint64_t intervalStart{0};
  for (auto& transcriptSeq : transcriptSeqs) {
      auto readLen = transcriptSeq.length();
      rapmap::utils::my_mer mer;
      mer.polyT();
      for (size_t b = 0; b < readLen; ++b) {
          int c = jellyfish::mer_dna::code(transcriptSeq[b]);
          mer.shift_left(c);
          if (b >= k) {
              bool isRC{false};
              auto canonicalMer = mer.get_canonical();
              auto key = canonicalMer.get_bits(0, 2*k);
              auto hIt = hmap.find(key);
              // if this is the first time we have seen this k-mer
              if (hIt == hmap.end()) {
                  isRC = (mer != canonicalMer);
                  auto occ = cmap[key];
                  hmap[key] = {intervalStart, occ};
                  intervalStart += occ;
                  cmap[key] = 0;
              }
              auto& interval = hmap[key];
              auto currentCount = cmap[key];
              tidVec[interval.offset + currentCount] = tid; 
              auto pos = b - k;
              // If this k-mer is from the rc of the txp, flip the high bit
              if (isRC) { pos |= 0x80000000; } 
              posVec[interval.offset + currentCount] = pos;
              cmap[key]++;
          }
      }
      if (tid % 10000 == 0) {
         std::cerr << "built index for " << tid << " transcripts\n";
      }
      ++tid;
  }
  
  std::ofstream txpStream("rapidx.txps.bin", std::ios::binary);
  {
    cereal::BinaryOutputArchive txpArchive(txpStream);
    txpArchive(transcriptNames);
  }
  txpStream.close();

  std::ofstream idxStream("rapidx.imap.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive indexArchive(idxStream);
      indexArchive(hmap);
  }
  idxStream.close();

  std::ofstream tidStream("rapidx.tid.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive tidArchive(tidStream);
      tidArchive(tidVec);
  }
  tidStream.close();

  std::ofstream posStream("rapidx.pos.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive posArchive(posStream);
      posArchive(posVec);
  }
  posStream.close();
}



int rapMapIndex(int argc, char* argv[]) {
    std::cerr << "indexing\n";
    // stupid parsing for now
    std::string reads1(argv[0]);
    std::vector<std::string> readFiles({ reads1 });

    int k = std::atoi(argv[1]);
    rapmap::utils::my_mer::k(static_cast<uint32_t>(k));

    size_t maxReadGroup{1000}; // Number of reads in each "job"
    size_t concurrentFile{2}; // Number of files to read simultaneously
    size_t numThreads{2};
    stream_manager streams(readFiles.begin(), readFiles.end(), concurrentFile);
    std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
    transcriptParserPtr.reset(new single_parser(4 * numThreads, maxReadGroup,
                              concurrentFile, streams));
    std::mutex iomutex;
    processTranscripts(transcriptParserPtr.get(), iomutex);
    return 0;
}


