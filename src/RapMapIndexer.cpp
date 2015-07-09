#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <fstream>

#include "tclap/CmdLine.h"

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

#include "RapMapUtils.hpp"
#include "RapMapFileSystem.hpp"
#include "ScopedTimer.hpp"

#include "jellyfish/file_header.hpp"
#include "jellyfish/binary_dumper.hpp"
#include "jellyfish/thread_exec.hpp"
#include "jellyfish/hash_counter.hpp"
#include "jellyfish/mer_overlap_sequence_parser.hpp"
#include "jellyfish/mer_iterator.hpp"
#include <chrono>

using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using MerMapT = jellyfish::cooperative::hash_counter<rapmap::utils::my_mer>;

uint64_t encode(uint32_t tid, uint32_t pos) {
    uint64_t res = tid;
    res = res << 32;
    res |= pos;
    return res;
}

constexpr uint32_t uint32HighestBitMask = 0x80000000;
constexpr uint32_t uint32LowBitMask = 0x7FFFFFFF;

constexpr uint32_t rcSetMask = 0x40000000;

// marks the second highest bit
void markRCBit(uint32_t& i) { i |= rcSetMask; }

// marks the highest bit
void markNewTxpBit(uint32_t& i) { i |= uint32HighestBitMask; }

bool wasSeen(uint32_t i) { return ((i & uint32HighestBitMask) >> 31) == 1; }

void markSeen(uint32_t& i) { i |= uint32HighestBitMask; }

uint32_t unmarked(uint32_t i) { return (i & uint32LowBitMask); }

class VectorHasher {
    public:
    size_t operator()(const std::vector<uint32_t>& vec) const {
        size_t hashVal{0};
        for (auto v : vec) {
            rapmap::utils::hashCombine(hashVal, v);
        }
        return hashVal;
    }
};

struct PosInfo {
    PosInfo(rapmap::utils::my_mer merIn, bool isRCIn, uint32_t posIn) :
        mer(merIn), isRC(isRCIn), pos(posIn) {}

    rapmap::utils::my_mer mer;
    bool isRC;
    uint32_t pos;
};

template <typename T>
void printVec(std::vector<T>& vec) {
    std::cerr << "{ ";
    for (auto& e : vec) {
        std::cerr << e << ", ";
    }
    std::cerr << "}";
}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT>//, typename CoverageCalculator>
void processTranscripts(ParserT* parser,
			std::string& outputDir,
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

    using KmerBinT = uint64_t;
    //create the hash
    size_t hashSize = 100000000;
    MerMapT merIntMap(hashSize, rapmap::utils::my_mer::k()*2, 32, 1, 126);
    std::vector<rapmap::utils::KmerInfo> kmerInfos;

    std::vector<std::string> transcriptSeqs;
    size_t numDistinctKmers{0};
    size_t numKmers{0};

    std::cerr << "\n[Step 1 of 4] : counting k-mers\n";

    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsedNS;
    {
        ScopedTimer timer;
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

                    uint64_t val;
                    start = std::chrono::high_resolution_clock::now();
                    auto found = merIntMap.ary()->get_val_for_key(canonicalMer, &val);
                    end = std::chrono::high_resolution_clock::now();
                    elapsedNS += end - start;
                    if (!found) {
                        merIntMap.add(canonicalMer, numDistinctKmers);
                        kmerInfos.emplace_back(txpIndex, 0, 1);
                        ++numDistinctKmers;
                    } else {
                        kmerInfos[val].count++;
                    }

                   /*
                    start = std::chrono::high_resolution_clock::now();
                    auto it = cmap.find(key);
                    end = std::chrono::high_resolution_clock::now();
                    elapsedNS += end - start;
                    // If we found the k-mer, increment the count
                    if (it != cmap.end()) {
                        it->second.count++;
                    } else { // Otherwise, add it
                        cmap[key] = KmerInfo(txpIndex, 0, 1);
                    }
                    */
                    // No matter what, our k-mer count increased
                    numKmers++;
                }
            }
            transcriptSeqs.push_back(j->data[i].seq);
            if (n % 10000 == 0) {
                std::cerr << "\r\rcounted k-mers for " << n << " transcripts";
            }
        }
    }
        std::cerr << "\n";
    }
    size_t vsize{0};
    using eager_iterator = MerMapT::array::eager_iterator;
    for (auto& kinfo : kmerInfos) {
        vsize += kinfo.count;
    }

    constexpr uint32_t uint32Invalid = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> transcriptIDs(vsize, uint32Invalid);

    std::cerr << "\n[Step 2 of 4] : marking k-mers\n";
    // Compute the equivalence classes for the k-mers
    uint32_t offset{0};
    uint32_t transcriptID{0};
    {
    ScopedTimer timer;
    for (auto& transcriptSeq : transcriptSeqs) {
        auto readLen = transcriptSeq.length();
        rapmap::utils::my_mer mer;
        mer.polyT();
        std::vector<KmerBinT> kmers;
        for (size_t b = 0; b < readLen; ++b) {
            int c = jellyfish::mer_dna::code(transcriptSeq[b]);
            mer.shift_left(c);
            if (b >= k) {
                auto canonicalMer = mer.get_canonical();

                uint64_t kmerIndex;
                auto found = merIntMap.ary()->get_val_for_key(canonicalMer, &kmerIndex);
                // Should ALWAYS find the key
                assert(found);

                auto& v = kmerInfos[kmerIndex];
                // use the highest bit to mark if we've seen this k-mer yet or not
                // If we haven't seen this k-mer yet
                if (!wasSeen(v.count)) {
                   // Where we start looking for transcripts for this k-mer
                   v.offset = offset;
                   offset += v.count;
                   // The number of transcripts we've currently added
                   v.count = 0;
                   markSeen(v.count);
                }

                // Note: We allow duplicate transcripts here --- they will always be adjacent
                auto lastOffset = v.offset + unmarked(v.count);
                transcriptIDs[lastOffset] = transcriptID;
                v.count++;
            }
        }
        if (transcriptID % 10000 == 0) {
            std::cerr << "\r\rmarked kmers for " << transcriptID << " transcripts";
        }
        ++transcriptID;
    }
    	std::cerr << "\n";
    }
    //printVec(transcriptIDs);
    // A hash to quickly and easily determine the equivalence classes
    std::unordered_map<std::vector<uint32_t>, uint32_t, VectorHasher> eqClassMap;
    // Holds the members of each equivalence class in the order in which classes
    // are assigned.  The final size should be \sum_{c \in eqclasses} |c|.
    std::vector<uint32_t> eqClassLabelVec;
    // Holds pointer information about a k-mer's equivalence class.
    // Specifically, where the label for the eq can be found
    // as an offset and length into eqClassTxpVec, and where the
    std::vector<rapmap::utils::EqClass> eqClasses;

    uint32_t eqClassVecSize{0};

    std::cerr << "\n[Step 3 of 4] : building k-mers equivalence classes\n";
    {
        ScopedTimer timer;
        const auto ary = merIntMap.ary();
        auto hashIt = ary->iterator_all<eager_iterator>();
        std::vector<uint32_t> tlist;
        while (hashIt.next()) {
            auto& key = hashIt.key();
            auto& val = kmerInfos[hashIt.val()];
            //auto& val = kv.second;
            auto offset = val.offset;
            auto num = unmarked(val.count);

	    tlist.clear();
            tlist.reserve(num);

            for (size_t idx = offset; idx < offset + num; ++idx) {
                auto tid = transcriptIDs[idx];
                // We won't consider duplicate transcript IDs when building the
                // equivalence classes
		//
		if (tlist.size() > 0 and tlist.back() > tid) {
		    std::cerr << "Non monotnoically increasing transcript id!\n";
		}
                if (tlist.size() == 0 or tlist.back() != tid) {
                    tlist.push_back(tid);
                }
            }

            auto eqIt = eqClassMap.find(tlist);
            uint32_t eqId{0};
            // If there is no such equivalence class yet, then add it
            if (eqIt == eqClassMap.end()) {
                eqId = eqClassMap.size();
                eqClassMap[tlist] = eqId;
                // The label of this eq-class starts at eqClassVecSize and is
                // tlist.size() elements long.
                eqClasses.emplace_back(eqClassVecSize, tlist.size());
                // Insert the label information into eqClassTxpVec
                eqClassLabelVec.insert(eqClassLabelVec.end(), tlist.begin(), tlist.end());
                eqClassVecSize += tlist.size();
            } else {
                eqId = eqIt->second;
            }
            // Set the equivalence class ID here for this transcript
            val.eqId = eqId;
            val.count = 0;
        }
    	std::cerr << "done! There were " << eqClassMap.size() << " classes\n";
    }
    // reuse the transcript IDs vector
    auto& posVec = transcriptIDs;


    std::cerr << "\n[Step 4 of 4] : finalizing index\n";
    transcriptID = 0;
    {
        ScopedTimer finalizeTimer;
        // Local vector to hold k-mers per transcript
        std::vector<PosInfo> posInfos;
        for (auto& transcriptSeq : transcriptSeqs) {
            posInfos.clear();
            auto readLen = transcriptSeq.length();
            rapmap::utils::my_mer mer;
            mer.polyT();
            for (size_t b = 0; b < readLen; ++b) {
                int c = jellyfish::mer_dna::code(transcriptSeq[b]);
                mer.shift_left(c);
                if (b >= k) {
                    auto canonicalMer = mer.get_canonical();
                    bool isRC = (mer != canonicalMer);
                    // Record the position of this k-mer in the transcript
                    auto pos = b - k;
                    posInfos.emplace_back(canonicalMer, isRC, pos);
                }
            }
            std::sort(posInfos.begin(), posInfos.end(),
                      [] (PosInfo& a, PosInfo& b) -> bool {
                        if (a.mer < b.mer) {
                            return true;
                        } else if (a.mer == b.mer) {
                            return a.pos < b.pos;
                        } else {
                            return false;
                        }
                      });

            PosInfo* prev{nullptr};
            for (auto& posInfo : posInfos) {
                    uint64_t kmerIndex;
                    auto found = merIntMap.ary()->get_val_for_key(
                                    posInfo.mer, &kmerIndex);
                    // Should ALWAYS find the key
                    assert(found);
                    auto& val = kmerInfos[kmerIndex];
                    auto offset = val.offset + val.count;
                    // Check if this offset is for a new transcript in the
                    // position list and, if so, set the proper flag bit.
                    if ( prev == nullptr or prev->mer != posInfo.mer) {
                        markNewTxpBit(posInfo.pos);
                    }
                    if ( posInfo.isRC ) {
                        markRCBit(posInfo.pos);
                    }

                    transcriptIDs[offset] = posInfo.pos;
                    ++val.count;
                }

            if (transcriptID % 10000 == 0) {
                std::cerr << "\r\rfinalized kmers for " << transcriptID << " transcripts";
            }
            ++transcriptID;
        }
	std::cerr << "\n";
    }

    // merIntMap --- jf hash
    // kmerInfos --- info for each k-mer
    // std::vector<uint32_t> eqClassLabelVec --- transcripts for each eq class
    // std::vector<EqClass> eqClasses --- where each eq class starts and ends
    // transcriptIDs --- position vec

    std::cerr << "Writing the index to " << outputDir << "\n";

    using JFFileHeader = jellyfish::file_header;
    using JFDumper = jellyfish::binary_dumper<MerMapT::array>;
    //typedef jellyfish::binary_reader<mer_dna, uint64_t> reader;
    //typedef jellyfish::binary_query_base<mer_dna, uint64_t> query;

    JFFileHeader fh;
    fh.fill_standard();
    fh.update_from_ary(*merIntMap.ary());
    JFDumper dumper(8*sizeof(uint32_t), k*2, 4, (outputDir + "rapidx.jfhash").c_str(), &fh);
    dumper.one_file(true);
    dumper.zero_array(false);
    dumper.dump(merIntMap.ary());

    std::ofstream kinfoStream(outputDir + "kinfo.bin", std::ios::binary);
    {
        cereal::BinaryOutputArchive kinfoArchive(kinfoStream);
        kinfoArchive(kmerInfos);
    }
    kinfoStream.close();

  std::ofstream eqLabelStream(outputDir + "eqlab.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive eqLabelArchive(eqLabelStream);
      eqLabelArchive(eqClassLabelVec);
  }
  eqLabelStream.close();

  std::ofstream eqClassStream(outputDir + "eqclass.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive eqClassArchive(eqClassStream);
      eqClassArchive(eqClasses);
  }
  eqClassStream.close();

  std::ofstream posStream(outputDir + "pos.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive posArchive(posStream);
      posArchive(transcriptIDs);
  }
  posStream.close();
  std::ofstream txpStream(outputDir + "txpnames.bin", std::ios::binary);
  {
    cereal::BinaryOutputArchive txpArchive(txpStream);
    txpArchive(transcriptNames);
  }
  txpStream.close();

  std::cerr << "transcriptIDs.size() = " << transcriptIDs.size() << "\n";
  std::cerr << "parsed " << transcriptNames.size() << " transcripts\n";
  std::cerr << "There were " << numDistinctKmers << " distinct k-mers (canonicalized)\n";

  // Data structure idea:
  // k-mer => equivalence class, position array offset
  // equivalence class = ordered (unique) list of transcripts
  // position array = *self-delimited* list of positions with the same order as txps in the equivalence class

    // using KmerBinT = uint64_t;

    /*
    // Compute the equivalence classes for the k-mers
    uint32_t transcriptID{0};
    for (auto& transcriptSeq : transcriptSeqs) {
        auto readLen = transcriptSeq.length();
        rapmap::utils::my_mer mer;
        mer.polyT();
        std::vector<KmerBinT> kmers;
        for (size_t b = 0; b < readLen; ++b) {
            int c = jellyfish::mer_dna::code(transcriptSeq[b]);
            mer.shift_left(c);
            if (b >= k) {
                bool isRC{false};
                auto canonicalMer = mer.get_canonical();
                auto key = canonicalMer.get_bits(0, 2*k);
                // push back the *ID* of this k-mer
                kmers.push_back(cmap[key]);
            }
        }
        partitioner.splitWith(kmers, tid);
        ++transcriptID;
    }
    // Compact the partition IDs
    partitioner.relabel();
    */

    /*
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
  */
}



int rapMapIndex(int argc, char* argv[]) {
    std::cerr << "RapMap Indexer\n";

    TCLAP::CmdLine cmd("RapMap Indexer");
    TCLAP::ValueArg<std::string> transcripts("t", "transcripts", "The transcript file to be indexed", true, "", "path");
    TCLAP::ValueArg<std::string> index("i", "index", "The location where the index should be written", true, "", "path");
    TCLAP::ValueArg<uint32_t> kval("k", "klen", "The length of k-mer to index", false, 31, "positive integer less than 32");
    cmd.add(transcripts);
    cmd.add(index);
    cmd.add(kval);

    cmd.parse(argc, argv);

    // stupid parsing for now
    std::string transcriptFile(transcripts.getValue());
    std::vector<std::string> transcriptFiles({ transcriptFile });

    uint32_t k = kval.getValue();
    rapmap::utils::my_mer::k(k);

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

    size_t maxReadGroup{1000}; // Number of reads in each "job"
    size_t concurrentFile{2}; // Number of files to read simultaneously
    size_t numThreads{2};
    stream_manager streams(transcriptFiles.begin(), transcriptFiles.end(), concurrentFile);
    std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
    transcriptParserPtr.reset(new single_parser(4 * numThreads, maxReadGroup,
                              concurrentFile, streams));
    std::mutex iomutex;
    processTranscripts(transcriptParserPtr.get(), indexDir, iomutex);
    return 0;
}


