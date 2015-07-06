#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <sstream>

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "SIMDCompressionAndIntersection/intersection.h"
#include "xxhash.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

extern "C" {
#include "kseq.h"
}

#include "PairSequenceParser.hpp"
#include "RapMapUtils.hpp"
#include "RapMapIndex.hpp"
#include "ScopedTimer.hpp"

// STEP 1: declare the type of file handler and the read() function  
KSEQ_INIT(int, read)

using paired_parser = pair_sequence_parser<char**>;
using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using TranscriptList = std::vector<uint32_t>;
using PositionList = std::vector<uint32_t>;
using KmerIndex = std::unordered_map<uint64_t, TranscriptList, rapmap::utils::KmerKeyHasher>;
using IntervalIndex = std::unordered_map<uint64_t, rapmap::utils::KmerInterval, rapmap::utils::KmerKeyHasher>;
using OccList = std::vector<uint64_t>;

constexpr char bases[] = {'A', 'C', 'G', 'T'};

inline int32_t tid(uint64_t x) { return static_cast<uint32_t>(x >> 32); }
inline int32_t pos(uint64_t x) { return static_cast<uint32_t>(x); }

// Seed with a real random value, if available
std::random_device rd;

// Create a random uniform distribution
std::default_random_engine eng(rd());

std::uniform_int_distribution<> dis(0, 3);

struct QuasiAlignment {
    QuasiAlignment(uint32_t tidIn, uint32_t posIn, bool fwdIn, uint32_t fragLenIn) :
           tid(tidIn), pos(posIn), fwd(fwdIn), fragLen(fragLenIn) {}
    QuasiAlignment(QuasiAlignment&& other) = default;
    QuasiAlignment& operator=(const QuasiAlignment&) = default;
    uint32_t tid;
    uint32_t pos;
    bool fwd;
    uint32_t fragLen;
};

void collectHits(RapMapIndex& rmi, std::string& readStr, std::vector<QuasiAlignment>& hits) {

    auto& idx = rmi.idx;
    auto& tidList = rmi.tidList;
    auto& posList = rmi.posList;

    rapmap::utils::my_mer mer;
    rapmap::utils::my_mer rcmer;
    auto k = rapmap::utils::my_mer::k();
    auto kbits = 2*k;
    auto readLen = readStr.length();
    uint32_t maxDist = static_cast<uint32_t>(readLen) * 1.25;
    auto endIt = idx.end();
    size_t leftQueryPos = std::numeric_limits<size_t>::max();
    size_t rightQueryPos = std::numeric_limits<size_t>::max();

    IntervalIndex::iterator miniLeftHits = endIt;
    IntervalIndex::iterator miniRightHits = endIt;

    bool leftFwd{true};
    bool rightFwd{true};

    for (size_t i = 0; i < readLen; ++i) {
        int c = jellyfish::mer_dna::code(readStr[i]);
        if (jellyfish::mer_dna::not_dna(c)) {
            c = jellyfish::mer_dna::code('G');
        }
        mer.shift_left(c);
        rcmer.shift_right(jellyfish::mer_dna::complement(c));
        if (i >= k) {
            auto merIt = idx.find(mer.get_bits(0, kbits));
            if (merIt != endIt) {
                miniLeftHits = merIt;
                leftQueryPos = i - k;
                break;
            } else {
                auto rcMerIt = idx.find(rcmer.get_bits(0, kbits));
                if (rcMerIt != endIt) {
                    miniLeftHits = rcMerIt;
                    leftQueryPos = i - k;
                    leftFwd = false;
                    break;
                }
            }
        } 
    }
    
    // found no hits in the entire read
    if (miniLeftHits == endIt) { return; } 
    /* Super-fast
       else {
        auto leftIt = miniLeftHits->second.offset;
        auto leftLen = miniLeftHits->second.length;
        auto leftEnd = leftIt + leftLen; 
        for (auto it = leftIt; it < leftEnd; ++it) {
            hits.push_back(tidList[it]);
        } 
        return;
    }
    */


    // Start from the right and move left
    for (size_t i = readLen - 1; i >= leftQueryPos; --i) {
        int c = jellyfish::mer_dna::code(readStr[i]);
        if (jellyfish::mer_dna::not_dna(c)) {
            c = jellyfish::mer_dna::code('G');
        }
        mer.shift_right(c);
        rcmer.shift_left(jellyfish::mer_dna::complement(c));
        if (readLen - i >= k) {
            auto merIt = idx.find(mer.get_bits(0, kbits));
            if (merIt != endIt) {
                miniRightHits = merIt;
                break;
            } else {
                auto rcMerIt = idx.find(rcmer.get_bits(0, kbits));
                if (rcMerIt != endIt) {
                    miniRightHits = rcMerIt;
                    rightFwd = false;
                    break;
                }
            }
        } 
    }

    // Take the intersection of these two hit lists
    // Adapted from : http://en.cppreference.com/w/cpp/algorithm/set_intersection
    if (miniLeftHits != endIt) {
        auto leftIt = miniLeftHits->second.offset;
        auto leftLen = miniLeftHits->second.length;
        auto leftEnd = leftIt + leftLen; 
        if (miniRightHits != endIt) {
            auto rightIt = miniRightHits->second.offset; 
            auto rightLen = miniRightHits->second.length; 
            auto rightEnd = rightIt + rightLen; 
            /*
            hits.resize(std::min(leftLen, rightLen));
            size_t intSize = SIMDCompressionLib::SIMDintersection(&tidList[leftIt], leftLen,
                                                                  &tidList[rightIt], rightLen,
                                                                  &hits[0]);
            hits.resize(intSize);                                          
            */
            hits.reserve(std::min(leftLen, rightLen));
            uint32_t leftTxp, rightTxp;
            while (leftIt != leftEnd && rightIt != rightEnd) {
                leftTxp = tidList[leftIt];
                rightTxp = tidList[rightIt];
                if (leftTxp < rightTxp) { 
                    ++leftIt;
                } else {
                    if (!(rightTxp < leftTxp)) { 
                        if (std::abs(static_cast<int32_t>(posList[leftIt]) - static_cast<int32_t>(posList[rightIt])) < maxDist) { 
                            hits.emplace_back(leftTxp, std::min(posList[leftIt], posList[rightIt]), leftFwd, 0);
                        }
                        ++leftIt;
                    }
                    ++rightIt;
                }
            }
        } else {
            hits.reserve(miniLeftHits->second.length);
            for (auto it = leftIt; it < leftEnd; ++it) {
                hits.emplace_back(tidList[it], posList[it], leftFwd, 0);
            } 
        }
    }
    
}

template <typename ParserT>//, typename CoverageCalculator>
void processReadsKSeq(ParserT* lseq,
                      ParserT* rseq,
                      RapMapIndex& rmi,
                      std::mutex& iomutex) {

    auto& txpNames = rmi.txpNames;
    uint32_t n{0};
    uint32_t k = rapmap::utils::my_mer::k();
    std::vector<std::string> transcriptNames;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    size_t batchSize{1000};
    std::vector<QuasiAlignment> leftHits;
    std::vector<QuasiAlignment> rightHits;
    std::vector<std::vector<QuasiAlignment>> jointHits(batchSize);
   
   int l1=0;
   int l2=0; 
   size_t peHits{0};
   size_t seHits{0};
   size_t readLen{0};
   while ( (l1 = kseq_read(lseq)) > 0 and (l2 = kseq_read(rseq)) > 0 ) {
       readLen = lseq->seq.l;
        n++;
        jointHits.clear();
        leftHits.clear();
        rightHits.clear();
        collectHits(rmi, lseq->seq.s, lseq->seq.l, leftHits);
        collectHits(rmi, rseq->seq.s, rseq->seq.l, rightHits);

        /*
        jointHits.resize(std::min(leftHits.size(), rightHits.size()));
        size_t intSize = SIMDCompressionLib::SIMDintersection(&leftHits[0], leftHits.size(),
                                                              &rightHits[0], rightHits.size(),
                                                              &jointHits[0]);
        jointHits.resize(intSize);                                          
        std::set_intersection(leftHits.begin(), leftHits.end(),
                              rightHits.begin(), rightHits.end(),
                              std::back_inserter(jointHits));
        */

        if (leftHits.size() > 0) {
            auto leftIt = leftHits.begin(); 
            auto leftEnd = leftHits.end(); 
            auto leftLen = std::distance(leftIt, leftEnd);
            if (rightHits.size() > 0) {
                auto rightIt = rightHits.begin(); 
                auto rightEnd = rightHits.end(); 
                auto rightLen = std::distance(rightIt, rightEnd);
                jointHits.reserve(std::min(leftLen, rightLen));
                while (leftIt != leftEnd && rightIt != rightEnd) {
                    uint32_t leftTxp = leftIt->tid; 
                    uint32_t rightTxp = rightIt->tid; 
                    if (leftTxp < rightTxp) { 
                        ++leftIt;
                    } else {
                        if (!(rightTxp < leftTxp)) { 
                            auto startPos = std::min(leftIt->pos, rightIt->pos);
                            auto endPos = std::max(leftIt->pos, rightIt->pos) + readLen;
                            jointHits.emplace_back(leftTxp, startPos, leftIt->fwd, static_cast<uint16_t>(endPos - startPos));
                            ++leftIt;
                        }
                        ++rightIt;
                    }
                }
            } 
        }

        if (jointHits.size() > 0) { 
           peHits += jointHits.size();
        } else if (leftHits.size() + rightHits.size() > 0) {
           seHits += leftHits.size() + rightHits.size();
        } 

        if (n % 1000000 == 0) {
            std::cerr << "saw " << n << " reads\n";
            std::cerr << "# pe hits per read = " << peHits / static_cast<float>(n) << "\n";
            std::cerr << "# se hits per read = " << seHits / static_cast<float>(n) << "\n";
        }
    }
  }



// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT>//, typename CoverageCalculator>
void processReads(ParserT* parser,
        RapMapIndex& rmi,
        std::mutex& iomutex) {

    auto& txpNames = rmi.txpNames;
    uint32_t n{0};
    uint32_t k = rapmap::utils::my_mer::k();
    std::vector<std::string> transcriptNames;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    std::stringstream sstream;
    size_t batchSize{1000};
    std::vector<QuasiAlignment> leftHits;
    std::vector<QuasiAlignment> rightHits;
    std::vector<QuasiAlignment> jointHits;

    size_t peHits{0};
    size_t seHits{0};
    size_t readLen{0};
    size_t maxNumHits{200};
    size_t nonDegenerateHits{0};
    while(true) {
        typename ParserT::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;           // If got nothing, quit
        for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
            readLen = j->data[i].first.seq.length();
            n++;
            jointHits.clear();
            leftHits.clear();
            rightHits.clear();
            collectHits(rmi, j->data[i].first.seq, leftHits);
            collectHits(rmi, j->data[i].second.seq, rightHits);  
            /* 
               std::set_intersection(leftHits.begin(), leftHits.end(),
               rightHits.begin(), rightHits.end(),
               std::back_inserter(jointHits));
               */

            if (leftHits.size() > 0) {
                auto leftIt = leftHits.begin(); 
                auto leftEnd = leftHits.end(); 
                auto leftLen = std::distance(leftIt, leftEnd);
                if (rightHits.size() > 0) {
                    auto rightIt = rightHits.begin(); 
                    auto rightEnd = rightHits.end(); 
                    auto rightLen = std::distance(rightIt, rightEnd);
                    size_t numHits{0};
                    jointHits.reserve(std::min(leftLen, rightLen));
                    while (leftIt != leftEnd && rightIt != rightEnd) {
                        uint32_t leftTxp = leftIt->tid; 
                        uint32_t rightTxp = rightIt->tid; 
                        if (leftTxp < rightTxp) { 
                            ++leftIt;
                        } else {
                            if (!(rightTxp < leftTxp)) { 
                                uint32_t startPos = std::min(leftIt->pos, rightIt->pos);
                                uint32_t endPos = std::max(leftIt->pos, rightIt->pos) + readLen;
                                uint32_t fragLen = endPos - startPos;
                                jointHits.emplace_back(leftTxp, startPos, leftIt->fwd, fragLen);
                                nonDegenerateHits += (fragLen > 0) ? 1 : 0; 
                                ++numHits;
                                if (numHits > maxNumHits) { break; }
                                ++leftIt;
                            }
                            ++rightIt;
                        }
                    }
                    if (numHits > maxNumHits) { jointHits.clear(); }
                } 
            }

            if (jointHits.size() > 0) { 
                peHits += jointHits.size();
            } else if (leftHits.size() + rightHits.size() > 0) {
                auto numHits = leftHits.size() + rightHits.size();
                seHits += numHits; 
                if (numHits < maxNumHits) {
                    jointHits.insert(jointHits.end(), 
                            std::make_move_iterator(leftHits.begin()),
                            std::make_move_iterator(leftHits.end()));
                    jointHits.insert(jointHits.end(), 
                            std::make_move_iterator(rightHits.begin()),
                            std::make_move_iterator(rightHits.end()));
                }
            } 

            if (jointHits.size() > 0) {
                sstream << "[AG]\t" << jointHits.size() << '\t' << j->data[i].first.header << '\n';
                for (const auto& qa : jointHits) {
                    sstream << txpNames[qa.tid] << '\t' << qa.pos << '\t' << qa.fwd << '\t' << qa.fragLen << '\n'; 
                }
            }

            if (n % 1000000 == 0) {
                if (n > 0) {
                    std::cerr << "\033[F\033[F\033[F";
                }
                std::cerr << "saw " << n << " reads\n";
                std::cerr << "# pe hits per read = " << peHits / static_cast<float>(n) << " (# non-degenerate = " << nonDegenerateHits << ")\n";
                std::cerr << "# se hits per read = " << seHits / static_cast<float>(n) << "\n";
            }
        } // for all reads in this job
        
        // DUMP OUTPUT
        iomutex.lock();
        std::cout << sstream.str();
        iomutex.unlock();
        sstream.str("");

    } // processed all reads
}



int rapMapMap(int argc, char* argv[]) {
    std::cerr << "mapping\n";
    // stupid parsing for now
    std::string indexPrefix(argv[0]);
    IntervalIndex idx;
    TranscriptList tidList;
    PositionList posList;
    std::vector<std::string> txpNames;
   
    RapMapIndex rmi;
    rmi.load(indexPrefix);

    std::cerr << "\n\n\n\n";

    int k = std::atoi(argv[1]);
    rapmap::utils::my_mer::k(static_cast<uint32_t>(k));

    char* reads1 = argv[2];
    char* reads2 = argv[3];
    char* reads[] = {reads1, reads2};

    /*
    FILE* fp1 = fopen(reads1, "r"); // STEP 2: open the file handler
    kseq_t *seq1 = kseq_init(fileno(fp1)); // STEP 3: initialize seq 
    FILE* fp2 = fopen(reads2, "r"); // STEP 2: open the file handler
    kseq_t *seq2 = kseq_init(fileno(fp2)); // STEP 3: initialize seq 
    */

    int nthread = std::atoi(argv[4]);

    size_t maxReadGroup{1000}; // Number of reads in each "job"
    size_t concurrentFile{2}; // Number of files to read simultaneously
    size_t numThreads{2};
    std::unique_ptr<paired_parser> readParserPtr{nullptr};
    readParserPtr.reset(new paired_parser(4 * numThreads, maxReadGroup,
                                         concurrentFile, reads, reads + 2));
    std::mutex iomutex;
    {
        ScopedTimer timer; 
        std::cerr << "mapping reads . . . ";
        std::vector<std::thread> threads;
        for (size_t i = 0; i < 10; ++i) { 
            threads.emplace_back(processReads<paired_parser>, readParserPtr.get(), std::ref(rmi), std::ref(iomutex));
        }
        for (auto& t : threads) { t.join(); }
        //processReads(readParserPtr.get(), rmi, iomutex);
        //processReadsKSeq(seq1, seq2, idx, occs, iomutex);
        std::cerr << "done mapping reads\n";
    }
    return 0;
}
