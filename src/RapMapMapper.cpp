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
#include "jellyfish/hash_counter.hpp"

#include "tclap/CmdLine.h"

extern "C" {
#include "kseq.h"
}

#include "PairSequenceParser.hpp"
#include "RapMapUtils.hpp"
#include "RapMapIndex.hpp"
#include "RapMapFileSystem.hpp"
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
using KmerInfoList = std::vector<rapmap::utils::KmerInfo>;
using EqClassList = std::vector<rapmap::utils::EqClass>;
using EqClassLabelVec = std::vector<uint32_t>;

constexpr char bases[] = {'A', 'C', 'G', 'T'};

inline int32_t tid(uint64_t x) { return static_cast<uint32_t>(x >> 32); }
inline int32_t pos(uint64_t x) { return static_cast<uint32_t>(x); }

// Seed with a real random value, if available
std::random_device rd;

// Create a random uniform distribution
std::default_random_engine eng(rd());

std::uniform_int_distribution<> dis(0, 3);

struct QuasiAlignment {
    QuasiAlignment(uint32_t tidIn, uint32_t posIn,
                   bool fwdIn, uint32_t fragLenIn, bool isPairedIn = false) :
           tid(tidIn), pos(posIn), fwd(fwdIn),
           fragLen(fragLenIn), isPaired(isPairedIn) {}
    QuasiAlignment(QuasiAlignment&& other) = default;
    QuasiAlignment& operator=(const QuasiAlignment&) = default;
    uint32_t tid;
    uint32_t pos;
    bool fwd;
    uint32_t fragLen;
    bool isPaired;
};

// Walks the position list for this transcript and puts all hits
// on the back of the hits vector.
bool collectAllHits(uint32_t tid,
					uint32_t readLen,
                    bool hitRC,
					PositionList::iterator& posIt,
                    PositionList::iterator posEnd,
					std::vector<QuasiAlignment>& hits){
	bool foundHit{false};
	bool canAdvance = posIt < posEnd;
	bool nextTxp;
	bool isRC;
	int32_t pos;

	rapmap::utils::decodePosition(*posIt, pos, nextTxp, isRC);

	// The first position should always be a nextTxp, but we don't care
	nextTxp = false;

	while (canAdvance) {
        bool isReadRC = (isRC != hitRC);
		hits.emplace_back(tid, pos, isReadRC, readLen);
		foundHit = true;
		// If we can't advance the left but we need to, we're done
		if (!canAdvance) { return foundHit; }
		++posIt;
		rapmap::utils::decodePosition(*posIt, pos, nextTxp, isRC);
		canAdvance = !nextTxp;
	}
	return foundHit;
}


// Finds if there are positions within a specific transcript given by
// leftPosIt and rightPosIt within the distance constraints such that
// abs(leftPos - rightPos) < maxDist.  If so, the hit is appended to
// hits and the function returns true --- otherwise hits remains unchanged
// and the function returns false;
bool collectHitsWithPositionConstraint(uint32_t tid,
									   uint32_t readLen,
									   bool leftHitRC,
									   bool rightHitRC,
                                       uint32_t leftQueryPos,
                                       uint32_t rightQueryPos,
									   PositionList::iterator& leftPosIt,
									   PositionList::iterator& rightPosIt,
                                       PositionList::iterator posEnd,
									   uint32_t maxDist,
									   std::vector<QuasiAlignment>& hits){
	bool foundHit{false};
	bool canAdvance = true, canAdvanceLeft = leftPosIt < posEnd, canAdvanceRight = rightPosIt < posEnd;
	bool nextTxpLeft, nextTxpRight;
	bool isRCLeft, isRCRight;
	// True if the k-mer thinks the read is from fwd, false if from rc
	bool readStrandLeft, readStrandRight;
	int32_t leftPos, rightPos;

	rapmap::utils::decodePosition(*leftPosIt, leftPos, nextTxpLeft, isRCLeft);
	rapmap::utils::decodePosition(*rightPosIt, rightPos, nextTxpRight, isRCRight);

	//readStrandLeft = (leftHitRC == isRCLeft);
	//readStrandRight = (rightHitRC == isRCRight);

	// The first position should always be a nextTxp, but we don't care
	nextTxpLeft = false;
	nextTxpRight = false;
    //size_t i = 0;
	while (canAdvance) {
        //std::cerr << "advanced " << i++ << "\n";
		int32_t posDiff = rightPos - leftPos;
		uint32_t fragLen = std::abs(posDiff);
		// We found a hit (potentially -- what do we do about RCs here?)
		// I think we need to know if the k-mer from the *read* was fw or rc
		if (fragLen < maxDist) {
            bool isRC = (leftHitRC != isRCLeft);
            int32_t hitPos = (leftPos < rightPos) ? leftPos - leftQueryPos :
                                                 rightPos - rightQueryPos;
			hits.emplace_back(tid, hitPos, isRC, readLen);
			foundHit = true;
            break;
		}
		// rightPos >= leftPos (advance left)
		if (posDiff > 0) {
			// If we can't advance the left but we need to, we're done
			if (!canAdvanceLeft) { break; }
			++leftPosIt;
			rapmap::utils::decodePosition(*leftPosIt, leftPos, nextTxpLeft, isRCLeft);
			canAdvanceLeft = !nextTxpLeft and (leftPosIt < posEnd);
		} else if (posDiff < 0) { // leftPos > rightPos (advance right)
			// If we can't advance the right but we need to, we're done
			if (!canAdvanceRight) { break; }
			++rightPosIt;
			rapmap::utils::decodePosition(*rightPosIt, rightPos, nextTxpRight, isRCRight);
			canAdvanceRight = !nextTxpRight and (rightPosIt < posEnd);
		} else { // posDiff == 0 (advance both)
			// If we can't advance the left but we need to, we're done
			if (!canAdvanceLeft) { break; }
			++leftPosIt;
			rapmap::utils::decodePosition(*leftPosIt, leftPos, nextTxpLeft, isRCLeft);
			canAdvanceLeft = !nextTxpLeft and (leftPosIt < posEnd);
			// If we can't advance the right but we need to, we're done
			if (!canAdvanceRight) { break; }
			++rightPosIt;
			rapmap::utils::decodePosition(*rightPosIt, rightPos, nextTxpRight, isRCRight);
			canAdvanceRight = !nextTxpRight and (rightPosIt < posEnd);
		}

		// We can continue if we can advance either the left or right position
		canAdvance = (canAdvanceLeft or canAdvanceRight);
	}

    // Advance left and right until next txp
    while ( canAdvanceLeft ) {
        ++leftPosIt;
        rapmap::utils::decodePosition(*leftPosIt, leftPos, nextTxpLeft, isRCLeft);
        canAdvanceLeft = !nextTxpLeft and (leftPosIt < posEnd);
    }
    while ( canAdvanceRight ) {
        ++rightPosIt;
        rapmap::utils::decodePosition(*rightPosIt, rightPos, nextTxpRight, isRCRight);
        canAdvanceRight = !nextTxpRight and (rightPosIt < posEnd);
    }

	return foundHit;

}

void collectHits(RapMapIndex& rmi, std::string& readStr,
                 std::vector<QuasiAlignment>& hits) {

    /*
    auto& idx = rmi.idx;
    auto& tidList = rmi.tidList;
    auto& posList = rmi.posList;
    */

    auto jfhash = rmi.merHash.get();
    auto& kmerInfos = rmi.kmerInfos;
    auto& eqClasses = rmi.eqClassList;
    auto& eqClassLabels = rmi.eqLabelList;
    auto& posList = rmi.posList;
    auto posEnd = posList.end();

    rapmap::utils::my_mer mer;
    rapmap::utils::my_mer rcmer;
    auto k = rapmap::utils::my_mer::k();
    auto kbits = 2*k;
    auto readLen = readStr.length();
    uint32_t maxDist = static_cast<uint32_t>(readLen) * 1.25;
    size_t leftQueryPos = std::numeric_limits<size_t>::max();
    size_t rightQueryPos = std::numeric_limits<size_t>::max();
    bool leftHitRC = false, rightHitRC = false;

    auto endIt = kmerInfos.end();

    //IntervalIndex::iterator miniLeftHits = endIt;
    //IntervalIndex::iterator miniRightHits = endIt;

    KmerInfoList::iterator miniLeftHits = endIt;
    KmerInfoList::iterator miniRightHits = endIt;

    bool leftFwd{true};
    bool rightFwd{true};

    uint64_t merID;
    size_t kID;
    rapmap::utils::my_mer searchBuffer;

    for (size_t i = 0; i < readLen; ++i) {
        int c = jellyfish::mer_dna::code(readStr[i]);
        if (jellyfish::mer_dna::not_dna(c)) {
            c = jellyfish::mer_dna::code('G');
        }
        mer.shift_left(c);
        rcmer.shift_right(jellyfish::mer_dna::complement(c));
        if (i >= k) {
            auto& searchMer = (mer < rcmer) ? mer : rcmer;
            bool foundMer = jfhash->get_val_for_key(searchMer, &merID,
                                                    searchBuffer, &kID);
            if (foundMer) {
                miniLeftHits = kmerInfos.begin() + merID;
                leftHitRC = (searchMer == rcmer);
                leftQueryPos = i - k;
                break;
            }
        }
    }

    // found no hits in the entire read
    if (miniLeftHits == endIt) { return; }

    /* Super-fast & dirty --- just map the read based on the
	 * first k-mer hit
	 */
    //   else {
    //    auto leftIt = miniLeftHits->second.offset;
    //    auto leftLen = miniLeftHits->second.length;
    //    auto leftEnd = leftIt + leftLen;
    //    for (auto it = leftIt; it < leftEnd; ++it) {
    //        hits.push_back(tidList[it]);
    //    }
    //    return;
    //}

    // Now, start from the right and move left
    for (size_t i = readLen - 1; i >= leftQueryPos; --i) {
        int c = jellyfish::mer_dna::code(readStr[i]);
        if (jellyfish::mer_dna::not_dna(c)) {
            c = jellyfish::mer_dna::code('G');
        }
        mer.shift_right(c);
        rcmer.shift_left(jellyfish::mer_dna::complement(c));
        if (readLen - i >= k) {
            auto& searchMer = (mer < rcmer) ? mer : rcmer;
            bool foundMer = jfhash->get_val_for_key(searchMer, &merID,
                                                    searchBuffer, &kID);
            if (foundMer) {
                miniRightHits = kmerInfos.begin() + merID;
                rightHitRC = (searchMer == rcmer);
                // distance from the right end
                rightQueryPos = readLen - (i + k);
                break;
            }
        }
    }

	// Take the intersection of these two hit lists
	// Adapted from : http://en.cppreference.com/w/cpp/algorithm/set_intersection
	if (miniLeftHits != endIt) {
		// Equiv. class for left hit
		auto& eqClassLeft = eqClasses[miniLeftHits->eqId];
		// Iterator into, length of and end of the positon list
		auto leftPosIt = posList.begin() + miniLeftHits->offset;
		auto leftPosLen = miniLeftHits->count;
		auto leftPosEnd = leftPosIt + leftPosLen;
		// Iterator into, length of and end of the transcript list
		auto leftTxpIt = eqClassLabels.begin() + eqClassLeft.txpListStart;
		auto leftTxpListLen = eqClassLeft.txpListLen;
		auto leftTxpEnd = leftTxpIt + leftTxpListLen;

		if (miniRightHits != endIt) {
			// Equiv. class for right hit
			auto& eqClassRight = eqClasses[miniRightHits->eqId];
			// Iterator into, length of and end of the positon list
			auto rightPosIt = posList.begin() + miniRightHits->offset;
			auto rightPosLen = miniRightHits->count;
			auto rightPosEnd = rightPosIt + rightPosLen;
			// Iterator into, length of and end of the transcript list
			auto rightTxpIt = eqClassLabels.begin() + eqClassRight.txpListStart;
			auto rightTxpListLen = eqClassRight.txpListLen;
			auto rightTxpEnd = rightTxpIt + rightTxpListLen;

			//hits.resize(std::min(leftLen, rightLen));
			//size_t intSize = SIMDCompressionLib::SIMDintersection(&tidList[leftIt], leftLen,
			//                                                      &tidList[rightIt], rightLen,
			//                                                      &hits[0]);
			//hits.resize(intSize);


			hits.reserve(std::min(leftPosLen, rightPosLen));
			uint32_t leftTxp, rightTxp;
			while (leftTxpIt != leftTxpEnd and rightTxpIt != rightTxpEnd) {
				// Get the current transcript ID for the left and right eq class
				leftTxp = *leftTxpIt;
				rightTxp = *rightTxpIt;
				// If we need to advance the left txp, do it
				if (leftTxp < rightTxp) {
					++leftTxpIt;
				} else {
					// If the transcripts are equal (i.e. leftTxp >= rightTxp and !(rightTxp < leftTxp))
					// Then see if there are any hits here.
					if (!(rightTxp < leftTxp)) {
						collectHitsWithPositionConstraint(leftTxp, readLen,
                                                          leftHitRC, rightHitRC,
                                                          leftQueryPos, rightQueryPos,
                                                          leftPosIt, rightPosIt,
                                                          posEnd,
                                                          maxDist, hits);
						++leftTxpIt;
					}
					++rightTxpIt;
				}
			}
		} else { // If we had only hits from the left, then map this as an orphan
			hits.reserve(miniLeftHits->count);
			for (auto it = leftTxpIt; it < leftTxpEnd; ++it) {
				collectAllHits(*it, readLen, leftHitRC, leftPosIt, posEnd, hits);
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

        if (n % 1000 == 0) {
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
                                jointHits.emplace_back(leftTxp, startPos, leftIt->fwd, fragLen, true);
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
                auto& readName = j->data[i].first.header;
                sstream << "[AG]\t" << jointHits.size() << '\t'
                        <<  readName.substr(0, readName.length()-2) << '\n';
                for (const auto& qa : jointHits) {
                    sstream << txpNames[qa.tid] << '\t' << qa.pos << '\t' << qa.fwd
                            << '\t' << qa.fragLen << '\t'
                            << (qa.isPaired ? "Paired" : "Orphan") << '\n';
                }
            }

            if (n % 500000 == 0) {
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
    std::cerr << "RapMap Mapper\n";

    TCLAP::CmdLine cmd("RapMap Mapper");
    TCLAP::ValueArg<std::string> index("i", "index", "The location where the index should be written", true, "", "path");
    TCLAP::ValueArg<std::string> read1("1", "leftMates", "The location of the left paired-end reads", true, "", "path");
    TCLAP::ValueArg<std::string> read2("2", "rightMates", "The location of the right paired-end reads", true, "", "path");
    TCLAP::ValueArg<uint32_t> numThreads("t", "numThreads", "Number of threads to use", false, 1, "positive integer");
    cmd.add(index);
    cmd.add(read1);
    cmd.add(read2);
    cmd.add(numThreads);

    cmd.parse(argc, argv);


    // stupid parsing for now
    std::string indexPrefix(index.getValue());
    if (indexPrefix.back() != '/') {
        indexPrefix += "/";
    }

    if (!rapmap::fs::DirExists(indexPrefix.c_str())) {
        std::cerr << "It looks like the index you provided ["
                  << indexPrefix << "] doesn't exist\n";
        std::exit(1);
    }

    RapMapIndex rmi;
    rmi.load(indexPrefix);

    std::cerr << "\n\n\n\n";

    char* reads1 = const_cast<char*>(read1.getValue().c_str());
    char* reads2 = const_cast<char*>(read2.getValue().c_str());
    char* reads[] = {reads1, reads2};

    uint32_t nthread = numThreads.getValue();

    size_t maxReadGroup{1000}; // Number of reads in each "job"
    size_t concurrentFile{2}; // Number of files to read simultaneously
    std::unique_ptr<paired_parser> readParserPtr{nullptr};
    readParserPtr.reset(new paired_parser(4 * nthread, maxReadGroup,
                                         concurrentFile, reads, reads + 2));
    std::mutex iomutex;
    {
        ScopedTimer timer;
        std::cerr << "mapping reads . . . ";
        std::vector<std::thread> threads;
        for (size_t i = 0; i < nthread; ++i) {
            threads.emplace_back(processReads<paired_parser>, readParserPtr.get(), std::ref(rmi), std::ref(iomutex));
        }
        for (auto& t : threads) { t.join(); }
        //processReads(readParserPtr.get(), rmi, iomutex);
        //processReadsKSeq(seq1, seq2, idx, occs, iomutex);
        std::cerr << "done mapping reads\n";
    }
    return 0;
}
