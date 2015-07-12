#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <sstream>

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "SIMDCompressionAndIntersection/intersection.h"
#include "xxhash.h"

#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/hash_counter.hpp"

#include "tclap/CmdLine.h"

/*extern "C" {
#include "kseq.h"
}
*/

#include "stringpiece.h"

#include "PairSequenceParser.hpp"
#include "RapMapUtils.hpp"
#include "RapMapIndex.hpp"
#include "RapMapFileSystem.hpp"
#include "RapMapConfig.hpp"
#include "ScopedTimer.hpp"
#include "SpinLock.hpp"

// STEP 1: declare the type of file handler and the read() function
// KSEQ_INIT(int, read)

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
#if defined __APPLE__
using SpinLockT = SpinLock;
#else
using SpinLockT = std::mutex;
#endif

constexpr char bases[] = {'A', 'C', 'G', 'T'};

inline int32_t tid(uint64_t x) { return static_cast<uint32_t>(x >> 32); }
inline int32_t pos(uint64_t x) { return static_cast<uint32_t>(x); }

// Seed with a real random value, if available
std::random_device rd;

// Create a random uniform distribution
std::default_random_engine eng(rd());

std::uniform_int_distribution<> dis(0, 3);

struct HitCounters {
    std::atomic<uint64_t> peHits{0};
    std::atomic<uint64_t> seHits{0};
    std::atomic<uint64_t> trueHits{0};
    std::atomic<uint64_t> totHits{0};
    std::atomic<uint64_t> numReads{0};
};

struct QuasiAlignment {
    QuasiAlignment(uint32_t tidIn, uint32_t posIn,
                   bool fwdIn, uint32_t readLenIn,
                   uint32_t fragLenIn = 0,
                   bool isPairedIn = false) :
           tid(tidIn), pos(posIn), fwd(fwdIn),
           readLen(readLenIn), fragLen(fragLenIn),
           isPaired(isPairedIn) {}
    QuasiAlignment(QuasiAlignment&& other) = default;
    QuasiAlignment& operator=(const QuasiAlignment&) = default;
    // Only 1 since the mate should have the same tid
    uint32_t tid;
    // Left-most position of the hit
    int32_t pos;
    // left-most position of the mate
    int32_t matePos;
    // Is the read from the forward strand
    bool fwd;
    // Is the mate from the forward strand
    bool mateIsFwd;
    // The fragment length (template length)
    uint32_t fragLen;
    // The read's length;
    uint32_t readLen;
    // The mate's length;
    uint32_t mateLen;
    // Is this a paired *alignment* or not
    bool isPaired;
    uint32_t mateStatus;
};

// from https://github.com/cppformat/cppformat/issues/105
class FixedBuffer : public fmt::Buffer<char> {
 public:
  FixedBuffer(char *array, std::size_t size)
    : fmt::Buffer<char>(array, size) {}

 protected:
  void grow(std::size_t size) {
    throw std::runtime_error("buffer overflow");
  }
};

class FixedWriter : public fmt::Writer {
 private:
  FixedBuffer buffer_;
 public:
  FixedWriter(char *array, std::size_t size)
    : fmt::Writer(buffer_), buffer_(array, size) {}
};

inline void adjustOverhang(int32_t& pos, uint32_t readLen,
                           uint32_t txpLen, FixedWriter& cigarStr) {
    cigarStr.clear();
    if (pos < 0) {
        int32_t clipLen = -pos;
        int32_t matchLen = readLen + pos;
        cigarStr.write("{}S{}M", clipLen, matchLen);
        // Now adjust the mapping position
        pos = 0;
    } else if (pos + readLen > txpLen) {
        int32_t matchLen = txpLen - pos;
        int32_t clipLen = pos + readLen - matchLen;
        cigarStr.write("{}M{}S", matchLen, clipLen);
    } else {
        cigarStr.write("{}M", readLen);
    }
}

inline void adjustOverhang(QuasiAlignment& qa, uint32_t txpLen,
                           FixedWriter& cigarStr1, FixedWriter& cigarStr2) {
    if (qa.isPaired) { // both mapped
        adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);
        adjustOverhang(qa.matePos, qa.mateLen, txpLen, cigarStr1);
    } else if (qa.mateStatus == 1) { // left read mapped
        adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);
    } else if (qa.mateStatus == 2) { // right read mapped
        adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr2);
    }
}


// get the sam flags for the quasialignment qaln.
// peinput is true if the read is paired in *sequencing*; false otherwise
// the sam flags for mate 1 are written into flags1 and for mate2 into flags2
inline void getSamFlags(const QuasiAlignment& qaln,
                        uint16_t& flags) {
    constexpr uint16_t pairedInSeq = 0x1;
    constexpr uint16_t properlyAligned = 0x2;
    constexpr uint16_t unmapped = 0x4;
    constexpr uint16_t mateUnmapped = 0x8;
    constexpr uint16_t isRC = 0x10;
    constexpr uint16_t mateIsRC = 0x20;
    constexpr uint16_t isRead1 = 0x40;
    constexpr uint16_t isRead2 = 0x80;
    constexpr uint16_t isSecondaryAlignment = 0x100;
    constexpr uint16_t failedQC = 0x200;
    constexpr uint16_t isPCRDup = 0x400;
    constexpr uint16_t supplementaryAln = 0x800;

    flags = 0;
    // Not paired in sequencing
    // flags1 = (peInput) ? pairedInSeq : 0;
    flags |= properlyAligned;
    // we don't output unmapped yet
    // flags |= unmapped
    // flags |= mateUnmapped
    flags |= (qaln.fwd) ? 0 : isRC;
    // Mate flag meaningless
    // flags1 |= (qaln.mateIsFwd) ? 0 : mateIsRC;
    flags |= isRead1;
    //flags2 |= isRead2;
}

// get the sam flags for the quasialignment qaln.
// peinput is true if the read is paired in *sequencing*; false otherwise
// the sam flags for mate 1 are written into flags1 and for mate2 into flags2
inline void getSamFlags(const QuasiAlignment& qaln,
                        bool peInput,
                        uint16_t& flags1,
                        uint16_t& flags2) {
    constexpr uint16_t pairedInSeq = 0x1;
    constexpr uint16_t properlyAligned = 0x2;
    constexpr uint16_t unmapped = 0x4;
    constexpr uint16_t mateUnmapped = 0x8;
    constexpr uint16_t isRC = 0x10;
    constexpr uint16_t mateIsRC = 0x20;
    constexpr uint16_t isRead1 = 0x40;
    constexpr uint16_t isRead2 = 0x80;
    constexpr uint16_t isSecondaryAlignment = 0x100;
    constexpr uint16_t failedQC = 0x200;
    constexpr uint16_t isPCRDup = 0x400;
    constexpr uint16_t supplementaryAln = 0x800;

    flags1 = flags2 = 0;
    flags1 = (peInput) ? pairedInSeq : 0;
    flags1 |= (qaln.isPaired) ? properlyAligned : 0;
    flags2 = flags1;
    // we don't output unmapped yet
    flags1 |= (qaln.mateStatus == 2) ? unmapped : 0;
    flags2 |= (qaln.mateStatus == 1) ? mateUnmapped : 0;
    flags1 |= (qaln.fwd) ? 0 : isRC;
    flags1 |= (qaln.mateIsFwd) ? 0 : mateIsRC;
    flags2 |= (qaln.mateIsFwd) ? 0 : isRC;
    flags2 |= (qaln.fwd) ? 0 : mateIsRC;
    flags1 |= isRead1;
    flags2 |= isRead2;
}


// Wraps the standard iterator of the Position list to provide
// some convenient functionality.  In the future, maybe this
// should be a proper iterator adaptor.
struct PositionListHelper{
	using PLIt = PositionList::iterator;

	PositionListHelper(PLIt itIn, PLIt endIn) :
		it_(itIn), end_(endIn) {}
        // The underlying iterator shouldn't be advanced further
	inline bool done() { return it_ == end_; }

	// The actual postion on the transcript
	int32_t pos() { return static_cast<int32_t>((*it_) & 0x3FFFFFFF); }

	// True if the position encoded was on the reverse complement strand
	// of the reference transcript, false otherwise.
	bool isRC() { return (*it_) & 0x40000000; }

	// True if we hit the position list for a new transcript, false otherwise
	bool isNewTxp() { return (*it_) & 0x80000000; }

	void advanceToNextTranscript() {
		if (it_ < end_) {
			do {
				++it_;
			} while (!isNewTxp() and it_ != end_);

		}
	}

	PLIt it_; // The underlying iterator
	PLIt end_; // The end of the container
};

// Walks the position list for this transcript and puts all hits
// on the back of the hits vector.
bool collectAllHits(uint32_t tid,
		uint32_t readLen,
		bool hitRC,
		PositionListHelper& ph,
		std::vector<QuasiAlignment>& hits,
        uint32_t mateStatus){
	bool foundHit{false};
	bool canAdvance = !ph.done();
	// The first position should always be a nextTxp, but we don't care
	bool nextTxp{false};
	bool isRC;
	int32_t pos;

	while (canAdvance) {
		isRC = ph.isRC();
		pos = ph.pos();
		bool isReadRC = (isRC != hitRC);
		hits.emplace_back(tid, pos, isReadRC, readLen);
        hits.back().mateStatus = mateStatus;
		foundHit = true;
		// If we can't advance the left but we need to, we're done
		if (!canAdvance) { return foundHit; }
		++(ph.it_);
		canAdvance = !ph.isNewTxp();
	}
    if (canAdvance) { ph.advanceToNextTranscript(); }
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
		PositionListHelper& leftPH,
		PositionListHelper& rightPH,
		uint32_t maxDist,
		std::vector<QuasiAlignment>& hits,
        uint32_t mateStatus){
	bool foundHit{false};
	bool canAdvance = true, canAdvanceLeft = !leftPH.done(), canAdvanceRight = !rightPH.done();
	// The first position should always be a nextTxp, but we don't care
	bool nextTxpLeft = false, nextTxpRight = false;
	bool isRCLeft, isRCRight;
	// True if the k-mer thinks the read is from fwd, false if from rc
	bool readStrandLeft, readStrandRight;
	int32_t leftPos, rightPos;
    std::vector<PositionListHelper> leftPosQueue, rightPosQueue;

#ifdef __DEBUG__
    if (!leftPH.isNewTxp()) {
        std::cerr << "leftPH = (" << leftPH.pos()
                  << ", " << leftPH.isNewTxp() << "), but should be start of "
                  "new txp list";
    }

    if (!rightPH.isNewTxp()) {
        std::cerr << "rightPH = (" << rightPH.pos()
                  << ", " << rightPH.isNewTxp() << "), but should be start of "
                  "new txp list";
    }
#endif // __DEBUG__

    while (canAdvance) {
		leftPos = leftPH.pos();
		rightPos = rightPH.pos();
		isRCLeft = leftPH.isRC();
		isRCRight = rightPH.isRC();

		int32_t posDiff = rightPos - leftPos;
		uint32_t fragLen = std::abs(posDiff);
		// We found a hit (potentially -- what do we do about RCs here?)
		// I think we need to know if the k-mer from the *read* was fw or rc
		if (fragLen < maxDist) {
			bool isRC = (leftHitRC != isRCLeft);
			int32_t hitPos = (leftPos < rightPos) ? leftPos - leftQueryPos :
								                    rightPos - rightQueryPos;
			hits.emplace_back(tid, hitPos, isRC, readLen);
            hits.back().mateStatus = mateStatus;
			foundHit = true;
			break;
		}
		// rightPos >= leftPos (advance left)
		if (posDiff > 0) {
			// If we can't advance the left but we need to, we're done
			if (!canAdvanceLeft) { break; }
			++(leftPH.it_);
            if (leftPH.isNewTxp() or leftPH.done()) {
                canAdvanceLeft = false;
                break;
            }
		} else if (posDiff < 0) { // leftPos > rightPos (advance right)
			// If we can't advance the right but we need to, we're done
			if (!canAdvanceRight) { break; }
			++(rightPH.it_);
            if (rightPH.isNewTxp() or rightPH.done()) {
                canAdvanceRight = false;
                break;
            }
		} else { // posDiff == 0 (advance both)
            /**
             * This is a strange case --- both k-mers (from different places)
             * in the read map to the same position.  First, this should
             * probably be a hit (i.e. < maxDist).  If not, is advancing
             * both the right thing to do?
             */
			// If we can't advance the left but we need to, we're done
			if (!canAdvanceLeft) { break; }
			++(leftPH.it_);
            if (leftPH.isNewTxp() or leftPH.done()) {
                canAdvanceLeft = false;
                break;
            }
            if (leftPH.isNewTxp()) { --(leftPH.it_); }
			// If we can't advance the right but we need to, we're done
			if (!canAdvanceRight) { break; }
			++(rightPH.it_);
            if (rightPH.isNewTxp() or rightPH.done()) {
                canAdvanceRight = false;
                break;
            }
		}

		// We can continue if we can advance either the left or right position
		canAdvance = (canAdvanceLeft or canAdvanceRight);
	}
    // Advance left and right until next txp
    if ( canAdvanceLeft ) {
        leftPH.advanceToNextTranscript();
    }
    if ( canAdvanceRight ) {
        rightPH.advanceToNextTranscript();
    }
    return foundHit;
}

// for mateStatus 0 is proper pair (paired),
// 1 is left of pair,
// 2 is right of pair,
// 3 is unpaired
void collectHits(RapMapIndex& rmi, std::string& readStr,
                 std::vector<QuasiAlignment>& hits,
                 uint32_t mateStatus) {

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
    uint32_t maxDist = static_cast<uint32_t>(readLen) * 1.5;
    size_t leftQueryPos = std::numeric_limits<size_t>::max();
    size_t rightQueryPos = std::numeric_limits<size_t>::max();
    bool leftHitRC = false, rightHitRC = false;

    auto endIt = kmerInfos.end();

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
		PositionListHelper leftPosHelper(leftPosIt, posList.end());
#ifdef __DEBUG__
        if (!leftPosHelper.isNewTxp()) {
            std::cerr << "\n Should definitely be new txp but "
                      << "leftPosHelper = ( "
                      << leftPosHelper.pos() << ", "
                      << leftPosHelper.isNewTxp() << ")\n";
        }
#endif
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
			PositionListHelper rightPosHelper(rightPosIt, posList.end());
#ifdef __DEBUG__
            if (!rightPosHelper.isNewTxp()) {
                std::cerr << "\n Should definitely be new txp but "
                    << "rightPosHelper = ( "
                    << rightPosHelper.pos() << ", "
                    << rightPosHelper.isNewTxp() << ")\n";
            }
#endif
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
                    // Advance to the next transcript in the
                    // equivalence class label
					++leftTxpIt;
					// Advance in the position array to the next ranscript
					leftPosHelper.advanceToNextTranscript();
				} else {
					// If the transcripts are equal (i.e. leftTxp >= rightTxp and !(rightTxp < leftTxp))
					// Then see if there are any hits here.
					if (!(rightTxp < leftTxp)) {
                        // If the hits are on the same transcript, look for
                        // a mapping position where they appear the appropriate
                        // distance apart.
                        // Note: The iterators into the *position* vector will
                        // be advanced, and should be at the start of the
                        // positions for the *next* transcript when this function
                        // returns.
						collectHitsWithPositionConstraint(leftTxp, readLen,
                              leftHitRC, rightHitRC,
                              leftQueryPos, rightQueryPos,
							  leftPosHelper, rightPosHelper,
                              maxDist, hits, mateStatus);
						++leftTxpIt;
						// advance pos
						// leftPosHelper.advanceToNextTranscript();
					} else {
                      // If the right transcript id was less than the left
                      // transcript id, then advance the right position
                      // iterator to the next transcript.
					  rightPosHelper.advanceToNextTranscript();
					}
                    // Advance the right transcript id regardless of whether
                    // we looked for a hit or not.
					++rightTxpIt;
				}
			}
		} else { // If we had only hits from the left, then map this as an orphan
			hits.reserve(miniLeftHits->count);
			for (auto it = leftTxpIt; it < leftTxpEnd; ++it) {
				collectAllHits(*it, readLen, leftHitRC, leftPosHelper, hits, mateStatus);
			}
		}
	}

}

void processReadsSingle(single_parser* parser,
        RapMapIndex& rmi,
        SpinLockT& iomutex,
        std::ostream& outStream,
        HitCounters& hctr) {

    auto& txpNames = rmi.txpNames;
    auto& txpLens = rmi.txpLens;
    uint32_t n{0};
    uint32_t k = rapmap::utils::my_mer::k();
    std::vector<std::string> transcriptNames;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    fmt::MemoryWriter sstream;
    size_t batchSize{1000};
    std::vector<QuasiAlignment> hits;

    size_t readLen{0};
    size_t maxNumHits{200};
    char cigarBuff[1000];
    FixedWriter cigarStr(cigarBuff, 1000);
    uint16_t flags;
    // 1000-bp reads are max here (get rid of hard limit later).
    std::string qualStr(1000, '~');
    while(true) {
        typename single_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;           // If got nothing, quit
        for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
            readLen = j->data[i].seq.length();
            ++hctr.numReads;
            hits.clear();
            collectHits(rmi, j->data[i].seq, hits, 3);
            /*
               std::set_intersection(leftHits.begin(), leftHits.end(),
               rightHits.begin(), rightHits.end(),
               std::back_inserter(jointHits));
               */
            auto numHits = hits.size();
            hctr.totHits += numHits;

            if (hits.size() > 0 and hits.size() < maxNumHits) {
                auto& readName = j->data[i].header;
#ifdef __DEBUG__
                auto before = readName.find_first_of(':');
                before = readName.find_first_of(':', before+1);
                auto after = readName.find_first_of(':', before+1);
                const auto& txpName = readName.substr(before+1, after-before-1);
#endif //__DEBUG__
                uint32_t alnCtr{0};
                for (auto& qa : hits) {
                    auto& transcriptName = txpNames[qa.tid];
                    // === SAM
                    getSamFlags(qa, flags);
                    if (alnCtr != 0) {
                        flags |= 0x900;
                    }

                    std::string* readSeq = &(j->data[i].seq);
                    std::string* qstr = &(j->data[i].qual);
                    adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], cigarStr);

                    sstream << readName << '\t' // QNAME
                        << flags << '\t' // FLAGS
                        << transcriptName << '\t' // RNAME
                        << qa.pos + 1 << '\t' // POS (1-based)
                        << 255 << '\t' // MAPQ
                        << cigarStr.c_str() << '\t' // CIGAR
                        << '*' << '\t' // MATE NAME
                        << 0 << '\t' // MATE POS
                        << qa.fragLen << '\t' // TLEN
                        << *readSeq << '\t' // SEQ
                        << *qstr << '\n';
                    ++alnCtr;
                    // === SAM
#ifdef __DEBUG__
                    if (txpNames[qa.tid] == txpName) { ++hctr.trueHits; }
#endif //__DEBUG__
                }
            }

            if (hctr.numReads % 1000000 == 0) {
                if (iomutex.try_lock()){
                    if (hctr.numReads > 0) {
#ifdef __DEBUG__
                        std::cerr << "\033[F\033[F\033[F";
#else
                        std::cerr << "\033[F\033[F";
#endif // __DEBUG__
                    }
                    std::cerr << "saw " << hctr.numReads << " reads\n";
                    std::cerr << "# hits per read = "
                        << hctr.totHits / static_cast<float>(hctr.numReads) << "\n";
#ifdef __DEBUG__
                    std::cerr << "The true hit was in the returned set of hits "
                        << 100.0 * (hctr.trueHits / static_cast<float>(hctr.numReads))
                        <<  "% of the time\n";
#endif // __DEBUG__
                    iomutex.unlock();
                }
            }
        } // for all reads in this job

        // DUMP OUTPUT
        iomutex.lock();
        outStream << sstream.str();
        iomutex.unlock();
        sstream.clear();

    } // processed all reads
}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
void processReadsPair(paired_parser* parser,
        RapMapIndex& rmi,
        SpinLockT& iomutex,
        std::ostream& outStream,
        HitCounters& hctr) {
    auto& txpNames = rmi.txpNames;
    std::vector<uint32_t>& txpLens = rmi.txpLens;
    uint32_t n{0};
    uint32_t k = rapmap::utils::my_mer::k();
    std::vector<std::string> transcriptNames;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    fmt::MemoryWriter sstream;
    size_t batchSize{1000};
    std::vector<QuasiAlignment> leftHits;
    std::vector<QuasiAlignment> rightHits;
    std::vector<QuasiAlignment> jointHits;

    size_t readLen{0};
    size_t maxNumHits{200};
    uint16_t flags1, flags2;
    // 1000-bp reads are max here (get rid of hard limit later).
    std::string qualStr(1000, '~');

    char buff1[1000];
    char buff2[1000];
    FixedWriter cigarStr1(buff1, 1000);
    FixedWriter cigarStr2(buff2, 1000);

    // 0 means properly aligned
    // 0x1 means only alignments for left read
    // 0x2 means only alignments for right read
    // 0x3 means "orphaned" alignments for left and right
    // (currently not treated as orphan).
    uint32_t orphanStatus{0};
    while(true) {
        typename paired_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
        if(j.is_empty()) break;           // If got nothing, quit
        for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
            readLen = j->data[i].first.seq.length();
            ++hctr.numReads;
            jointHits.clear();
            leftHits.clear();
            rightHits.clear();
            collectHits(rmi, j->data[i].first.seq, leftHits, 1);
        	collectHits(rmi, j->data[i].second.seq, rightHits, 2);
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
                                int32_t startRead1 = leftIt->pos;
                                int32_t startRead2 = rightIt->pos;
                                int32_t fragStartPos = std::min(leftIt->pos, rightIt->pos);
                                int32_t fragEndPos = std::max(leftIt->pos, rightIt->pos) + readLen;
                                uint32_t fragLen = fragEndPos - fragStartPos;
                                jointHits.emplace_back(leftTxp,
                                                       startRead1,
                                                       leftIt->fwd,
                                                       leftIt->readLen,
                                                       fragLen, true);
                                // Fill in the mate info
                                auto& qaln = jointHits.back();
                                qaln.mateLen = rightIt->readLen;
                                qaln.matePos = startRead2;
                                qaln.mateIsFwd = rightIt->fwd;
                                jointHits.back().mateStatus = 0;

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
                hctr.peHits += jointHits.size();
                orphanStatus = 0;
            } else if (leftHits.size() + rightHits.size() > 0) {
                auto numHits = leftHits.size() + rightHits.size();
                hctr.seHits += numHits;
                orphanStatus = 0;
                orphanStatus |= (leftHits.size() > 0) ? 0x1 : 0;
                orphanStatus |= (rightHits.size() > 0) ? 0x2 : 0;
                jointHits.insert(jointHits.end(),
                        std::make_move_iterator(leftHits.begin()),
                        std::make_move_iterator(leftHits.end()));
                jointHits.insert(jointHits.end(),
                        std::make_move_iterator(rightHits.begin()),
                        std::make_move_iterator(rightHits.end()));
            }

            if (jointHits.size() > 0) {
                auto& readName = j->data[i].first.header;
                auto& mateName = j->data[i].second.header;
                // trim /1 and /2 from pe read names
                if (readName.length() > 2 and
                        readName[readName.length() - 2] == '/') {
                    readName[readName.length() - 2] = '\0';
                }
                if (mateName.length() > 2 and
                        mateName[mateName.length() - 2] == '/') {
                    mateName[mateName.length() - 2] = '\0';
                }


#ifdef __DEBUG__
                auto before = readName.find_first_of(':');
                before = readName.find_first_of(':', before+1);
                auto after = readName.find_first_of(':', before+1);
                const auto& txpName = readName.substr(before+1, after-before-1);
#endif //__DEBUG__
                uint32_t alnCtr{0};
                for (auto& qa : jointHits) {
                    auto& transcriptName = txpNames[qa.tid];
                    // === SAM
                    if (qa.isPaired) {
                        getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            flags2 |= 0x900;
                        }
                        adjustOverhang(qa, txpLens[qa.tid], cigarStr1, cigarStr2);

                        sstream << readName.c_str() << '\t' // QNAME
                                << flags1 << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.pos + 1 << '\t' // POS (1-based)
                                << 255 << '\t' // MAPQ
                                << cigarStr1.c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << qa.matePos + 1 << '\t' // PNEXT
                                << qa.fragLen << '\t' // TLEN
                                << j->data[i].first.seq << '\t' // SEQ
                                << j->data[i].first.qual << '\n';

                        sstream << mateName.c_str() << '\t' // QNAME
                                << flags2 << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.matePos + 1 << '\t' // POS (1-based)
                                << 255 << '\t' // MAPQ
                                << cigarStr2.c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << qa.pos + 1 << '\t' // PNEXT
                                << qa.fragLen << '\t' // TLEN
                                << j->data[i].second.seq << '\t' // SEQ
                                << j->data[i].first.qual << '\n';
                    } else {
                        getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            if (qa.mateStatus == 1) {
                                flags2 |= 0x900;
                            } else {
                                flags1 |= 0x900;
                            }
                        }

                        std::string* readSeq{nullptr};
                        std::string* unalignedSeq{nullptr};

                        uint32_t flags, unalignedFlags;
                        std::string* qstr{nullptr};
                        std::string* unalignedQstr{nullptr};
                        std::string* unalignedName{nullptr};
                        FixedWriter* cigarStr;
                        if (qa.mateStatus == 1) { // left read
                            readName = j->data[i].first.header;
                            unalignedName = &j->data[i].second.header;

                            readSeq = &(j->data[i].first.seq);
                            unalignedSeq = &(j->data[i].second.seq);

                            qstr = &(j->data[i].first.qual);
                            unalignedQstr = &(j->data[i].second.qual);

                            flags = flags1;
                            unalignedFlags = flags2;

                            cigarStr = &cigarStr1;
                        } else { // right read
                            readName = j->data[i].second.header;
                            unalignedName = &(j->data[i].first.header);

                            readSeq = &(j->data[i].second.seq);
                            unalignedSeq = &(j->data[i].first.seq);

                            qstr = &(j->data[i].second.qual);
                            unalignedQstr = &(j->data[i].first.qual);

                            flags = flags2;
                            unalignedFlags = flags1;

                            cigarStr = &cigarStr2;
                        }

                        // If this is the first alignment of the group, then
                        // output the info for the unaligned mate.
                        if ( alnCtr == 0 and orphanStatus < 3) {
                            sstream << unalignedName->c_str() << '\t' // QNAME
                                << unalignedFlags << '\t' // FLAGS
                                << '*' << '\t' // RNAME
                                << 0 << '\t' // POS (1-based)
                                << 0 << '\t' // MAPQ
                                << readSeq->length() << 'M' << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << 0 << '\t' // PNEXT (only 1 read in template)
                                << 0 << '\t' // TLEN (spec says 0, not read len)
                                << *unalignedSeq << '\t' // SEQ
                                << *unalignedQstr << '\n';
                        }

                        adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], *cigarStr);
                        sstream << readName.c_str() << '\t' // QNAME
                                << flags << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.pos + 1 << '\t' // POS (1-based)
                                << 255 << '\t' // MAPQ
                                << cigarStr->c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << 0 << '\t' // PNEXT (only 1 read in templte)
                                << 0 << '\t' // TLEN (spec says 0, not read len)
                                << *readSeq << '\t' // SEQ
                                << *qstr << '\n';
                    }
                    ++alnCtr;
                    // === SAM
                    /*
                       sstream << "[AG]\t" << jointHits.size() << '\t'
                       <<  readName.substr(0, readName.length()-2) << '\n';

                    sstream << txpNames[qa.tid] << '\t' << qa.pos << '\t' << qa.fwd
                        << '\t' << qa.fragLen << '\t'
                        << (qa.isPaired ? "Paired" : "Orphan") << '\n';
                    */
#ifdef __DEBUG__
                    if (txpNames[qa.tid] == txpName) { ++hctr.trueHits; }
#endif //__DEBUG__
                }
            }

            if (hctr.numReads % 1000000 == 0) {
                if (iomutex.try_lock()) {
                    if (hctr.numReads > 0) {
#ifdef __DEBUG__
                        std::cerr << "\033[F\033[F\033[F\033[F";
#else
                        std::cerr << "\033[F\033[F\033[F";
#endif // __DEBUG__
                    }
                    std::cerr << "saw " << hctr.numReads << " reads\n";
                    std::cerr << "# pe hits per read = "
                        << hctr.peHits / static_cast<float>(hctr.numReads) << "\n";
                    std::cerr << "# se hits per read = "
                        << hctr.seHits / static_cast<float>(hctr.numReads) << "\n";
#ifdef __DEBUG__
                    std::cerr << "The true hit was in the returned set of hits "
                        << 100.0 * (hctr.trueHits / static_cast<float>(hctr.numReads))
                        <<  "% of the time\n";
#endif // __DEBUG__
                    iomutex.unlock();
                }
            }
        } // for all reads in this job

        // DUMP OUTPUT
        iomutex.lock();
        outStream << sstream.str();
        iomutex.unlock();
        sstream.clear();

    } // processed all reads
}

void writeSAMHeader(RapMapIndex& rmi, std::ostream& outStream) {
    fmt::MemoryWriter hd;
    hd.write("@HD\tVN:0.1\tSO:unknown\n");

    auto& txpNames = rmi.txpNames;
    auto& txpLens = rmi.txpLens;

    auto numRef = txpNames.size();
    for (size_t i = 0; i < numRef; ++i) {
        hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
    }
    // Eventuall output a @PG line
    //hd.format("@PG\t");
    outStream << hd.str();
}

// from http://stackoverflow.com/questions/9435385/split-a-string-using-c11
std::vector<std::string> tokenize(const std::string &s, char delim) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> elems;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

int rapMapMap(int argc, char* argv[]) {
    std::cerr << "RapMap Mapper\n";

    std::string versionString = rapmap::version;
    TCLAP::CmdLine cmd(
            "RapMap Mapper",
            ' ',
            versionString);
    cmd.getProgramName() = "rapmap";

    TCLAP::ValueArg<std::string> index("i", "index", "The location where the index should be written", true, "", "path");
    TCLAP::ValueArg<std::string> read1("1", "leftMates", "The location of the left paired-end reads", true, "", "path");
    TCLAP::ValueArg<std::string> read2("2", "rightMates", "The location of the right paired-end reads", false, "", "path");
    TCLAP::ValueArg<std::string> unmatedReads("r", "unmatedReads", "The location of single-end reads", true, "", "path");
    TCLAP::ValueArg<uint32_t> numThreads("t", "numThreads", "Number of threads to use", false, 1, "positive integer");
    TCLAP::ValueArg<std::string> outname("o", "output", "The output file (default: stdout)", false, "", "path");
    cmd.add(index);

    std::vector<TCLAP::Arg*> xorList({&read1, &unmatedReads});
    cmd.xorAdd(xorList);
    cmd.add(read2);
    cmd.add(outname);
    cmd.add(numThreads);
    cmd.parse(argc, argv);

    bool pairedEnd = (read1.isSet() or read2.isSet());
    if (pairedEnd and (read1.isSet() != read2.isSet())) {
        std::cerr << "You must set both the -1 and -2 arguments to align "
                  << "paired end reads!\n";
        std::exit(1);
    }

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

    std::cerr << "\n\n";

    // from: http://stackoverflow.com/questions/366955/obtain-a-stdostream-either-from-stdcout-or-stdofstreamfile
    // set either a file or cout as the output stream
    std::streambuf* outBuf;
    std::ofstream outFile;
    bool haveOutputFile{false};
    if (outname.getValue() == "") {
        outBuf = std::cout.rdbuf();
    } else {
        outFile.open(outname.getValue());
        outBuf = outFile.rdbuf();
        haveOutputFile = true;
    }
    // Now set the output stream to the buffer, which is
    // either std::cout, or a file.
    std::ostream outStream(outBuf);

    uint32_t nthread = numThreads.getValue();
    std::unique_ptr<paired_parser> pairParserPtr{nullptr};
    std::unique_ptr<single_parser> singleParserPtr{nullptr};

    writeSAMHeader(rmi, outStream);

    SpinLockT iomutex;
    {
        ScopedTimer timer;
        HitCounters hctrs;
        std::cerr << "mapping reads . . . ";
        if (pairedEnd) {
            std::vector<std::thread> threads;
            std::vector<std::string> read1Vec = tokenize(read1.getValue(), ',');
            std::vector<std::string> read2Vec = tokenize(read2.getValue(), ',');

            if (read1Vec.size() != read2Vec.size()) {
                std::cerr << "The number of provided files for -1 and -2 must "
                    << "be the same!\n";
                std::exit(1);
            }

            size_t numFiles = read1Vec.size() + read2Vec.size();
            char** pairFileList = new char*[numFiles];
            for (size_t i = 0; i < read1Vec.size(); ++i) {
                pairFileList[2*i] = const_cast<char*>(read1Vec[i].c_str());
                pairFileList[2*i+1] = const_cast<char*>(read2Vec[i].c_str());
            }
            size_t maxReadGroup{1000}; // Number of reads in each "job"
            size_t concurrentFile{2}; // Number of files to read simultaneously
            pairParserPtr.reset(new paired_parser(4 * nthread, maxReadGroup,
                        concurrentFile,
                        pairFileList, pairFileList+numFiles));

            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsPair,
                        pairParserPtr.get(),
                        std::ref(rmi),
                        std::ref(iomutex),
                        std::ref(outStream),
                        std::ref(hctrs));
            }

            for (auto& t : threads) { t.join(); }
            delete [] pairFileList;
        } else {
            std::vector<std::thread> threads;
            std::vector<std::string> unmatedReadVec = tokenize(unmatedReads.getValue(), ',');
            size_t maxReadGroup{1000}; // Number of reads in each "job"
            size_t concurrentFile{1};
            stream_manager streams( unmatedReadVec.begin(), unmatedReadVec.end(),
                    concurrentFile);
            singleParserPtr.reset(new single_parser(4 * nthread,
                        maxReadGroup,
                        concurrentFile,
                        streams));

            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsSingle,
                        singleParserPtr.get(),
                        std::ref(rmi),
                        std::ref(iomutex),
                        std::ref(outStream),
                        std::ref(hctrs));
            }
            for (auto& t : threads) { t.join(); }
        }
        std::cerr << "done mapping reads\n";
    }

    if (haveOutputFile) {
        outFile.close();
    }
    return 0;
}


/*
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

//        jointHits.resize(std::min(leftHits.size(), rightHits.size()));
//        size_t intSize = SIMDCompressionLib::SIMDintersection(&leftHits[0], leftHits.size(),
//                                                              &rightHits[0], rightHits.size(),
//                                                              &jointHits[0]);
//        jointHits.resize(intSize);
//        std::set_intersection(leftHits.begin(), leftHits.end(),
//                              rightHits.begin(), rightHits.end(),
//                              std::back_inserter(jointHits));
//
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
*/
