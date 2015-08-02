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
#include <tuple>
#include <sstream>

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "HitManager.hpp"
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

// #define __DEBUG__
// #define __TRACK_CORRECT__


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
using PositionListHelper = rapmap::utils::PositionListHelper;
#if defined __APPLE__
using SpinLockT = SpinLock;
#else
using SpinLockT = std::mutex;
#endif

// "Fake" mutex for single-threaded exceuction that does nothing;
class NullMutex {
    public:
	void lock() { return; }
	bool try_lock() { return true; }
	void unlock() { return; }
};



constexpr char bases[] = {'A', 'C', 'G', 'T'};

inline int32_t tid(uint64_t x) { return static_cast<uint32_t>(x >> 32); }
inline int32_t pos(uint64_t x) { return static_cast<uint32_t>(x); }

// Seed with a real random value, if available
std::random_device rd;

// Create a random uniform distribution
std::default_random_engine eng(rd());

std::uniform_int_distribution<> dis(0, 3);

using HitCounters = rapmap::utils::HitCounters;
using MateStatus = rapmap::utils::MateStatus;
using HitInfo = rapmap::utils::HitInfo;
using ProcessedHit = rapmap::utils::ProcessedHit;
using QuasiAlignment = rapmap::utils::QuasiAlignment;
using FixedWriter = rapmap::utils::FixedWriter;




// Walks the position list for this transcript and puts all hits
// on the back of the hits vector.
bool collectAllHits(uint32_t tid,
		uint32_t readLen,
		bool hitRC,
		PositionListHelper& ph,
		std::vector<QuasiAlignment>& hits,
        MateStatus mateStatus){
	bool foundHit{false};
	bool canAdvance = !ph.done();
	// The first position should always be a nextTxp, but we don't care
	bool nextTxp{false};
	bool isRC;
	int32_t pos;

//	while (canAdvance) {
    // only return the first hit for now
    if (canAdvance) {
		isRC = ph.isRC();
		pos = ph.pos();
		bool isReadRC = (isRC != hitRC);
		hits.emplace_back(tid, pos, isReadRC, readLen);
        hits.back().mateStatus = mateStatus;
		foundHit = true;
		// If we can't advance the left but we need to, we're done
		/*if (!canAdvance) { return foundHit; }
		++(ph.it_);
		canAdvance = !ph.isNewTxp();
		*/
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
		MateStatus mateStatus){
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


class SkippingKmerSearcher{
    private:
	std::string* qstr;
	const char* qCharArray;
	uint32_t qlen;
	uint32_t k;
	uint32_t klen;
	uint32_t startPos;
	uint32_t nextBaseIndex;
	rapmap::utils::my_mer mer;
	rapmap::utils::my_mer rcmer;
	rapmap::utils::my_mer tempMer;
	static constexpr uint32_t invalidIndex = std::numeric_limits<uint32_t>::max();

    public:
	SkippingKmerSearcher(std::string& queryStr) :
		qstr(&queryStr),
		qCharArray(queryStr.c_str()),
		qlen(queryStr.length()),
		k(rapmap::utils::my_mer::k()),
		klen(0),
		startPos(0),
		nextBaseIndex(0) {
		    next();
	}


	// return the index of the current k-mer (start) in the query string
	uint32_t queryIndex() { return startPos; }

	rapmap::utils::my_mer getMer(bool& isRC) {
	    tempMer = mer.get_canonical();
	    isRC = (mer != tempMer);
	    return tempMer;
	}

	bool backwardSearch(uint32_t searchPos) {
	    // Perform a backward search starting from searchPos

	    // If we're not at least k bases in, we can't do a
	    // backward search
	    if (searchPos < k) {
	    	return false;
	    }

	    // We can't search off the end
	    if (searchPos + k > qlen) {
		searchPos = qlen - k;
	    }


	    // otherwise start a new k-mer at the jump position
	    while (!tempMer.from_chars(&qCharArray[searchPos])) {
		// If it wasn't a valid k-mer, find the offending base
		// and try to start k bases before it
		uint32_t invalidLoc = qstr->find_last_not_of("aAcCgGtT", searchPos + k);
		// Make sure we don't fall off the end
		if (invalidLoc < k + 1) {
		    return false;
		}
		searchPos = invalidLoc - k - 1;
	    }

	    // we found a hit, so make it the current k-mer
	    klen = k;
	    startPos = searchPos;
	    nextBaseIndex = searchPos + k + 1;
	    mer = tempMer;
	    return true;
	}

	bool isOutsideQuery(uint32_t jumpLen) {
	    uint32_t jumpPos = startPos + jumpLen;
	    return (jumpPos > qlen - k);
	}


	// move to the next valid k-mer that is at least skipVal past
	// the current position.  If skipVal positions forward is past
	// the end of the query, try and move to the last k-mer.
	//
	// If a valid k-mer is found, make it the current k-mer and
	// return true and the k-mer position. Otherwise, return
	// false.
	std::tuple<bool, uint32_t> skipForward(uint32_t skipVal) {
	   tempMer = mer;
	   uint32_t tempKLen = klen;
      	   uint32_t tempStartPos = startPos;
	   uint32_t tempNextBaseIndex = nextBaseIndex;

	   uint32_t jumpPos = startPos + skipVal;
	   // Would we jump past the end of the read?
	   // If so, just jump to the end.
	   if (jumpPos > qlen - k) {
	     jumpPos = qlen - k;
	   }
	   uint32_t initJumpPos = jumpPos;
	   if (jumpPos == startPos) { return std::make_pair(false, 0); }

	   // if the jumpPos is < k away, it's more efficient to just
	   // walk there b/c it overlaps the current k-mer
	   bool reachedGoal{false};
	   if (jumpPos - startPos < k) {
	       while (!reachedGoal and next()) {
		   reachedGoal = (startPos >= jumpPos);
	       }
	   } else {
	      // otherwise start a new k-mer at the jump position
	      while (! (reachedGoal = mer.from_chars(&qCharArray[jumpPos])) ) {
		  // If it wasn't a valid k-mer, find the offending base
		  // and try to start after it
		  jumpPos = qstr->find_first_not_of("aAcCgGtT", jumpPos) + 1;
		  // Make sure we don't fall off the end
		  if (jumpPos > qlen - k) {
		      startPos = invalidIndex;
		      reachedGoal = false;
		      break;
		  }
	      }
	      // If the search was successful
	      if (reachedGoal) {
		  // set startPos to the position we ended up jumping to
		  klen = k;
		  startPos = jumpPos;
		  nextBaseIndex = startPos + k + 1;
	      }
	   }
	   // If the search was un-successful, return the searcher to it's previous state
	   // and report the failure and the position where the backward search should begin.
	   if (!reachedGoal) {
	       mer = tempMer;
	       klen = tempKLen;
	       startPos = tempStartPos;
	       nextBaseIndex = tempNextBaseIndex;
	       return std::make_pair(false, initJumpPos);
	   }
	   return std::make_pair(true, startPos);
	}

	// Move to the next *valid* k-mer.  If we found a k-mer, return true,
	// If we can't move forward anymore, then return false.
	bool next() {
	    bool valid{false};
	    while (!valid and nextBaseIndex < qlen) {
		int c = jellyfish::mer_dna::code(qCharArray[nextBaseIndex]);
		// If the next base isn't a valid nucleotide
		if (jellyfish::mer_dna::not_dna(c)) {
		    // reset the k-mer
		    klen = 0;
		    ++nextBaseIndex;
		    if (qlen - nextBaseIndex < k) {
			startPos = invalidIndex;
			return false;
		    } else {
			continue;
		    }
		}
		mer.shift_left(c);
		//rcmer.shift_right(jellyfish::mer_dna::complement(c));
		++klen;
	        ++nextBaseIndex;
		if (klen >= k) {
		    startPos = nextBaseIndex - k;
		    valid = true;
		}
	    }
	    if (!valid) { startPos = invalidIndex; }
	    return valid;
	}

	bool isValid() { return startPos != invalidIndex; }
};

struct JumpStats {
    std::atomic<uint64_t> jumpSizes{0};
    std::atomic<uint64_t> numJumps{0};
};


class SkippingCollector {
    private:
	RapMapIndex* rmi_;
    public:
	SkippingCollector(RapMapIndex* rmiIn) : rmi_(rmiIn) {}

	void operator()(std::string& readStr,
		std::vector<QuasiAlignment>& hits,
		MateStatus mateStatus) {

	    auto jfhash = rmi_->merHash.get();
	    auto& kmerInfos = rmi_->kmerInfos;
	    auto& eqClasses = rmi_->eqClassList;
	    auto& eqClassLabels = rmi_->eqLabelList;
	    auto& posList = rmi_->posList;
	    auto& fwdJumpTable = rmi_->fwdJumpTable;
	    auto& revJumpTable = rmi_->revJumpTable;
	    auto posEnd = posList.end();

	    auto k = rapmap::utils::my_mer::k();
	    auto readLen = readStr.length();
	    uint32_t maxDist = static_cast<uint32_t>(readLen) * 1.5;

	    auto endIt = kmerInfos.end();

	    std::vector<HitInfo> kmerHits;
	    uint64_t merID;
	    size_t kID;
	    rapmap::utils::my_mer searchBuffer;

	    bool terminateSearch{false}; // terminate search after checking the *next* hit
	    bool validKmer{false};
	    uint32_t numAnchors{0};
	    uint32_t searchPos{0};
	    uint32_t jumpLen{0};
	    bool isRC;
	    SkippingKmerSearcher ksearch(readStr);

	    while (ksearch.isValid()) {
		auto searchMer = ksearch.getMer(isRC);
		bool foundMer = jfhash->get_val_for_key(searchMer, &merID,
			searchBuffer, &kID);
		// OK --- we found a hit
		if (foundMer) {

		    kmerHits.emplace_back(kmerInfos.begin() + merID,
			    merID, ksearch.queryIndex(), !isRC);
		    // Increment the # of anchor k-mers we found
		    ++numAnchors;

		    // If we decided to terminate the search in the last loop, then we're done.
		    if (terminateSearch) { break; }

		    jumpLen = isRC ? revJumpTable[merID] :
			      fwdJumpTable[merID];
		    //js.jumpSizes += jumpLen;
		    //++js.numJumps;

		    if (jumpLen > 1) { // only jump if it's worth it
			if (ksearch.isOutsideQuery(jumpLen)) {
			    if (ksearch.backwardSearch(searchPos)) {
				terminateSearch = true;
				continue;
			    } else {
				break;
			    }
			}

			std::tie(validKmer, searchPos) = ksearch.skipForward(jumpLen);
			if (!validKmer) {
			    // There was no valid k-mer from the skip position to the end
			    // of the read --- execute reverse search

			    // But -- only do so if we don't have at least 2 anchors
			    if (numAnchors >= 2) { break; }

			    if (ksearch.backwardSearch(searchPos)) {
				// If we find something in the reverse search, it will
				// be the last thing we check.
				terminateSearch = true;
				continue;
			    } else {
				// Otherwise, if we don't find anything in the reverse
				// search --- just give up and take what we have.
				break;
			    }
			} else {
			    if (searchPos == readLen - k) {
				terminateSearch = true;
			    }
			    // The skip was successful --- continue the search normally
			    // from here.
			    continue;
			}
		    }
		} else {
		    if (terminateSearch) { break; }
		}
		ksearch.next();
	    }

	    // found no hits in the entire read
	    if (kmerHits.size() == 0) { return; }

	    if (kmerHits.size() > 0) {
		if (kmerHits.size() > 1) {
		    //std::cerr << "kmerHits.size() = " << kmerHits.size() << "\n";
		    auto processedHits = rapmap::hit_manager::intersectHits(kmerHits, *rmi_);
		    rapmap::hit_manager::collectHitsSimple(processedHits, readLen, maxDist, hits, mateStatus);
		} else {
		    // std::cerr << "kmerHits.size() = " << kmerHits.size() << "\n";
		    auto& kinfo = *kmerHits[0].kinfo;
		    hits.reserve(kinfo.count);
		    // Iterator into, length of and end of the transcript list
		    auto& eqClassLeft = eqClasses[kinfo.eqId];
		    // Iterator into, length of and end of the positon list
		    auto leftPosIt = posList.begin() + kinfo.offset;
		    auto leftPosLen = kinfo.count;
		    auto leftPosEnd = leftPosIt + leftPosLen;
		    PositionListHelper leftPosHelper(leftPosIt, posList.end());
		    bool leftHitRC = kmerHits[0].queryRC;

		    auto leftTxpIt = eqClassLabels.begin() + eqClassLeft.txpListStart;
		    auto leftTxpListLen = eqClassLeft.txpListLen;
		    auto leftTxpEnd = leftTxpIt + leftTxpListLen;

		    for (auto it = leftTxpIt; it < leftTxpEnd; ++it) {
			collectAllHits(*it, readLen, leftHitRC, leftPosHelper, hits, mateStatus);
		    }
		}
	    }

	}
};



class EndCollector {
    private:
	RapMapIndex* rmi_;
    public:
	EndCollector(RapMapIndex* rmiIn) : rmi_(rmiIn) {}

    void operator()(std::string& readStr,
	    std::vector<QuasiAlignment>& hits,
	    MateStatus mateStatus) {

	auto jfhash = rmi_->merHash.get();
	auto& kmerInfos = rmi_->kmerInfos;
	auto& eqClasses = rmi_->eqClassList;
	auto& eqClassLabels = rmi_->eqLabelList;
	auto& posList = rmi_->posList;
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

	std::vector<HitInfo> kmerHits;
	bool leftFwd{true};
	bool rightFwd{true};

	uint64_t merID;
	size_t kID;
	rapmap::utils::my_mer searchBuffer;

	size_t klen{0};
	for (size_t i = 0; i < readLen; ++i) {
	    int c = jellyfish::mer_dna::code(readStr[i]);
	    // If the next base isn't a valid nucleotide
	    if (jellyfish::mer_dna::not_dna(c)) {
		// reset the k-mer
		klen = 0;
		continue;
	    }
	    mer.shift_left(c);
	    rcmer.shift_right(jellyfish::mer_dna::complement(c));
	    ++klen;
	    if (klen >= k) {
		auto& searchMer = (mer < rcmer) ? mer : rcmer;
		bool foundMer = jfhash->get_val_for_key(searchMer, &merID,
			searchBuffer, &kID);
		if (foundMer) {
		    kmerHits.emplace_back(kmerInfos.begin() + merID,
			    merID,
			    i - k,
			    searchMer == rcmer);
		    leftQueryPos = i - k;
		    break;
		}
	    }
	}

	// found no hits in the entire read
	if (kmerHits.size() == 0) { return; }

	// Now, start from the right and move left
	klen = 0;
	for (size_t i = readLen - 1; i > leftQueryPos; --i) {
	    int c = jellyfish::mer_dna::code(readStr[i]);
	    // If the next base isn't a valid nucleotide
	    if (jellyfish::mer_dna::not_dna(c)) {
		klen = 0;
		continue;
	    }
	    mer.shift_right(c);
	    rcmer.shift_left(jellyfish::mer_dna::complement(c));
	    ++klen;
	    if (klen >= k) {
		auto& searchMer = (mer < rcmer) ? mer : rcmer;
		bool foundMer = jfhash->get_val_for_key(searchMer, &merID,
			searchBuffer, &kID);
		if (foundMer) {
		    kmerHits.emplace_back(kmerInfos.begin() + merID,
			    merID,
			    readLen - (i + k),
			    searchMer == rcmer);
		    break;
		}
	    }
	}

	if (kmerHits.size() > 0) {
	    if (kmerHits.size() > 1) {
		//std::cerr << "kmerHits.size() = " << kmerHits.size() << "\n";
		auto processedHits = rapmap::hit_manager::intersectHits(kmerHits, *rmi_);
		rapmap::hit_manager::collectHitsSimple(processedHits, readLen, maxDist, hits, mateStatus);
	    } else {
		//std::cerr << "kmerHits.size() = " << kmerHits.size() << "\n";
		auto& kinfo = *kmerHits[0].kinfo;
		hits.reserve(kinfo.count);
		// Iterator into, length of and end of the transcript list
		auto& eqClassLeft = eqClasses[kinfo.eqId];
		// Iterator into, length of and end of the positon list
		auto leftPosIt = posList.begin() + kinfo.offset;
		auto leftPosLen = kinfo.count;
		auto leftPosEnd = leftPosIt + leftPosLen;
		PositionListHelper leftPosHelper(leftPosIt, posList.end());
		leftHitRC = kmerHits[0].queryRC;

		auto leftTxpIt = eqClassLabels.begin() + eqClassLeft.txpListStart;
		auto leftTxpListLen = eqClassLeft.txpListLen;
		auto leftTxpEnd = leftTxpIt + leftTxpListLen;

		for (auto it = leftTxpIt; it < leftTxpEnd; ++it) {
		    collectAllHits(*it, readLen, leftHitRC, leftPosHelper, hits, mateStatus);
		}
	    }
	}

    }
};

template <typename CollectorT, typename MutexT>
void processReadsSingle(single_parser* parser,
        RapMapIndex& rmi,
	CollectorT& hitCollector,
        MutexT* iomutex,
        std::ostream& outStream,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput) {

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

    char cigarBuff[1000];
    std::string readTemp(1000, 'A');
    std::string qualTemp(1000, '~');

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
            hitCollector(j->data[i].seq, hits, MateStatus::SINGLE_END);
            /*
               std::set_intersection(leftHits.begin(), leftHits.end(),
               rightHits.begin(), rightHits.end(),
               std::back_inserter(jointHits));
               */
            auto numHits = hits.size();
            hctr.totHits += numHits;

            if (hits.size() > 0 and hits.size() < maxNumHits) {
                auto& readName = j->data[i].header;
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                auto before = readName.find_first_of(':');
                before = readName.find_first_of(':', before+1);
                auto after = readName.find_first_of(':', before+1);
                const auto& txpName = readName.substr(before+1, after-before-1);
#endif //__DEBUG__
                uint32_t alnCtr{0};
                for (auto& qa : hits) {
                    auto& transcriptName = txpNames[qa.tid];
                    // === SAM
                    rapmap::utils::getSamFlags(qa, flags);
                    if (alnCtr != 0) {
                        flags |= 0x900;
                    }

                    std::string* readSeq = &(j->data[i].seq);
                    std::string* qstr = &(j->data[i].qual);
                    rapmap::utils::adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], cigarStr);

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
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    if (txpNames[qa.tid] == txpName) { ++hctr.trueHits; }
#endif //__DEBUG__
                }
            }

            if (hctr.numReads > hctr.lastPrint + 1000000) {
		hctr.lastPrint.store(hctr.numReads.load());
                if (iomutex->try_lock()){
                    if (hctr.numReads > 0) {
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                        std::cerr << "\033[F\033[F\033[F";
#else
                        std::cerr << "\033[F\033[F";
#endif // __DEBUG__
                    }
                    std::cerr << "saw " << hctr.numReads << " reads\n";
                    std::cerr << "# hits per read = "
                        << hctr.totHits / static_cast<float>(hctr.numReads) << "\n";
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    std::cerr << "The true hit was in the returned set of hits "
                        << 100.0 * (hctr.trueHits / static_cast<float>(hctr.numReads))
                        <<  "% of the time\n";
#endif // __DEBUG__
                    iomutex->unlock();
                }
            }
        } // for all reads in this job

        // DUMP OUTPUT
        iomutex->lock();
        outStream << sstream.str();
        iomutex->unlock();
        sstream.clear();

    } // processed all reads
}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename CollectorT, typename MutexT>
void processReadsPair(paired_parser* parser,
        RapMapIndex& rmi,
	CollectorT& hitCollector,
        MutexT* iomutex,
        std::ostream& outStream,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput) {
    auto& txpNames = rmi.txpNames;
    std::vector<uint32_t>& txpLens = rmi.txpLens;
    uint32_t n{0};
    uint32_t k = rapmap::utils::my_mer::k();
    std::vector<std::string> transcriptNames;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    auto logger = spdlog::get("stderrLog");

    fmt::MemoryWriter sstream;
    size_t batchSize{1000};
    std::vector<QuasiAlignment> leftHits;
    std::vector<QuasiAlignment> rightHits;
    std::vector<QuasiAlignment> jointHits;

    size_t readLen{0};
	bool tooManyHits{false};
    uint16_t flags1, flags2;
    // 1000-bp reads are max here (get rid of hard limit later).
    std::string read1Temp(1000, 'A');
    std::string qual1Temp(1000, '~');
    std::string read2Temp(1000, 'A');
    std::string qual2Temp(1000, '~');

    char buff1[1000];
    char buff2[1000];
    FixedWriter cigarStr1(buff1, 1000);
    FixedWriter cigarStr2(buff2, 1000);

    JumpStats js;
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
		    tooManyHits = false;
            readLen = j->data[i].first.seq.length();
            ++hctr.numReads;
            jointHits.clear();
            leftHits.clear();
            rightHits.clear();
	    hitCollector(j->data[i].first.seq,
                        leftHits, MateStatus::PAIRED_END_LEFT);
            hitCollector(j->data[i].second.seq,
                        rightHits, MateStatus::PAIRED_END_RIGHT);
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
					uint32_t leftTxp, rightTxp;
                    while (leftIt != leftEnd && rightIt != rightEnd) {
                        leftTxp = leftIt->tid;
                        rightTxp = rightIt->tid;
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
                                jointHits.back().mateStatus = MateStatus::PAIRED_END_PAIRED;

                                ++numHits;
                                if (numHits > maxNumHits) { tooManyHits = true; break; }
                                ++leftIt;
                            }
                            ++rightIt;
                        }
                    }
                }
				if (tooManyHits) { jointHits.clear(); ++hctr.tooManyHits; }
            }

			// If we had proper paired hits
            if (jointHits.size() > 0) {
                hctr.peHits += jointHits.size();
                orphanStatus = 0;
            } else if (leftHits.size() + rightHits.size() > 0 and !tooManyHits) {
		    // If there weren't proper paired hits, then either
			// there were too many hits, and we forcibly discarded the read
			// or we take the single end hits.
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

            if (jointHits.size() > 0 and !noOutput) {
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


#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                auto before = readName.find_first_of(':');
                before = readName.find_first_of(':', before+1);
                auto after = readName.find_first_of(':', before+1);
                const auto& trueTxpName = readName.substr(before+1, after-before-1);
#endif //__DEBUG__
                uint32_t alnCtr{0};
				uint32_t trueHitCtr{0};
				QuasiAlignment* firstTrueHit{nullptr};
                for (auto& qa : jointHits) {
                    auto& transcriptName = txpNames[qa.tid];
                    // === SAM
                    if (qa.isPaired) {
                        rapmap::utils::getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            flags2 |= 0x900;
                        }
                        rapmap::utils::adjustOverhang(qa, txpLens[qa.tid], cigarStr1, cigarStr2);

                        // Reverse complement the read and reverse
                        // the quality string if we need to
                        if (!qa.fwd) {
                            rapmap::utils::reverseRead(j->data[i].first.seq,
                                        j->data[i].first.qual,
                                        read1Temp,
                                        qual1Temp);
                        }
                        if (!qa.mateIsFwd) {
                            rapmap::utils::reverseRead(j->data[i].second.seq,
                                        j->data[i].second.qual,
                                        read2Temp,
                                        qual2Temp);
                        }

                        sstream << readName.c_str() << '\t' // QNAME
                                << flags1 << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.pos + 1 << '\t' // POS (1-based)
                                << 1 << '\t' // MAPQ
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
                                << 1 << '\t' // MAPQ
                                << cigarStr2.c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << qa.pos + 1 << '\t' // PNEXT
                                << qa.fragLen << '\t' // TLEN
                                << j->data[i].second.seq << '\t' // SEQ
                                << j->data[i].first.qual << '\n';
                    } else {
                        rapmap::utils::getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            // If this is the first alignment for this read
                            // If the left end is mapped, set 0x900 on the right end
                            if (qa.mateStatus == MateStatus::PAIRED_END_LEFT) {
                                flags2 |= 0x900;
                            } else {
                            // Otherwise, set 0x900 on the left end
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
                        if (qa.mateStatus == MateStatus::PAIRED_END_LEFT) { // left read
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

                        // Reverse complement the read and reverse
                        // the quality string if we need to
                        if (!qa.fwd) {
                            rapmap::utils::reverseRead(*readSeq, *qstr,
                                        read1Temp, qual1Temp);
                        }

                        rapmap::utils::adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], *cigarStr);
                        sstream << readName.c_str() << '\t' // QNAME
                                << flags << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.pos + 1 << '\t' // POS (1-based)
                                << 1 << '\t' // MAPQ
                                << cigarStr->c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << 0 << '\t' // PNEXT (only 1 read in templte)
                                << 0 << '\t' // TLEN (spec says 0, not read len)
                                << *readSeq << '\t' // SEQ
                                << *qstr << '\n';

                        // Output the info for the unaligned mate.
                        sstream << unalignedName->c_str() << '\t' // QNAME
                            << unalignedFlags << '\t' // FLAGS
                            << transcriptName << '\t' // RNAME (same as mate)
                            << qa.pos + 1 << '\t' // POS (same as mate)
                            << 0 << '\t' // MAPQ
                            << readSeq->length() << 'M' << '\t' // CIGAR
                            << '=' << '\t' // RNEXT
                            << 0 << '\t' // PNEXT (only 1 read in template)
                            << 0 << '\t' // TLEN (spec says 0, not read len)
                            << *unalignedSeq << '\t' // SEQ
                            << *unalignedQstr << '\n';
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
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    if (transcriptName == trueTxpName) {
							if (trueHitCtr == 0) {
									++hctr.trueHits;
									++trueHitCtr;
									firstTrueHit = &qa;
							} else {
									++trueHitCtr;
									std::cerr << "Found true hit " << trueHitCtr << " times!\n";
									std::cerr << transcriptName << '\t' << firstTrueHit->pos
											<< '\t' << firstTrueHit->fwd << '\t' << firstTrueHit->fragLen
											<< '\t' << (firstTrueHit->isPaired ? "Paired" : "Orphan") << '\t';
								    printMateStatus(firstTrueHit->mateStatus);
								    std::cerr << '\n';
									std::cerr << transcriptName << '\t' << qa.pos
											  << '\t' << qa.fwd << '\t' << qa.fragLen
										      << '\t' << (qa.isPaired ? "Paired" : "Orphan") << '\t';
								    printMateStatus(qa.mateStatus);
								    std::cerr << '\n';

							}
					}
#endif //__DEBUG__
                }
            }

            if (hctr.numReads > hctr.lastPrint + 1000000) {
		hctr.lastPrint.store(hctr.numReads.load());
                if (iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
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
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    std::cerr << "The true hit was in the returned set of hits "
                        << 100.0 * (hctr.trueHits / static_cast<float>(hctr.numReads))
                        <<  "% of the time\n";
		    /*
                    std::cerr << "Average jump size = "
                              << js.jumpSizes / static_cast<double>(js.numJumps) << "\n";
			      */
#endif // __DEBUG__
                    iomutex->unlock();
                }
            }
        } // for all reads in this job

        // DUMP OUTPUT
        if (!noOutput) {
            iomutex->lock();
            outStream << sstream.str();
            iomutex->unlock();
            sstream.clear();
        }

    } // processed all reads

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
    TCLAP::ValueArg<std::string> read1("1", "leftMates", "The location of the left paired-end reads", false, "", "path");
    TCLAP::ValueArg<std::string> read2("2", "rightMates", "The location of the right paired-end reads", false, "", "path");
    TCLAP::ValueArg<std::string> unmatedReads("r", "unmatedReads", "The location of single-end reads", false, "", "path");
    TCLAP::ValueArg<uint32_t> numThreads("t", "numThreads", "Number of threads to use", false, 1, "positive integer");
    TCLAP::ValueArg<uint32_t> maxNumHits("m", "maxNumHits", "Reads mapping to more than this many loci are discarded", false, 200, "positive integer");
    TCLAP::ValueArg<std::string> outname("o", "output", "The output file (default: stdout)", false, "", "path");
    TCLAP::SwitchArg endCollectorSwitch("e", "endCollector", "Use the simpler (and faster) \"end\" collector as opposed to the more sophisticated \"skipping\" collector", false);
    TCLAP::SwitchArg noout("n", "noOutput", "Don't write out any alignments (for speed testing purposes)", false);
    cmd.add(index);
    cmd.add(noout);

    cmd.add(read1);
    cmd.add(read2);
    cmd.add(unmatedReads);
    cmd.add(outname);
    cmd.add(numThreads);
    cmd.add(maxNumHits);
    cmd.add(endCollectorSwitch);

    auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
    auto consoleLog = spdlog::create("stderrLog", {consoleSink});

    try {

	cmd.parse(argc, argv);
	bool pairedEnd = (read1.isSet() or read2.isSet());
	if (pairedEnd and (read1.isSet() != read2.isSet())) {
	    consoleLog->error("You must set both the -1 and -2 arguments to align "
		    "paired end reads!");
	    std::exit(1);
	}

	if (pairedEnd and unmatedReads.isSet()) {
	    consoleLog->error("You cannot specify both paired-end and unmated "
		    "reads in the input!");
	    std::exit(1);
	}

	if (!pairedEnd and !unmatedReads.isSet()) {
	    consoleLog->error("You must specify input; either both paired-end "
			      "or unmated reads!");
	    std::exit(1);

	}

	std::string indexPrefix(index.getValue());
	if (indexPrefix.back() != '/') {
	    indexPrefix += "/";
	}

	if (!rapmap::fs::DirExists(indexPrefix.c_str())) {
	    consoleLog->error("It looks like the index you provided [{}] "
		    "doesn't exist", indexPrefix);
	    std::exit(1);
	}

	RapMapIndex rmi;
	rmi.load(indexPrefix);

	std::cerr << "\n\n\n\n";

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

	if (!noout.getValue()) {
        rapmap::utils::writeSAMHeader(rmi, outStream);
	}

	SpinLockT iomutex;
	{
	    ScopedTimer timer;
	    HitCounters hctrs;
	    consoleLog->info("mapping reads . . . \n\n\n");
	    if (pairedEnd) {
		std::vector<std::thread> threads;
		std::vector<std::string> read1Vec = rapmap::utils::tokenize(read1.getValue(), ',');
		std::vector<std::string> read2Vec = rapmap::utils::tokenize(read2.getValue(), ',');

		if (read1Vec.size() != read2Vec.size()) {
		    consoleLog->error("The number of provided files for "
			    "-1 and -2 must be the same!");
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

		/** Create the threads depending on the collector type **/
		if (endCollectorSwitch.getValue()) {
		    EndCollector endCollector(&rmi);
		    for (size_t i = 0; i < nthread; ++i) {
			threads.emplace_back(processReadsPair<EndCollector, SpinLockT>,
				pairParserPtr.get(),
				std::ref(rmi),
				std::ref(endCollector),
				&iomutex,
				std::ref(outStream),
				std::ref(hctrs),
				maxNumHits.getValue(),
				noout.getValue());
		    }
		} else {
		    SkippingCollector skippingCollector(&rmi);
		    for (size_t i = 0; i < nthread; ++i) {
			threads.emplace_back(processReadsPair<SkippingCollector, SpinLockT>,
				pairParserPtr.get(),
				std::ref(rmi),
				std::ref(skippingCollector),
				&iomutex,
				std::ref(outStream),
				std::ref(hctrs),
				maxNumHits.getValue(),
				noout.getValue());
		    }
		}

		for (auto& t : threads) { t.join(); }
		delete [] pairFileList;
	    } else {
		std::vector<std::thread> threads;
		std::vector<std::string> unmatedReadVec = rapmap::utils::tokenize(unmatedReads.getValue(), ',');
		size_t maxReadGroup{1000}; // Number of reads in each "job"
		size_t concurrentFile{1};
		stream_manager streams( unmatedReadVec.begin(), unmatedReadVec.end(),
			concurrentFile);
		singleParserPtr.reset(new single_parser(4 * nthread,
			    maxReadGroup,
			    concurrentFile,
			    streams));

		/** Create the threads depending on the collector type **/
		if (endCollectorSwitch.getValue()) {
		    EndCollector endCollector(&rmi);
		    for (size_t i = 0; i < nthread; ++i) {
			threads.emplace_back(processReadsSingle<EndCollector, SpinLockT>,
				singleParserPtr.get(),
				std::ref(rmi),
				std::ref(endCollector),
				&iomutex,
				std::ref(outStream),
				std::ref(hctrs),
				maxNumHits.getValue(),
				noout.getValue());
		    }
		} else {
		    SkippingCollector skippingCollector(&rmi);
		    for (size_t i = 0; i < nthread; ++i) {
			threads.emplace_back(processReadsSingle<SkippingCollector, SpinLockT>,
				singleParserPtr.get(),
				std::ref(rmi),
				std::ref(skippingCollector),
				&iomutex,
				std::ref(outStream),
				std::ref(hctrs),
				maxNumHits.getValue(),
				noout.getValue());
		    }
		}
		for (auto& t : threads) { t.join(); }
	    }
	    consoleLog->info("done mapping reads.");
	    consoleLog->info("Discarded {} reads because they had > {} alignments",
		    hctrs.tooManyHits, maxNumHits.getValue());

	}

	if (haveOutputFile) {
	    outFile.close();
	}
	return 0;
    } catch (TCLAP::ArgException& e) {
	consoleLog->error("Exception [{}] when parsing argument {}", e.error(), e.argId());
	return 1;
    }

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


void collectHits(RapMapIndex& rmi, std::string& readStr,
                 std::vector<QuasiAlignment>& hits,
                 MateStatus mateStatus) {

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
	uint32_t klen{0};
    rapmap::utils::my_mer searchBuffer;

    for (size_t i = 0; i < readLen; ++i) {
        int c = jellyfish::mer_dna::code(readStr[i]);
        if (jellyfish::mer_dna::not_dna(c)) {
		    klen = 0;
			continue;
        }
        mer.shift_left(c);
        rcmer.shift_right(jellyfish::mer_dna::complement(c));
		++klen;
        if (klen >= k) {
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

    // Now, start from the right and move left
	klen = 0;
    for (size_t i = readLen - 1; i > leftQueryPos; --i) {
        int c = jellyfish::mer_dna::code(readStr[i]);
        if (jellyfish::mer_dna::not_dna(c)) {
		  //c = jellyfish::mer_dna::code('G');
		  klen = 0;
		  continue;
        }
        mer.shift_right(c);
        rcmer.shift_left(jellyfish::mer_dna::complement(c));
		++klen;
        if (klen >= k) {
            auto& searchMer = (mer < rcmer) ? mer : rcmer;
            bool foundMer = jfhash->get_val_for_key(searchMer, &merID,
                                                    searchBuffer, &kID);
            if (foundMer) {
                miniRightHits = kmerInfos.begin() + merID;
				//if (miniLeftHits == miniRightHits) { continue; }
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

*/
