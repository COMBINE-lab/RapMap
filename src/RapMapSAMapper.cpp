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
#include <fstream>
#include <iostream>
#include <tuple>
#include <memory>
#include <cstring>

#include "ScopedTimer.hpp"

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

#include "HitManager.hpp"
//#include "SIMDCompressionAndIntersection/intersection.h"
#include "xxhash.h"

#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"
#include "jellyfish/hash_counter.hpp"

#include "tclap/CmdLine.h"

#include <RapMapAligner.hpp>

#include "stringpiece.h"

#include "PairSequenceParser.hpp"
#include "PairAlignmentFormatter.hpp"
#include "SingleAlignmentFormatter.hpp"
#include "RapMapUtils.hpp"
#include "RapMapSAIndex.hpp"
#include "RapMapFileSystem.hpp"
#include "RapMapConfig.hpp"
#include "ScopedTimer.hpp"
#include "SpinLock.hpp"
#include "IndexHeader.hpp"
#include "SASearcher.hpp"
#include "SACollector.hpp"

//#define __TRACK_CORRECT__
//#define CIGAR_VERIFY 1

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

using HitCounters = rapmap::utils::HitCounters;
using MateStatus = rapmap::utils::MateStatus;
using HitInfo = rapmap::utils::HitInfo;
using ProcessedHit = rapmap::utils::ProcessedHit;
using QuasiAlignment = rapmap::utils::QuasiAlignment;
using FixedWriter = rapmap::utils::FixedWriter;

// This assumes the readLength < transcriptLength
template <typename RapMapIndexT>
inline int alignHit(RapMapIndexT& rmi,
                    RapMapAligner& beforeAligner,
                    RapMapAligner& afterAligner,
                    std::string& read,
                    uint32_t tid,
                    int32_t& pos,
                    uint32_t queryPos,
                    uint32_t matchLen,
                    std::string& cigar,
                    bool align) {
    auto readLen = read.length();
    CigarString cigarString;
    auto extraAlignFactor = 2;
    int score = 0;

    //     txpHitStart
    //       |txpStart                   txpHitEnd             txpEnd
    //       |   |txpMatchStart  txpMatchEnd |                    |
    //       v   v    v               v      v                    v
    // ----------|----================----------------------------|---------|
    //       |---|----================|------|
    //                ^               ^
    //            matchStart       matchEnd

    // Get copy of transcript
    int64_t matchStart = queryPos;
    int64_t matchEnd = matchStart + matchLen;
    int64_t txpLen = rmi.txpLens[tid];
    int64_t txpStart = rmi.txpOffsets[tid];
    int64_t txpEnd = txpStart + txpLen;
    int64_t txpHitStart = txpStart + pos;
    int64_t txpHitEnd = txpHitStart + readLen;
    int64_t txpMatchStart = txpHitStart + matchStart;
    int64_t txpMatchEnd = txpHitStart + matchEnd;

    // If the read is over hanging before/after the transcript
    int64_t clipBefore = std::max(0, -pos);
    int64_t clipAfter = std::max(0L, txpHitEnd - txpEnd);
    int64_t alignBeforeLen = std::max(0L, matchStart - clipBefore);
    int64_t txpAlignBeforeStart = txpMatchStart - alignBeforeLen;
    int64_t readAlignBeforeStart = matchStart - alignBeforeLen;
    int64_t alignAfterLen = std::max(0L, txpHitEnd - clipAfter - txpMatchEnd);
    int64_t txpAlignAfterStart = txpMatchEnd;
    int64_t readAlignAfterStart = matchEnd;

    // Assumptions
    assert(txpMatchEnd <= txpEnd);
    assert(txpMatchStart >= txpStart);

    if (clipAfter > 0) {
        cigarString.emplace_front(CigarOp::S, clipAfter);
    }

    if (alignAfterLen > 0) {
        // align after with some extra space for aligning to transcript
        int64_t extraLen = std::min(extraAlignFactor * alignAfterLen, txpEnd - txpAlignAfterStart);
        // Align with *free end gaps* in the transcript
        if (align) {
            score += afterAligner.align(rmi.seq, txpAlignAfterStart, extraLen + alignAfterLen,
                                        read, readAlignAfterStart, alignAfterLen,
                                        cigarString);
        } else {
            cigarString.emplace_front(CigarOp::M, alignAfterLen);
            score += alignAfterLen;
        }
    }

    // align match
    score += matchLen * beforeAligner.match;
    cigarString.emplace_front(CigarOp::EQ, matchLen);

    if (alignBeforeLen > 0) {
        // align before with some extra space for aligning to transcript
        int64_t extraLen = std::min(extraAlignFactor * alignBeforeLen, txpAlignBeforeStart - txpStart);
        int64_t txpExtraStart = txpAlignBeforeStart - extraLen;
        // Align with *free begin gaps* in the transcript
        if (align) {
            score += beforeAligner.align(rmi.seq, txpExtraStart, extraLen + alignBeforeLen,
                                         read, readAlignBeforeStart, alignBeforeLen,
                                         cigarString);
        } else {
            cigarString.emplace_front(CigarOp::M, alignBeforeLen);
            score += alignBeforeLen;
        }

    }

    // Clip *before* the match
    if (clipBefore > 0) {
        cigarString.emplace_front(CigarOp::S, clipBefore);
    }
    cigarString.toString(cigar);

    //Verify CIGAR
#ifdef CIGAR_VERIFY

    std::string tempTxp = rmi.seq.substr(txpMatchStart); //do something here
    int txpIndex = 0;
    std::string tempRead = read;
    int readIndex = 0;
    std::string reconstructedTxp;

    //clip before
    CigarElement f = cigarString.cigar.front();
    if(f.op==CigarOp::S){
    	//std::cout << "F: "<< (char)(f.op) << ":" << f.count << std::endl;
    	reconstructedTxp+=tempRead.substr(readIndex,f.count);
    	readIndex+=f.count;
    }

    for(CigarElement c: cigarString.cigar){
    	int len = c.count;
    	if(c.op==CigarOp::EQ){ //sequence match
    		reconstructedTxp+=tempRead.substr(readIndex,len);
    		readIndex+=len;
    		txpIndex+=len;
    	} else if(c.op==CigarOp::X){ //sequence mismatch
    		reconstructedTxp+=tempRead.substr(readIndex,len);
    		readIndex+=len;
    		txpIndex+=len;
    	} else if(c.op==CigarOp::I){ // insertion into the reference
    		reconstructedTxp+=tempRead.substr(readIndex,len);
    		readIndex+=len;
    	} else if(c.op==CigarOp::D){ //deletion from the reference
    		txpIndex+=len;
    	} else {
    		continue;
    	}
    }

    //clip after
    CigarElement b = cigarString.cigar.back();
    if(b.op==CigarOp::S){
        reconstructedTxp+=tempRead.substr(readIndex,b.count);
        readIndex+=b.count;
    }

    tempRead = read.substr(matchStart,matchLen);
    if(reconstructedTxp!=read){
    	std::cout << reconstructedTxp << "\n" << read << std::endl;
    	std::cout << rmi.seq.substr(txpMatchStart,readLen) << std::endl;
    	std::cout << matchStart << "\t" << matchEnd << "\t" << txpLen << "\t" << txpStart << "\t" << txpEnd << "\t" << txpHitStart << "\t" << txpHitEnd << "\t" << txpMatchStart << "\t" << txpMatchEnd << std::endl;
    } else {
    	std::cout << ".";
    }

    assert(reconstructedTxp==read);
#endif


    return score;
}

template <typename RapMapIndexT>
void alignSingleAll(RapMapIndexT& rmi,
                    RapMapAligner& beforeAligner,
                    RapMapAligner& afterAligner,
                    std::string& fwdRead,
                    std::vector<rapmap::utils::QuasiAlignment>& hits,
                    bool align) {
    int64_t readLen = fwdRead.length();
    std::string revRead;
    for (int64_t i = readLen - 1; i >= 0; --i) {
        revRead += static_cast<char>(rapmap::utils::rc_table[(int8_t) fwdRead[i]]);
    }
    for (auto& qa : hits) {
        // Forward or reverse?
        std::string& read = qa.fwd ? fwdRead : revRead;
        qa.score = alignHit(rmi, beforeAligner, afterAligner, read, qa.tid,
                            qa.pos, qa.queryPos, qa.matchLen, qa.cigar, align);
    }
}

template <typename RapMapIndexT>
void alignPairedAll(RapMapIndexT& rmi,
                    RapMapAligner& beforeAligner,
                    RapMapAligner& afterAligner,
                    std::string& fwdLeftRead,
                    std::string& fwdRightRead,
                    std::vector<rapmap::utils::QuasiAlignment>& hits,
                    bool align) {
    std::string revLeftRead, revRightRead;
    for (int64_t i = fwdLeftRead.length() - 1; i >= 0; --i) {
        revLeftRead += static_cast<char>(rapmap::utils::rc_table[(int8_t) fwdLeftRead[i]]);
    }
    for (int64_t i = fwdRightRead.length() - 1; i >= 0; --i) {
        revRightRead += static_cast<char>(rapmap::utils::rc_table[(int8_t) fwdRightRead[i]]);
    }
    for (auto& qa : hits) {
        // Forward or reverse?
        std::string& leftRead = qa.fwd ? fwdLeftRead : revLeftRead;
        std::string& rightRead = qa.mateIsFwd ? fwdRightRead : revRightRead;
        switch (qa.mateStatus) {
            case MateStatus::PAIRED_END_PAIRED:
                qa.score = alignHit(rmi, beforeAligner, afterAligner,
                                    leftRead, qa.tid, qa.pos, qa.queryPos,
                                    qa.matchLen, qa.cigar, align);
                qa.mateScore = alignHit(rmi, beforeAligner, afterAligner,
                                        rightRead, qa.tid, qa.matePos,
                                        qa.mateQueryPos, qa.mateMatchLen,
                                        qa.mateCigar, align);
                break;
            case MateStatus::PAIRED_END_LEFT:
                qa.score = alignHit(rmi, beforeAligner, afterAligner,
                                    leftRead, qa.tid, qa.pos, qa.queryPos,
                                    qa.matchLen, qa.cigar, align);
                break;
            case MateStatus::PAIRED_END_RIGHT:
                qa.score = alignHit(rmi, beforeAligner, afterAligner,
                                    rightRead, qa.tid, qa.pos, qa.queryPos,
                                    qa.matchLen, qa.cigar, align);
                break;
            default:
                break;
        }
    }
}


template <typename RapMapIndexT, typename CollectorT, typename MutexT>
void processReadsSingleSA(single_parser * parser,
        RapMapIndexT& rmi,
    	CollectorT& hitCollector,
        MutexT* iomutex,
    	std::shared_ptr<spdlog::logger> outQueue,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput,
        bool strictCheck,
        bool doAlign) {

    using OffsetT = typename RapMapIndexT::IndexType;
    auto& txpNames = rmi.txpNames;
    auto& txpLens = rmi.txpLens;
    uint32_t n{0};
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    auto logger = spdlog::get("stderrLog");

    fmt::MemoryWriter sstream;
    size_t batchSize{2500};
    std::vector<QuasiAlignment> hits;

    size_t readLen{0};
	bool tooManyHits{false};
    uint16_t flags;

    SingleAlignmentFormatter<RapMapIndexT*> formatter(&rmi);

    SASearcher<RapMapIndexT> saSearcher(&rmi);

    uint32_t orphanStatus{0};

    RapMapAligner beforeAligner(true, false);
    RapMapAligner afterAligner(false, true);

    while(true) {
        typename single_parser::job j(*parser); // Get a job from the parser: a bunch of reads (at most max_read_group)
        if(j.is_empty()) break;                 // If we got nothing, then quit.
        for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
            //readLen = j->data[i].seq.length();
            ++hctr.numReads;
            hits.clear();
            hitCollector(j->data[i].seq, hits, saSearcher, MateStatus::SINGLE_END, strictCheck);
            auto numHits = hits.size();
            hctr.totHits += numHits;

	    if (hits.size() > 0 and !noOutput and hits.size() <= maxNumHits) {
                /*
                std::sort(hits.begin(), hits.end(),
                            [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                return a.tid < b.tid;
                            });
                */
            alignSingleAll(rmi, beforeAligner, afterAligner, j->data[i].seq, hits, doAlign);
            if (!noOutput) {
                rapmap::utils::writeAlignmentsToStream(j->data[i], formatter,
                                                       hctr, hits, sstream);
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
                    std::cerr << "# matches per hit = "
                        << hctr.totMatchLens / static_cast<float>(hctr.totHits) << "\n";
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
        if (!noOutput) {
            std::string outStr(sstream.str());
            // Get rid of last newline
            if (!outStr.empty()) {
                outStr.pop_back();
                outQueue->info() << std::move(outStr);
            }
            sstream.clear();
            /*
             iomutex->lock();
             outStream << sstream.str();
             iomutex->unlock();
             sstream.clear();
             */
        }

    } // processed all reads


}

/**
 *  Map reads from a collection of paired-end files.
 */
template <typename RapMapIndexT, typename CollectorT, typename MutexT>
void processReadsPairSA(paired_parser* parser,
        RapMapIndexT& rmi,
    	CollectorT& hitCollector,
        MutexT* iomutex,
	    std::shared_ptr<spdlog::logger> outQueue,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput,
        bool strictCheck,
        bool nonStrictMerge,
        bool doAlign) {

    using OffsetT = typename RapMapIndexT::IndexType;

    auto& txpNames = rmi.txpNames;
    auto& txpLens = rmi.txpLens;
    uint32_t n{0};
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

    // Create a formatter for alignments
    PairAlignmentFormatter<RapMapIndexT*> formatter(&rmi);

    SASearcher<RapMapIndexT> saSearcher(&rmi);

    uint32_t orphanStatus{0};
    RapMapAligner beforeAligner(true, false);
    RapMapAligner afterAligner(false, true);

    while(true) {
        typename paired_parser::job j(*parser); // Get a job from the parser: a bunch of reads (at most max_read_group)
        if(j.is_empty()) break;                 // If we got nothing, quit
        for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
		    tooManyHits = false;
            readLen = j->data[i].first.seq.length();
            ++hctr.numReads;
            jointHits.clear();
            leftHits.clear();
            rightHits.clear();

            bool lh = hitCollector(j->data[i].first.seq,
                        leftHits, saSearcher,
                        MateStatus::PAIRED_END_LEFT,
                        strictCheck);
            bool rh = hitCollector(j->data[i].second.seq,
                        rightHits, saSearcher,
                        MateStatus::PAIRED_END_RIGHT,
                        strictCheck);

            if (nonStrictMerge) {
                rapmap::utils::mergeLeftRightHitsFuzzy(
                        lh, rh,
                        leftHits, rightHits, jointHits,
                        readLen, maxNumHits, tooManyHits, hctr);

            } else {
                rapmap::utils::mergeLeftRightHits(
                        leftHits, rightHits, jointHits,
                        readLen, maxNumHits, tooManyHits, hctr);
            }

            for (auto& h : jointHits) {
                hctr.totMatchLens += h.matchLen;
            }
            // If we have reads to output, and we're writing output.
            if (jointHits.size() > 0 and !noOutput and jointHits.size() <= maxNumHits) {
                alignPairedAll(rmi, beforeAligner, afterAligner,
                               j->data[i].first.seq, j->data[i].second.seq, jointHits, doAlign);
                if (!noOutput) {
                    rapmap::utils::writeAlignmentsToStream(j->data[i], formatter,
                                                           hctr, jointHits, sstream);
                }
            }

            if (hctr.numReads > hctr.lastPrint + 1000000) {
        		hctr.lastPrint.store(hctr.numReads.load());
                if (iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
                        std::cerr << "\r\r";
                    }
                    std::cerr << "saw " << hctr.numReads << " reads : "
                              << "pe / read = " << hctr.peHits / static_cast<float>(hctr.numReads)
                              << " : se / read = " << hctr.seHits / static_cast<float>(hctr.numReads) << " "
                              << " : matches / hit = "
                              << hctr.totMatchLens / static_cast<float>(hctr.peHits + hctr.seHits) << " ";
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    std::cerr << ": true hit \% = "
                        << (100.0 * (hctr.trueHits / static_cast<float>(hctr.numReads)));
#endif // __DEBUG__
                    iomutex->unlock();
                }
            }
        } // for all reads in this job

        // DUMP OUTPUT
        if (!noOutput) {
            std::string outStr(sstream.str());
            // Get rid of last newline
            if (!outStr.empty()) {
                outStr.pop_back();
                outQueue->info() << std::move(outStr);
            }
            sstream.clear();
	        /*
            iomutex->lock();
            outStream << sstream.str();
            iomutex->unlock();
            sstream.clear();
	        */
        }

    } // processed all reads

}

template <typename RapMapIndexT, typename MutexT>
bool spawnProcessReadsThreads(
        uint32_t nthread,
        paired_parser* parser,
        RapMapIndexT& rmi,
        MutexT& iomutex,
	      std::shared_ptr<spdlog::logger> outQueue,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput,
        bool strictCheck,
        bool fuzzy,
        bool align) {

            std::vector<std::thread> threads;
            SACollector<RapMapIndexT> saCollector(&rmi);
            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsPairSA<RapMapIndexT, SACollector<RapMapIndexT>, MutexT>,
                        parser,
                        std::ref(rmi),
                        std::ref(saCollector),
                        &iomutex,
            			      outQueue,
                        std::ref(hctr),
                        maxNumHits,
                        noOutput,
                        strictCheck,
                        fuzzy,
                        align);
            }

            for (auto& t : threads) { t.join(); }
            return true;
        }

template <typename RapMapIndexT, typename MutexT>
bool spawnProcessReadsThreads(
        uint32_t nthread,
        single_parser* parser,
        RapMapIndexT& rmi,
        MutexT& iomutex,
	      std::shared_ptr<spdlog::logger> outQueue,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput,
        bool strictCheck,
        bool align) {

            std::vector<std::thread> threads;
            SACollector<RapMapIndexT> saCollector(&rmi);
            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsSingleSA<RapMapIndexT, SACollector<RapMapIndexT>, MutexT>,
                        parser,
                        std::ref(rmi),
                        std::ref(saCollector),
                        &iomutex,
            			      outQueue,
                        std::ref(hctr),
                        maxNumHits,
                        noOutput,
                        strictCheck,
                        align);
            }
            for (auto& t : threads) { t.join(); }
            return true;
        }

int rapMapSAMap(int argc, char* argv[]) {
    std::cerr << "RapMap Mapper (SA-based)\n";

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
    TCLAP::SwitchArg noout("n", "noOutput", "Don't write out any alignments (for speed testing purposes)", false);
    TCLAP::SwitchArg strict("s", "strictCheck", "Perform extra checks to try and assure that only equally \"best\" mappings for a read are reported", false);
    TCLAP::SwitchArg fuzzy("f", "fuzzyIntersection", "Find paired-end mapping locations using fuzzy intersection", false);
    TCLAP::SwitchArg align("a", "align", "Compute the optimal alignments (CIGAR) for each mapping position", false);
    cmd.add(index);
    cmd.add(noout);

    cmd.add(read1);
    cmd.add(read2);
    cmd.add(unmatedReads);
    cmd.add(outname);
    cmd.add(numThreads);
    cmd.add(maxNumHits);
    cmd.add(strict);
    cmd.add(fuzzy);
    cmd.add(align);

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

	IndexHeader h;
	std::ifstream indexStream(indexPrefix + "header.json");
	{
		cereal::JSONInputArchive ar(indexStream);
		ar(h);
	}
	indexStream.close();

	if (h.indexType() != IndexType::QUASI) {
	    consoleLog->error("The index {} does not appear to be of the "
			    "appropriate type (quasi)", indexPrefix);
	    std::exit(1);
	}

  std::unique_ptr<RapMapSAIndex<int32_t>> SAIdxPtr{nullptr};
  std::unique_ptr<RapMapSAIndex<int64_t>> BigSAIdxPtr{nullptr};

  if (h.bigSA()) {
    std::cerr << "Loading 64-bit suffix array index: \n";
    BigSAIdxPtr.reset(new RapMapSAIndex<int64_t>);
	  BigSAIdxPtr->load(indexPrefix);
  } else {
    std::cerr << "Loading 32-bit suffix array index: \n";
    SAIdxPtr.reset(new RapMapSAIndex<int32_t>);
	  SAIdxPtr->load(indexPrefix);
  }

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

	// Must be a power of 2
	size_t queueSize{268435456};
	spdlog::set_async_mode(queueSize);
	auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(outStream);
	std::shared_ptr<spdlog::logger> outLog = std::make_shared<spdlog::logger>("outLog", outputSink);
	outLog->set_pattern("%v");

	uint32_t nthread = numThreads.getValue();
	std::unique_ptr<paired_parser> pairParserPtr{nullptr};
	std::unique_ptr<single_parser> singleParserPtr{nullptr};

	if (!noout.getValue()) {
    if (h.bigSA()) {
      rapmap::utils::writeSAMHeader(*BigSAIdxPtr, outLog);
    } else {
      rapmap::utils::writeSAMHeader(*SAIdxPtr, outLog);
    }
	}

    bool strictCheck = strict.getValue();
    bool fuzzyIntersection = fuzzy.getValue();
	SpinLockT iomutex;
	    ScopedTimer timer;
	    HitCounters hctrs;
	    consoleLog->info("mapping reads . . . \n\n\n");
        if (pairedEnd) {
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

            if (h.bigSA()) {
              spawnProcessReadsThreads(nthread, pairParserPtr.get(), *BigSAIdxPtr, iomutex,
                outLog, hctrs, maxNumHits.getValue(), noout.getValue(), strictCheck, fuzzyIntersection, align.getValue());
            } else {
              spawnProcessReadsThreads(nthread, pairParserPtr.get(), *SAIdxPtr, iomutex,
                outLog, hctrs, maxNumHits.getValue(), noout.getValue(), strictCheck, fuzzyIntersection, align.getValue());
            }
            hctrs.totHits += hctrs.peHits + hctrs.seHits;
            delete [] pairFileList;
        } else {
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
            if (h.bigSA()) {
              spawnProcessReadsThreads(nthread, singleParserPtr.get(), *BigSAIdxPtr, iomutex,
                                      outLog, hctrs, maxNumHits.getValue(), noout.getValue(), strictCheck, align.getValue());
            } else {
              spawnProcessReadsThreads(nthread, singleParserPtr.get(), *SAIdxPtr, iomutex,
                                      outLog, hctrs, maxNumHits.getValue(), noout.getValue(), strictCheck, align.getValue());
            }
        }
	std::cerr << "\n\n";


    consoleLog->info("Done mapping reads.");
    consoleLog->info("In total saw {} reads.", hctrs.numReads);
    consoleLog->info("Final # hits per read = {}", hctrs.totHits / static_cast<float>(hctrs.numReads));
	consoleLog->info("flushing output queue.");
	outLog->flush();
	/*
	    consoleLog->info("Discarded {} reads because they had > {} alignments",
		    hctrs.tooManyHits, maxNumHits.getValue());
		    */


	if (haveOutputFile) {
	    outFile.close();
	}
	return 0;
    } catch (TCLAP::ArgException& e) {
	consoleLog->error("Exception [{}] when parsing argument {}", e.error(), e.argId());
	return 1;
    }

}
