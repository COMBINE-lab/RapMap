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

/*extern "C" {
#include "kseq.h"
}
*/
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



template <typename CollectorT, typename MutexT>
void processReadsSingleSA(single_parser * parser,
        RapMapSAIndex& rmi,
    	CollectorT& hitCollector,
        MutexT* iomutex,
    	std::shared_ptr<spdlog::logger> outQueue,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput) {

    auto& txpNames = rmi.txpNames;
    std::vector<uint32_t>& txpOffsets = rmi.txpOffsets;
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

    SingleAlignmentFormatter<RapMapSAIndex*> formatter(&rmi);

    SASearcher saSearcher(&rmi);

    uint32_t orphanStatus{0};
    while(true) {
        typename single_parser::job j(*parser); // Get a job from the parser: a bunch of reads (at most max_read_group)
        if(j.is_empty()) break;                 // If we got nothing, then quit.
        for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
            readLen = j->data[i].seq.length();
            ++hctr.numReads;
            hits.clear();
            hitCollector(j->data[i].seq, hits, saSearcher, MateStatus::SINGLE_END);
            auto numHits = hits.size();
            hctr.totHits += numHits;

            if (hits.size() > 0 and hits.size() <= maxNumHits) {
                /*
                std::sort(hits.begin(), hits.end(),
                            [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                return a.tid < b.tid;
                            });
                */
                rapmap::utils::writeAlignmentsToStream(j->data[i], formatter,
                                                       hctr, hits, sstream);
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
template <typename CollectorT, typename MutexT>
void processReadsPairSA(paired_parser* parser,
        RapMapSAIndex& rmi,
    	CollectorT& hitCollector,
        MutexT* iomutex,
	    std::shared_ptr<spdlog::logger> outQueue,
        HitCounters& hctr,
        uint32_t maxNumHits,
        bool noOutput) {
    auto& txpNames = rmi.txpNames;
    std::vector<uint32_t>& txpOffsets = rmi.txpOffsets;
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
    PairAlignmentFormatter<RapMapSAIndex*> formatter(&rmi);

    SASearcher saSearcher(&rmi);

    uint32_t orphanStatus{0};
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
                        MateStatus::PAIRED_END_LEFT);
            bool rh = hitCollector(j->data[i].second.seq,
                        rightHits, saSearcher,
                        MateStatus::PAIRED_END_RIGHT);

            rapmap::utils::mergeLeftRightHits(
                    leftHits, rightHits, jointHits,
                    readLen, maxNumHits, tooManyHits, hctr);

            // If we have reads to output, and we're writing output.
            if (jointHits.size() > 0 and !noOutput and jointHits.size() <= maxNumHits) {
                rapmap::utils::writeAlignmentsToStream(j->data[i], formatter,
                                                       hctr, jointHits, sstream);
            }

            if (hctr.numReads > hctr.lastPrint + 1000000) {
        		hctr.lastPrint.store(hctr.numReads.load());
                if (iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
                        std::cerr << "\r\r";
                    }
                    std::cerr << "saw " << hctr.numReads << " reads : "
                              << "pe / read = " << hctr.peHits / static_cast<float>(hctr.numReads)
                              << " : se / read = " << hctr.seHits / static_cast<float>(hctr.numReads) << ' ';
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

	RapMapSAIndex rmi;
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
	    rapmap::utils::writeSAMHeader(rmi, outLog);
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

            SACollector saCollector(&rmi);
            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsPairSA<SACollector, SpinLockT>,
                        pairParserPtr.get(),
                        std::ref(rmi),
                        std::ref(saCollector),
                        &iomutex,
            			outLog,
                        std::ref(hctrs),
                        maxNumHits.getValue(),
                        noout.getValue());
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
            SACollector saCollector(&rmi);
            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsSingleSA<SACollector, SpinLockT>,
                        singleParserPtr.get(),
                        std::ref(rmi),
                        std::ref(saCollector),
                        &iomutex,
            			outLog,
                        std::ref(hctrs),
                        maxNumHits.getValue(),
                        noout.getValue());
            }
            for (auto& t : threads) { t.join(); }
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
bool verbose{false};
if (verbose) {
    if (lb == lbLeftFwd and ub == ubLeftFwd) {
	//std::cerr << "couldn't narrow interval at all!\n";
    } else {
	diffCount += 1;
	int lbNaive, ubNaive, blah;

	//std::cerr << "before naive search . . . ";
	std::tie(lbNaive, ubNaive, blah) =
	    saSearcher.extendSearch(0, SA.size(), 0, rb, readEndIt);
	//std::cerr << "after naive search\n";

	if (ubNaive - lbNaive > ubLeftFwd - lbLeftFwd and blah >= matchedLen) {//lbLeftFwd != lbNaive or ubLeftFwd != ubNaive) {
	    diffCount += 1;
	//}
	//if (ubNaive - lbNaive > ubLeftFwd - lbLeftFwd and !verbose) {
	    std::cerr << "narrowed interval from "
		<< "[ " << lb << ", " << ub << ") to "
		<< "[ "<< lbLeftFwd << ", " << ubLeftFwd << ")\n";

	    std::cerr << "mer is    : " << mer << "\n";
	    std::cerr << "string is : ";
	    for (auto it = rb; it != rb + matchedLen; ++it) {
		std::cerr << *it;
	    } std::cerr << "\n";
	    std::cerr << "read is   : " << read << "\n";

	    std::cerr << "narrowed len = " << matchedLen << "\n";
	    std::cerr << "naive len = " << blah << "\n";

	    std::cerr << "entries at narrowed interval:\n";
	    for (auto x = lbLeftFwd - 1; x < ubLeftFwd + 1; ++x) {
		std::cerr << "\tT[SA[" << x << "] = "
		    << rmi_->seq.substr(SA[x], matchedLen) << "\n";
	    }

	    std::cerr << "naive interval is [" << lbNaive << ", " << ubNaive << ")\n";

	    std::cerr << "entries at naive interval:\n";
	    for (auto x = lbNaive - 1; x < ubNaive + 1; ++x) {
		std::cerr << "\tT[SA[" << x << "] = "
		    << rmi_->seq.substr(SA[x], blah) << "\n";
	    }
	    std::exit(1);
	}
    }

}
*/
/* Old read output code
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
	// The left and right transcipt ids
	leftTxp = leftIt->tid;
	rightTxp = rightIt->tid;

	// They don't point to the same transcript
	if (leftTxp < rightTxp) {
		++leftIt;
	} else {

		// The left and right iterators point to the same transcript
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
*/





