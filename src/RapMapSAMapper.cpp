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
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"

// Jellyfish 2 include
#include "FastxParser.hpp"
#include "jellyfish/mer_dna.hpp"

#include "tclap/CmdLine.h"

/*extern "C" {
#include "kseq.h"
}
*/
#include "stringpiece.h"
#include "BooMap.hpp"
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

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;
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



template <typename RapMapIndexT, typename MutexT>
void processReadsSingleSA(single_parser * parser,
                          RapMapIndexT& rmi,
                          //CollectorT& hitCollector,
                          MutexT* iomutex,
                          std::shared_ptr<spdlog::logger> outQueue,
                          HitCounters& hctr,
                          uint32_t maxNumHits,
                          bool noOutput,
                          bool strictCheck,
                          bool consistentHits) {

    using OffsetT = typename RapMapIndexT::IndexType;

    SACollector<RapMapIndexT> hitCollector(&rmi);
    auto& txpNames = rmi.txpNames;
    auto& txpLens = rmi.txpLens;
    uint32_t n{0};

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
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();

    while (parser->refill(rg)) {
      //while(true) {
      //  typename single_parser::job j(*parser); // Get a job from the parser: a bunch of reads (at most max_read_group)
      //  if(j.is_empty()) break;                 // If we got nothing, then quit.
      //  for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
      for (auto& read : rg) {
	    readLen = read.seq.length();//j->data[i].seq.length();
            ++hctr.numReads;
            hits.clear();
            hitCollector(read.seq, hits, saSearcher, MateStatus::SINGLE_END, strictCheck, consistentHits);
            auto numHits = hits.size();
            hctr.totHits += numHits;

	    if (hits.size() > 0 and !noOutput and hits.size() <= maxNumHits) {
                /*
                std::sort(hits.begin(), hits.end(),
                            [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                return a.tid < b.tid;
                            });
                */
                rapmap::utils::writeAlignmentsToStream(read, formatter,
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
                outQueue->info(std::move(outStr));
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
template <typename RapMapIndexT, typename MutexT>
void processReadsPairSA(paired_parser* parser,
                        RapMapIndexT& rmi,
                        //CollectorT& hitCollector,
                        MutexT* iomutex,
                        std::shared_ptr<spdlog::logger> outQueue,
                        HitCounters& hctr,
                        uint32_t maxNumHits,
                        bool noOutput,
                        bool strictCheck,
                        bool nonStrictMerge,
                        bool consistentHits) {

    using OffsetT = typename RapMapIndexT::IndexType;

    SACollector<RapMapIndexT> hitCollector(&rmi);
    auto& txpNames = rmi.txpNames;
    auto& txpLens = rmi.txpLens;
    uint32_t n{0};

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

    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();

    while (parser->refill(rg)) {
      //while(true) {
      //typename paired_parser::job j(*parser); // Get a job from the parser: a bunch of reads (at most max_read_group)
      //if(j.is_empty()) break;                 // If we got nothing, quit
      //  for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
      for (auto& rpair : rg) {
	tooManyHits = false;
	    readLen = rpair.first.seq.length();
            ++hctr.numReads;
            jointHits.clear();
            leftHits.clear();
            rightHits.clear();

            bool lh = hitCollector(rpair.first.seq,
                                   leftHits, saSearcher,
                                   MateStatus::PAIRED_END_LEFT,
                                   strictCheck,
                                   consistentHits);

            bool rh = hitCollector(rpair.second.seq,
                                   rightHits, saSearcher,
                                   MateStatus::PAIRED_END_RIGHT,
                                   strictCheck,
                                   consistentHits);

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

            hctr.totHits += jointHits.size();

            // If we have reads to output, and we're writing output.
            if (jointHits.size() > 0 and !noOutput and jointHits.size() <= maxNumHits) {
                rapmap::utils::writeAlignmentsToStream(rpair, formatter,
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
                outQueue->info(std::move(outStr));
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
                              bool consistentHits) {

            std::vector<std::thread> threads;
            SACollector<RapMapIndexT> saCollector(&rmi);
            //saCollector.disableNIP();
            //saCollector.setCoverageRequirement(0.5);

            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsPairSA<RapMapIndexT, MutexT>,
                                     parser,
                                     std::ref(rmi),
                                     //std::ref(saCollector),
                                     &iomutex,
                                     outQueue,
                                     std::ref(hctr),
                                     maxNumHits,
                                     noOutput,
                                     strictCheck,
                                     fuzzy,
                                     consistentHits);
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
                              bool consistentHits) {

            std::vector<std::thread> threads;
            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsSingleSA<RapMapIndexT, MutexT>,
                                     parser,
                                     std::ref(rmi),
                                     //std::ref(saCollector),
                                     &iomutex,
                                     outQueue,
                                     std::ref(hctr),
                                     maxNumHits,
                                     noOutput,
                                     strictCheck,
                                     consistentHits);
            }
            for (auto& t : threads) { t.join(); }
            return true;
        }

template <typename RapMapIndexT>
bool mapReads(RapMapIndexT& rmi,
	      std::shared_ptr<spdlog::logger> consoleLog,
	      TCLAP::ValueArg<std::string>& index,
	      TCLAP::ValueArg<std::string>& read1,
	      TCLAP::ValueArg<std::string>& read2,
	      TCLAP::ValueArg<std::string>& unmatedReads,
	      TCLAP::ValueArg<uint32_t>& numThreads,
	      TCLAP::ValueArg<uint32_t>& maxNumHits,
	      TCLAP::ValueArg<std::string>& outname,
	      TCLAP::SwitchArg& noout,
	      TCLAP::SwitchArg& strict,
          TCLAP::SwitchArg& fuzzy,
          TCLAP::SwitchArg& consistent) {

	std::cerr << "\n\n\n\n";

	bool pairedEnd = (read1.isSet() or read2.isSet());
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

    bool strictCheck = strict.getValue();
    bool fuzzyIntersection = fuzzy.getValue();
    bool consistentHits = consistent.getValue();
	SpinLockT iomutex;
	{
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

	    uint32_t nprod = (read1Vec.size() > 1) ? 2 : 1; 
	    pairParserPtr.reset(new paired_parser(read1Vec, read2Vec, nthread, nprod));
	    pairParserPtr->start();
            spawnProcessReadsThreads(nthread, pairParserPtr.get(), rmi, iomutex,
                                     outLog, hctrs, maxNumHits.getValue(), noout.getValue(), strictCheck,
                                     fuzzyIntersection, consistentHits);
        } else {
            std::vector<std::string> unmatedReadVec = rapmap::utils::tokenize(unmatedReads.getValue(), ',');


	    uint32_t nprod = (unmatedReadVec.size() > 1) ? 2 : 1; 
	    singleParserPtr.reset(new single_parser(unmatedReadVec, nthread, nprod));
	    singleParserPtr->start();
            /** Create the threads depending on the collector type **/
            spawnProcessReadsThreads(nthread, singleParserPtr.get(), rmi, iomutex,
                                      outLog, hctrs, maxNumHits.getValue(), noout.getValue(),
                                     strictCheck, consistentHits);
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

  TCLAP::ValueArg<std::string> index("i", "index", "The location of the quasiindex", true, "", "path");
  TCLAP::ValueArg<std::string> read1("1", "leftMates", "The location of the left paired-end reads", false, "", "path");
  TCLAP::ValueArg<std::string> read2("2", "rightMates", "The location of the right paired-end reads", false, "", "path");
  TCLAP::ValueArg<std::string> unmatedReads("r", "unmatedReads", "The location of single-end reads", false, "", "path");
  TCLAP::ValueArg<uint32_t> numThreads("t", "numThreads", "Number of threads to use", false, 1, "positive integer");
  TCLAP::ValueArg<uint32_t> maxNumHits("m", "maxNumHits", "Reads mapping to more than this many loci are discarded", false, 200, "positive integer");
  TCLAP::ValueArg<std::string> outname("o", "output", "The output file (default: stdout)", false, "", "path");
  TCLAP::SwitchArg noout("n", "noOutput", "Don't write out any alignments (for speed testing purposes)", false);
  TCLAP::SwitchArg strict("s", "strictCheck", "Perform extra checks to try and assure that only equally \"best\" mappings for a read are reported", false);
  TCLAP::SwitchArg fuzzy("f", "fuzzyIntersection", "Find paired-end mapping locations using fuzzy intersection", false);
  TCLAP::SwitchArg consistent("c", "consistentHits", "Ensure that the hits collected are consistent (co-linear)", false);
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
  cmd.add(consistent);

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

    //std::unique_ptr<RapMapSAIndex<int32_t>> SAIdxPtr{nullptr};
    //std::unique_ptr<RapMapSAIndex<int64_t>> BigSAIdxPtr{nullptr};

    bool success{false};
    if (h.bigSA()) {
        //std::cerr << "Loading 64-bit suffix array index: \n";
      //BigSAIdxPtr.reset(new RapMapSAIndex<int64_t>);
      //BigSAIdxPtr->load(indexPrefix, h.kmerLen());
      if (h.perfectHash()) {
          RapMapSAIndex<int64_t, BooMap<uint64_t, rapmap::utils::SAInterval<int64_t>>> rmi;
          rmi.load(indexPrefix);
          success = mapReads(rmi, consoleLog, index, read1, read2,
                             unmatedReads, numThreads, maxNumHits,
                             outname, noout, strict, fuzzy, consistent);
      } else {
          RapMapSAIndex<int64_t,
                        RegHashT<uint64_t, rapmap::utils::SAInterval<int64_t>,
                                               rapmap::utils::KmerKeyHasher>> rmi;
          rmi.load(indexPrefix);
          success = mapReads(rmi, consoleLog, index, read1, read2,
                             unmatedReads, numThreads, maxNumHits,
                             outname, noout, strict, fuzzy, consistent);
      }
    } else {
        //std::cerr << "Loading 32-bit suffix array index: \n";
      //SAIdxPtr.reset(new RapMapSAIndex<int32_t>);
      //SAIdxPtr->load(indexPrefix, h.kmerLen());
        if (h.perfectHash()) {
            RapMapSAIndex<int32_t, BooMap<uint64_t, rapmap::utils::SAInterval<int32_t>>> rmi;
            rmi.load(indexPrefix);
            success = mapReads(rmi, consoleLog, index, read1, read2,
                               unmatedReads, numThreads, maxNumHits,
                               outname, noout, strict, fuzzy, consistent);
        } else {
            RapMapSAIndex<int32_t,
                          RegHashT<uint64_t, rapmap::utils::SAInterval<int32_t>,
                                                 rapmap::utils::KmerKeyHasher>> rmi;
            rmi.load(indexPrefix);
            success = mapReads(rmi, consoleLog, index, read1, read2,
                               unmatedReads, numThreads, maxNumHits,
                               outname, noout, strict, fuzzy, consistent);
        }
    }

    return success ? 0 : 1;
  } catch (TCLAP::ArgException& e) {
    consoleLog->error("Exception [{}] when parsing argument {}", e.error(), e.argId());
    return 1;
  }

}
