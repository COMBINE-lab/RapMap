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
//#include "jellyfish/mer_dna.hpp"

#include "tclap/CmdLine.h"

/*extern "C" {
#include "kseq.h"
}
*/

#include "stringpiece.h"
#include "RapMapUtils.hpp"
#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
//#include "PairSequenceParser.hpp"
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
#include "SelectiveAlignment.hpp"
#include "SACollectorPair.hpp"
#include "HitManagerSA.hpp"

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
//using SelectiveAlignmentFormat = rapmap::utils::SelectiveAlignmentFormat;
using FixedWriter = rapmap::utils::FixedWriter;

struct MappingOpts {
    std::string index;
    std::string read1;
    std::string read2;
    std::string unmatedReads;
    uint32_t numThreads{1};
    uint32_t maxNumHits{200};
    std::string outname;
    double quasiCov{0.0};
    bool pairedEnd{false};
    bool noOutput{true};
    bool sensitive{false};
    bool strictCheck{false};
    bool fuzzy{false};
    bool consistentHits{false};
    bool quiet{false};
    bool remap{false};
    uint32_t mmpThreshold{15};
    int32_t editThreshold{10};
    uint32_t maxInsertSize{1000}; 
};

template <typename RapMapIndexT, typename MutexT>
void processReadsSingleSA(single_parser * parser,
                          RapMapIndexT& rmi,
                          MutexT* iomutex,
                          std::shared_ptr<spdlog::logger> outQueue,
                          HitCounters& hctr,
                          MappingOpts* mopts) {
    using OffsetT = typename RapMapIndexT::IndexType;
    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT> ;

    //Normal Call to SACollector
    SACollector<RapMapIndexT, rapmap::utils::my_mer> hitCollector(&rmi);
    //Remap caqll to SACollector
    //SACollector<RapMapIndexT, rapmap::utils::my_mer9> hitCollector9(&rmi);

    SACollectorPair<RapMapIndexT, rapmap::utils::my_mer> hitCollectorPair(&rmi);

    SECollector<RapMapIndexT> hitSECollector(&rmi);

    if (mopts->sensitive) {
        hitCollector.disableNIP();
    }
    hitCollector.setStrictCheck(mopts->strictCheck);
    if (mopts->quasiCov > 0.0) {
        hitCollector.setCoverageRequirement(mopts->quasiCov);
    }

    auto& txpNames = rmi.txpNames;
    auto& txpLens = rmi.txpLens;
    uint32_t n{0};

    auto logger = spdlog::get("stderrLog");

    fmt::MemoryWriter sstream;
    size_t batchSize{2500};

    std::vector<QuasiAlignment> hits;
    //std::vector<QuasiAlignment> selectedHits ;

    size_t readLen{0};
	bool tooManyHits{false};
    uint16_t flags;

    SingleAlignmentFormatter<RapMapIndexT*> formatter(&rmi);

    SASearcher<RapMapIndexT> saSearcher(&rmi);

    std::vector<uint32_t> dummy ;

    uint32_t orphanStatus{0};
    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();


    std::vector<SAIntervalHit> leftFwdSAInts;
    std::vector<SAIntervalHit> leftRcSAInts;

    while (parser->refill(rg)) {
      //while(true) {
      //  typename single_parser::job j(*parser); // Get a job from the parser: a bunch of reads (at most max_read_group)
      //  if(j.is_empty()) break;                 // If we got nothing, then quit.
      //  for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
      for (auto& read : rg) {
	    readLen = read.seq.length();//j->data[i].seq.length();
            ++hctr.numReads;
            hits.clear();
            //hitCollector(read.seq, hits, saSearcher, MateStatus::SINGLE_END, dummy, false, mopts->mmpThreshold, mopts->consistentHits);

            leftFwdSAInts.clear();
            leftRcSAInts.clear();

            bool lh = hitCollectorPair(read.seq,
                    leftFwdSAInts, leftRcSAInts, saSearcher,
                    MateStatus::SINGLE_END,
                    dummy,
                    false,
                    mopts->mmpThreshold,
                    mopts->consistentHits);


 	    bool res = rapmap::hit_manager_sa::mergeSAInts(
                    read,
                    lh,
                    leftFwdSAInts, leftRcSAInts,
                    hits,
                    rmi,
                    false,
                    mopts->consistentHits,
                    hctr
                    );


            /*if(hits.size() == 0 && mopts->remap){
                hitCollector9(read.seq, hits, saSearcher, MateStatus::SINGLE_END, dummy, true, mopts->mmpThreshold, mopts->consistentHits);
            }*/
            // @hirak
            // Here I collected all the QuasiALignments in
            // QuasiAlignment Object vector hits
            // which right now contains lcpLengths too
            // time to call the new function
            hitSECollector(read, hits, mopts->editThreshold);
                hits.erase(std::remove_if(hits.begin(), hits.end(),
                      [](QuasiAlignment& a) {
                      return !a.toAlign;
                      }), hits.end());



            auto numHits = hits.size();
            hctr.totHits += numHits;

	    if (hits.size() > 0 and !mopts->noOutput and hits.size() <= mopts->maxNumHits) {
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
                if (!mopts->quiet and iomutex->try_lock()){
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
        if (!mopts->noOutput) {
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
                        MutexT* iomutex,
                        std::shared_ptr<spdlog::logger> outQueue,
                        HitCounters& hctr,
                        MappingOpts* mopts) {
    using OffsetT = typename RapMapIndexT::IndexType;
    using SAIntervalHit = rapmap::utils::SAIntervalHit<OffsetT> ;

    //Normal Call to SACollector
    SACollector<RapMapIndexT, rapmap::utils::my_mer> hitCollector(&rmi);
    //Remap caqll to SACollector
    //SACollector<RapMapIndexT, rapmap::utils::my_mer9> hitCollector9(&rmi);
    //Rammap call to SACollectorPair
    SACollectorPair<RapMapIndexT, rapmap::utils::my_mer> hitCollectorPair(&rmi);

    SECollector<RapMapIndexT> hitSECollector(&rmi);

    if (mopts->sensitive) {
        hitCollector.disableNIP();
    }
    hitCollector.setStrictCheck(mopts->strictCheck);
    if (mopts->quasiCov > 0.0) {
        hitCollector.setCoverageRequirement(mopts->quasiCov);
    }

    auto& txpNames = rmi.txpNames;
    auto& txpLens = rmi.txpLens;
    uint32_t n{0};

    auto logger = spdlog::get("stderrLog");

    fmt::MemoryWriter sstream;
    size_t batchSize{1000};
    std::vector<QuasiAlignment> leftHits;
    std::vector<QuasiAlignment> rightHits;
    std::vector<QuasiAlignment> jointHits;

    //create things for SAIntervalHit s
    std::vector<SAIntervalHit> leftFwdSAInts ;
    std::vector<SAIntervalHit> leftRcSAInts ;
    std::vector<SAIntervalHit> rightFwdSAInts ;
    std::vector<SAIntervalHit> rightRcSAInts ;



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

    std::vector<uint32_t> dummy ;
    //auto txpNames = formatter.index->txpNames ;

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

            leftFwdSAInts.clear();
            leftRcSAInts.clear();
            rightFwdSAInts.clear();
            rightRcSAInts.clear();



            int32_t minLDist{mopts->editThreshold};
            int32_t minRDist{mopts->editThreshold};

            bool lhp = hitCollectorPair(rpair.first.seq,
                                   leftFwdSAInts, leftRcSAInts, saSearcher,
                                   MateStatus::PAIRED_END_LEFT,
                                   dummy,
                                   false,
                                   mopts->mmpThreshold,
                                   mopts->consistentHits);

            //comment out the stuff for now
            /*
            bool lh = hitCollector(rpair.first.seq,
                                   leftHits, saSearcher,
                                   MateStatus::PAIRED_END_LEFT,
                                   dummy,
                                   false,
                                   mopts->mmpThreshold,
                                   mopts->consistentHits);


            //call to hitCollectorPair

            if(leftHits.size() == 0 && mopts->remap){
                //std::cout <<"Before remap "  <<leftHits.size() << "\n";
                hitCollector9(rpair.first.seq,
                                   leftHits, saSearcher,
                                   MateStatus::PAIRED_END_LEFT,
                                   true,
                                   mopts->mmpThreshold,
                                   mopts->consistentHits);
                //std::cout << leftHits.size() << "\n";
            }


             if(leftHits.size() > 0) {
                  hitSECollector(rpair.first, leftHits, mopts->editThreshold);

            //@debug purpose
            auto readName = rapmap::utils::getReadName(rpair.first) ;


            if(readName == "read272375/ENST00000373816.4;mate1:2693-2792;mate2:2873-2972"){
                for(auto leftHitsIt = leftHits.begin(); leftHitsIt != leftHits.end(); ++leftHitsIt)
                    std::cout << "\nLeft After one call to SECollector "<< txpNames[leftHitsIt->tid] <<" edit distance: "<< leftHitsIt->editD<<" pos: "<< leftHitsIt->pos <<" Length: " << txpLens[leftHitsIt->tid] << "\n" ;
            }

            //end @debug purpose


                  leftHits.erase(std::remove_if(leftHits.begin(), leftHits.end(),
                              [](QuasiAlignment& a) {
                              return !a.toAlign;
                              }), leftHits.end());


                  auto old_size = leftHits.size() ;
                  if(leftHits.size() > 0){
                        std::for_each(leftHits.begin(), leftHits.end(),
                                [&minLDist](QuasiAlignment& a) {
                                               if (a.editD < minLDist) { minLDist = a.editD; }
                                              });
                  leftHits.erase(std::remove_if(leftHits.begin(), leftHits.end(),
                              [&minLDist](QuasiAlignment& a) {
                              return a.editD > minLDist;
                              }), leftHits.end());
                  }

                  auto new_size = leftHits.size() ;
                  if(old_size != 0 and new_size == 0) std::cout << "WTF \n" ;

             }
             */
            //comment out the right hits old mechanism
            bool rhp = hitCollectorPair(rpair.second.seq,
                                   rightFwdSAInts, rightRcSAInts, saSearcher,
                                   MateStatus::PAIRED_END_RIGHT,
                                   dummy,
                                   false,
                                   mopts->mmpThreshold,
                                   mopts->consistentHits);


            /*
            bool rh = hitCollector(rpair.second.seq,
                                   rightHits, saSearcher,
                                   MateStatus::PAIRED_END_RIGHT,
                                   dummy,
                                   false,
                                   mopts->mmpThreshold,
                                   mopts->consistentHits);


            if(rightHits.size() == 0 && mopts->remap){

                hitCollector9(rpair.second.seq,
                                   rightHits, saSearcher,
                                   MateStatus::PAIRED_END_RIGHT,
                                   true,
                                   mopts->mmpThreshold,
                                   mopts->consistentHits);
            }




           if(rightHits.size() > 0){

              hitSECollector(rpair.second, rightHits, mopts->editThreshold);


            //@debug purpose
            auto readName = rapmap::utils::getReadName(rpair.first) ;


            if(readName == "read272375/ENST00000373816.4;mate1:2693-2792;mate2:2873-2972"){
                for(auto leftHitsIt = rightHits.begin(); leftHitsIt != rightHits.end(); ++leftHitsIt)
                    std::cout << "\nRight After one call to SECollector "<< txpNames[leftHitsIt->tid] <<" edit distance: "<< leftHitsIt->editD<<" pos: "<< leftHitsIt->pos <<" Length: " << txpLens[leftHitsIt->tid] << "\n" ;
            }

              rightHits.erase(std::remove_if(rightHits.begin(), rightHits.end(),
                          [](QuasiAlignment& a) {
                          return !a.toAlign;
                          }), rightHits.end());

              auto old_size = rightHits.size();

              if(rightHits.size() > 0){
                  std::for_each(rightHits.begin(), rightHits.end(),
                            [&minRDist](QuasiAlignment& a) {
                                           if (a.editD < minRDist) { minRDist = a.editD; }
                                          });

                  rightHits.erase(std::remove_if(rightHits.begin(), rightHits.end(),
                          [&minRDist](QuasiAlignment& a) {
                          return a.editD > minRDist;
                          }), rightHits.end());
              }
              auto new_size = rightHits.size() ;

              if(old_size != 0 and new_size == 0) std::cout <<minRDist <<" WTF\n" ;
            }
           */

            //Now I have both forward and reverse checks from both the
            //ends, we can call some function from here which merges
            //them with respect to
            bool res = rapmap::hit_manager_sa::mergeLeftRightSAInts(
                    rpair,
                    lhp, rhp,
                    leftFwdSAInts, leftRcSAInts,
                    rightFwdSAInts, rightRcSAInts,
                    jointHits,
                    rmi,
                    mopts->maxNumHits,
                    mopts->consistentHits,
                    hctr,
                    mopts->editThreshold, 
		    mopts->maxInsertSize);

	   hitSECollector(rpair.first,rpair.second, jointHits, mopts->editThreshold);

           jointHits.erase(std::remove_if(jointHits.begin(), jointHits.end(),
                      [](QuasiAlignment& a) {
                      return !a.toAlign;
                      }), jointHits.end());

            /*
            if (mopts->fuzzy) {
            AI
                        lh, rh,
                        leftHits, rightHits, jointHits,
                        readLen, mopts->maxNumHits, tooManyHits, hctr);

            } else {
                rapmap::utils::mergeLeftRightHits(
                        leftHits, rightHits, jointHits,
                        readLen, mopts->maxNumHits, tooManyHits, hctr);
            }*/


            /*
            if(mopts->remap && jointHits.size() > 0){
                auto& startHit = jointHits.front();
                if(startHit.mateStatus == rapmap::utils::MateStatus::PAIRED_END_LEFT){
                    std::vector<uint32_t> goldenTids ;
                    for(auto& qa : jointHits){
                        goldenTids.push_back(qa.tid);
                    }
                rightHits.clear();
                hitCollector9(rpair.second.seq,
                                   rightHits, saSearcher,
                                   MateStatus::PAIRED_END_RIGHT,
                                   goldenTids,
                                   true,
                                   mopts->mmpThreshold,
                                   mopts->consistentHits);
                    if(rightHits.size() > 0){
                        jointHits.clear();
                        rapmap::utils::mergeLeftRightHitsFuzzy(
                                lh, rh,
                                leftHits, rightHits, jointHits,
                                readLen, mopts->maxNumHits, tooManyHits, hctr);
                    }
                }else if(startHit.mateStatus == rapmap::utils::MateStatus::PAIRED_END_RIGHT){
                    std::vector<uint32_t> goldenTids ;
                    for(auto& qa : jointHits){
                        goldenTids.push_back(qa.tid);
                    }
                    leftHits.clear();
                hitCollector9(rpair.first.seq,
                                   leftHits, saSearcher,
                                   MateStatus::PAIRED_END_LEFT,
                                   goldenTids,
                                   true,
                                   mopts->mmpThreshold,
                                   mopts->consistentHits);
                    if(leftHits.size() > 0){
                        jointHits.clear();
                        rapmap::utils::mergeLeftRightHitsFuzzy(
                                lh, rh,
                                leftHits, rightHits, jointHits,
                                readLen, mopts->maxNumHits, tooManyHits, hctr);
                    }
                }
            }*/

            hctr.totHits += jointHits.size();

            /** NOTE : the following line only holds until we fix orphans **/
            hctr.peHits += jointHits.size();

            // If we have reads to output, and we're writing output.
            //if (jointHits.size() > 0 and !mopts->noOutput and jointHits.size() <= mopts->maxNumHits) {
            if (jointHits.size() > 0 and !mopts->noOutput and jointHits.size() <= mopts->maxNumHits) {
                rapmap::utils::writeAlignmentsToStream(rpair, formatter,
                                                       hctr, jointHits, sstream);
            }

            if (hctr.numReads > hctr.lastPrint + 1000000) {
        		hctr.lastPrint.store(hctr.numReads.load());
                if (!mopts->quiet and iomutex->try_lock()) {
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
        if (!mopts->noOutput) {
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
                              MappingOpts* mopts) {

            std::vector<std::thread> threads;
            /*
            std::ofstream lcpLengthDumpFile;
            lcpLengthDumpFile.open("~/Projects/RapMap/build/lcpLength.list");
            //@debug with khash
            auto& khash = rmi.khash ;
            for(auto const &el : khash){
                lcpLengthDumpFile<<el.second.lcpLength << "\n";
            }
            lcpLengthDumpFile.close();
            */

            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsPairSA<RapMapIndexT, MutexT>,
                                     parser,
                                     std::ref(rmi),
                                     &iomutex,
                                     outQueue,
                                     std::ref(hctr),
                                     mopts);
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
                              MappingOpts* mopts) {
            std::vector<std::thread> threads;
            for (size_t i = 0; i < nthread; ++i) {
                threads.emplace_back(processReadsSingleSA<RapMapIndexT, MutexT>,
                                     parser,
                                     std::ref(rmi),
                                     &iomutex,
                                     outQueue,
                                     std::ref(hctr),
                                     mopts);
            }
            for (auto& t : threads) { t.join(); }
            return true;
        }

template <typename RapMapIndexT>
bool mapReads(RapMapIndexT& rmi,
	      std::shared_ptr<spdlog::logger> consoleLog,
          MappingOpts* mopts) {
	if (!mopts->quiet) { std::cerr << "\n\n\n\n"; }

	bool pairedEnd = mopts->pairedEnd;//(read1.isSet() or read2.isSet());
	// from: http://stackoverflow.com/questions/366955/obtain-a-stdostream-either-from-stdcout-or-stdofstreamfile
	// set either a file or cout as the output stream
	std::streambuf* outBuf;
	std::ofstream outFile;
	bool haveOutputFile{false};
	if (mopts->outname == "") {
	    outBuf = std::cout.rdbuf();
	} else {
	    outFile.open(mopts->outname);
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
	std::shared_ptr<spdlog::logger> outLog = std::make_shared<spdlog::logger>("rapmap::outLog", outputSink);
	outLog->set_pattern("%v");

	uint32_t nthread = mopts->numThreads;
	std::unique_ptr<paired_parser> pairParserPtr{nullptr};
	std::unique_ptr<single_parser> singleParserPtr{nullptr};

	if (!mopts->noOutput) {
	  rapmap::utils::writeSAMHeader(rmi, outLog);
	}

    //for the parser
    size_t chunkSize{10000};
	SpinLockT iomutex;
	{
	    ScopedTimer timer(!mopts->quiet);
	    HitCounters hctrs;
	    consoleLog->info("mapping reads . . . \n\n\n");
        if (pairedEnd) {
            std::vector<std::string> read1Vec = rapmap::utils::tokenize(mopts->read1, ',');
            std::vector<std::string> read2Vec = rapmap::utils::tokenize(mopts->read2, ',');

            if (read1Vec.size() != read2Vec.size()) {
                consoleLog->error("The number of provided files for "
                                  "-1 and -2 must be the same!");
                std::exit(1);
            }

	    uint32_t nprod = (read1Vec.size() > 1) ? 2 : 1;
	    pairParserPtr.reset(new paired_parser(read1Vec, read2Vec, nthread, nprod, chunkSize));
	    pairParserPtr->start();
            spawnProcessReadsThreads(nthread, pairParserPtr.get(), rmi, iomutex,
                                     outLog, hctrs, mopts);
        } else {
            std::vector<std::string> unmatedReadVec = rapmap::utils::tokenize(mopts->unmatedReads, ',');


	    uint32_t nprod = (unmatedReadVec.size() > 1) ? 2 : 1;
	    singleParserPtr.reset(new single_parser(unmatedReadVec, nthread, nprod, chunkSize));
	    singleParserPtr->start();
            /** Create the threads depending on the collector type **/
            spawnProcessReadsThreads(nthread, singleParserPtr.get(), rmi, iomutex,
                                      outLog, hctrs, mopts);
        }
	if (!mopts->quiet) { std::cerr << "\n\n"; }


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

void displayOpts(MappingOpts& mopts, spdlog::logger* log) {
        fmt::MemoryWriter optWriter;
        optWriter.write("\ncommand line options\n"
                        "====================\n");
        optWriter.write("index: {}\n", mopts.index);
        if (mopts.pairedEnd) {
            optWriter.write("read(s) 1: {}\n", mopts.read1);
            optWriter.write("read(s) 2: {}\n", mopts.read2);
        } else {
            optWriter.write("unmated read(s): {}\n", mopts.unmatedReads);
        }
        optWriter.write("output: {}\n", mopts.outname);
        optWriter.write("num. threads: {}\n", mopts.numThreads);
        optWriter.write("max num. hits: {}\n", mopts.maxNumHits);
        optWriter.write("quasi-coverage: {}\n", mopts.quasiCov);
        optWriter.write("no output: {}\n", mopts.noOutput);
        optWriter.write("sensitive: {}\n", mopts.sensitive);
        optWriter.write("strict check: {}\n", mopts.strictCheck);
        optWriter.write("fuzzy intersection: {}\n", mopts.fuzzy);
        optWriter.write("consistent hits: {}\n", mopts.consistentHits);
        optWriter.write("====================");
        log->info(optWriter.str());
}


int rapMapSAMap(int argc, char* argv[]) {
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
  TCLAP::ValueArg<double> quasiCov("z", "quasiCoverage", "Require that this fraction of a read is covered by MMPs before it is considered mappable.", false, 0.0, "double in [0,1]");
  TCLAP::SwitchArg noout("n", "noOutput", "Don't write out any alignments (for speed testing purposes)", false);
  TCLAP::SwitchArg sensitive("e", "sensitive", "Perform a more sensitive quasi-mapping by disabling NIP skipping", false);
  TCLAP::SwitchArg strict("s", "strictCheck", "Perform extra checks to try and assure that only equally \"best\" mappings for a read are reported", false);
  TCLAP::SwitchArg fuzzy("f", "fuzzyIntersection", "Find paired-end mapping locations using fuzzy intersection", false);
  TCLAP::SwitchArg consistent("c", "consistentHits", "Ensure that the hits collected are consistent (co-linear)", false);
  TCLAP::SwitchArg quiet("q", "quiet", "Disable all console output apart from warnings and errors", false);
  //flag if you want remap
  TCLAP::SwitchArg remap("g", "remap", "In case of unmapped reads/mates try again with 9 mer hash and recover", false);
  TCLAP::ValueArg<uint32_t> mmpThreshold("M", "mmpLen", "In case of unmapped reads/mates try again with 9 mer hash and recover", false, 15, "positive integer");
  TCLAP::ValueArg<int32_t> editThreshold("E", "editThreshold", "Edit score to allow", false, 10, "positive integer");
  TCLAP::ValueArg<uint32_t> maxInsertSize("I", "maxInsertSize", "Insert size to allow", false, 1000, "positive integer");

  cmd.add(index);
  cmd.add(noout);

  cmd.add(read1);
  cmd.add(read2);
  cmd.add(unmatedReads);
  cmd.add(outname);
  cmd.add(numThreads);
  cmd.add(maxNumHits);
  cmd.add(quasiCov);
  cmd.add(sensitive);
  cmd.add(strict);
  cmd.add(fuzzy);
  cmd.add(consistent);
  cmd.add(quiet);
  cmd.add(remap);
  cmd.add(mmpThreshold);
  cmd.add(editThreshold);
  cmd.add(maxInsertSize);

  auto rawConsoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
  auto consoleSink =
      std::make_shared<spdlog::sinks::ansicolor_sink>(rawConsoleSink);
  auto consoleLog = spdlog::create("stderrLog", {consoleSink});

  try {

    cmd.parse(argc, argv);
    // If we're supposed to be quiet, only print out warnings and above
    if (quiet.getValue()) {
        consoleLog->set_level(spdlog::level::warn);
    }

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

    MappingOpts mopts;
    if (pairedEnd) {
        mopts.read1 = read1.getValue();
        mopts.read2 = read2.getValue();
        mopts.pairedEnd = true;
    } else {
        mopts.unmatedReads = unmatedReads.getValue();
    }
    mopts.numThreads = numThreads.getValue();
    mopts.maxNumHits = maxNumHits.getValue();
    mopts.outname = (outname.isSet()) ? outname.getValue() : "";
    mopts.quasiCov = quasiCov.getValue();
    mopts.noOutput = noout.getValue();
    mopts.sensitive = sensitive.getValue();
    mopts.strictCheck = strict.getValue();
    mopts.consistentHits = consistent.getValue();
    mopts.fuzzy = fuzzy.getValue();
    mopts.quiet = quiet.getValue();
    mopts.remap = remap.getValue();
    mopts.mmpThreshold = mmpThreshold.getValue();
    mopts.editThreshold = editThreshold.getValue();
    mopts.maxInsertSize = maxInsertSize.getValue();

    if (quasiCov.isSet() and !sensitive.isSet()) {
        consoleLog->info("The --quasiCoverage option is set to {}, but the --sensitive flag was not set. The former implies the later. Enabling sensitive mode.", quasiCov.getValue());
        mopts.sensitive = true;
    }

    displayOpts(mopts, consoleLog.get());

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
          //RapMapSAIndex<int64_t, PerfectHashT<uint64_t, rapmap::utils::SAInterval<int64_t>>> rmi;
          RapMapSAIndex<int64_t, PerfectHashT<uint64_t, rapmap::utils::kmerVal<int64_t>>> rmi;
          rmi.load(indexPrefix);
          success = mapReads(rmi, consoleLog, &mopts);
      } else {
    	  /*
          RapMapSAIndex<int64_t,
                        RegHashT<uint64_t, rapmap::utils::SAInterval<int64_t>,
                                               rapmap::utils::KmerKeyHasher>> rmi_1;
          */
          RapMapSAIndex<int64_t,
                        RegHashT<uint64_t, rapmap::utils::kmerVal<int64_t>,
                                               rapmap::utils::KmerKeyHasher>> rmi;
          rmi.load(indexPrefix);
          success = mapReads(rmi, consoleLog, &mopts);
      }
    } else {
        //std::cerr << "Loading 32-bit suffix array index: \n";
      //SAIdxPtr.reset(new RapMapSAIndex<int32_t>);
      //SAIdxPtr->load(indexPrefix, h.kmerLen());
        if (h.perfectHash()) {
            //RapMapSAIndex<int32_t, PerfectHashT<uint64_t, rapmap::utils::SAInterval<int32_t>>> rmi;
            RapMapSAIndex<int32_t, PerfectHashT<uint64_t, rapmap::utils::kmerVal<int32_t>>> rmi;
            rmi.load(indexPrefix);
            success = mapReads(rmi, consoleLog, &mopts);
        } else {
            RapMapSAIndex<int32_t,
                          RegHashT<uint64_t, rapmap::utils::kmerVal<int32_t>,
                                                 rapmap::utils::KmerKeyHasher>> rmi;
            rmi.load(indexPrefix);
            success = mapReads(rmi, consoleLog, &mopts);
        }
    }

    return success ? 0 : 1;
  } catch (TCLAP::ArgException& e) {
    consoleLog->error("Exception [{}] when parsing argument {}", e.error(), e.argId());
    return 1;
  }

}
