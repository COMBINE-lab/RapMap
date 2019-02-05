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

#include "FastxParser.hpp"

#include "tclap/CmdLine.h"
#include "zstr/zstr.hpp"

/*extern "C" {
#include "kseq.h"
}
*/

#include "stringpiece.h"
#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
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
#include "ksw2pp/KSW2Aligner.hpp"
#include "SelectiveAlignmentUtils.hpp"
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
  /** DEPRECATED --- should remain false
   *  until we can remove it from the codebase
   **/
  bool consistentHits{false};
 
  bool doChaining{false};
  bool compressedOutput{false};
  bool quiet{false};
  int32_t maxMMPExtension{7};
  bool selAln{false};
  float consensusSlack{0.0};
  double minScoreFraction{0.65};
  int16_t matchScore{2};
  int16_t mismatchPenalty{-4};
  int16_t gapOpenPenalty{4};
  int16_t gapExtendPenalty{2};
  int32_t dpBandwidth{15};
  //int32_t consensusSlack{0.2};
  bool hardFilter{false};
  selective_alignment::utils::AlignmentPolicy ap{selective_alignment::utils::AlignmentPolicy::DEFAULT};
  bool noOrphans{false};
  bool noDovetail{false};
  bool recoverOrphans{false};
  bool writeUnmapped{false};
};

using AlnCacheMap = selective_alignment::utils::AlnCacheMap;

template <typename RapMapIndexT, typename MutexT>
void processReadsSingleSA(single_parser * parser,
                          RapMapIndexT& rmi,
                          MutexT* iomutex,
                          std::shared_ptr<spdlog::logger> outQueue,
                          HitCounters& hctr,
                          MappingOpts* mopts) {
    using OffsetT = typename RapMapIndexT::IndexType;

    SACollector<RapMapIndexT> hitCollector(&rmi);
    if (mopts->sensitive) {
        hitCollector.disableNIP();
    }
    hitCollector.setStrictCheck(mopts->strictCheck);
    if (mopts->quasiCov > 0.0) {
        hitCollector.setCoverageRequirement(mopts->quasiCov);
    }

    auto logger = spdlog::get("stderrLog");

    fmt::MemoryWriter sstream;
    std::vector<QuasiAlignment> hits;

    rapmap::hit_manager::HitCollectorInfo<rapmap::utils::SAIntervalHit<OffsetT>> hcInfo;
    rapmap::utils::MappingConfig mc;
    mc.consistentHits = false;//mopts->consistentHits;
    mc.doChaining = mopts->selAln;
    if (mc.doChaining) {
      auto consensusSlack = mopts->consensusSlack;
      mc.consensusFraction = (consensusSlack == 0.0) ? 1.0 : (1.0 - consensusSlack);
      mc.considerMultiPos = true;
      hitCollector.enableChainScoring();
      hitCollector.setMaxMMPExtension(mopts->maxMMPExtension);
    }
    bool hardFilter = mopts->hardFilter;
    std::string rc1; rc1.reserve(300);

    // Start: Setup alignment engine
    using ksw2pp::KSW2Aligner;
    using ksw2pp::KSW2Config;
    using ksw2pp::EnumToType;
    using ksw2pp::KSW2AlignmentType;
    KSW2Config config;
    config.dropoff = -1;
    config.gapo = mopts->gapOpenPenalty;
    config.gape = mopts->gapExtendPenalty;
    config.bandwidth = mopts->dpBandwidth;
    config.flag = 0;
    config.flag |= KSW_EZ_SCORE_ONLY;
    int8_t a = static_cast<int8_t>(mopts->matchScore);
    int8_t b = static_cast<int8_t>(mopts->mismatchPenalty);
    KSW2Aligner aligner(static_cast<int8_t>(a), static_cast<int8_t>(b));
    aligner.config() = config;
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    auto ap = mopts->ap;

    size_t numMappingsDropped{0};
    size_t numFragsDropped{0};

    AlnCacheMap alnCache; alnCache.reserve(16);
    // Done: Setup alignment engine

    SingleAlignmentFormatter<RapMapIndexT*> formatter(&rmi);
    SASearcher<RapMapIndexT> saSearcher(&rmi);

    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();

    while (parser->refill(rg)) {
      for (auto& read : rg) {
            read.seq.length();
            ++hctr.numReads;
            hcInfo.clear();
            hits.clear();
            hitCollector(read.seq, saSearcher, hcInfo);

            rapmap::hit_manager::hitsToMappingsSimple(rmi, mc,
                                                      MateStatus::SINGLE_END,
                                                      hcInfo, hits);

            auto numHits = hits.size();
            hctr.totHits += numHits;

            if(hits.size() > mopts->maxNumHits) {
              hits.clear();
            }

            if (mopts->selAln) {
              alnCache.clear();
              auto* r1 = read.seq.data();
              auto l1 = static_cast<int32_t>(read.seq.length());

              char* r1rc = nullptr;
              int32_t bestScore{std::numeric_limits<int32_t>::lowest()};
              std::vector<decltype(bestScore)> scores(hits.size(), bestScore);
              size_t idx{0};
              double optFrac{mopts->minScoreFraction};
              int32_t maxReadScore = a * read.seq.length();
              bool multiMapping{hits.size() > 1};

              for (auto& h : hits) {
                int32_t score{std::numeric_limits<int32_t>::min()};
                auto txpID = h.tid;
                char* tseq = const_cast<char*>(rmi.seq.data()) + rmi.txpOffsets[txpID];
                const int32_t tlen = static_cast<int32_t>(rmi.txpLens[txpID]);
                const uint32_t buf{20};

                // compute the reverse complement only if we
                // need it and don't have it
                if (!h.fwd and !r1rc) {
                  rapmap::utils::reverseRead(read.seq, rc1);
                  // we will not break the const promise
                  r1rc = const_cast<char*>(rc1.data());
                }

                auto* rptr = h.fwd ? r1 : r1rc;
                int32_t s =
                  selective_alignment::utils::getAlnScore(aligner, ez, h.pos, rptr, l1, tseq, tlen, a, b, maxReadScore, h.chainStatus.getLeft(),
                                                          multiMapping, ap, buf, alnCache);
                if (s < (optFrac * maxReadScore)) {
                  score = std::numeric_limits<decltype(score)>::min();
                } else {
                  score = s;
                }
                bestScore = (score > bestScore) ? score : bestScore;
                scores[idx] = score;
                h.score(score);
                ++idx;
              }

              uint32_t ctr{0};
              if (bestScore > std::numeric_limits<int32_t>::min()) {
                // Note --- with soft filtering, only those hits that are given the minimum possible
                // score are filtered out.
                hits.erase(
                                std::remove_if(hits.begin(), hits.end(),
                                               [&ctr, &scores, &numMappingsDropped, bestScore, hardFilter] (const QuasiAlignment& qa) -> bool {
                                                 // if soft filtering, we only drop things with an invalid score
                                                 // if hard filtering, we drop everything with a sub-optimal score.
                                                 bool rem = hardFilter ? (scores[ctr] < bestScore) :
                                                   (scores[ctr] == std::numeric_limits<int32_t>::min());
                                                 ++ctr;
                                                 numMappingsDropped += rem ? 1 : 0;
                                                 return rem;
                                               }),
                                hits.end()
                                );
                // for soft filter
                double bestScoreD = static_cast<double>(bestScore);
                std::for_each(hits.begin(), hits.end(),
                              [bestScoreD, hardFilter](QuasiAlignment& qa) -> void {
                                qa.alnScore(static_cast<int32_t>(qa.score())); 
                                double v = bestScoreD - qa.score();
                                qa.score( (hardFilter ? -1.0 : std::exp(-v)) );
                              });
              } else {
                ++numFragsDropped;
                hits.clear();
              }
            }

            if (!mopts->noOutput){
              if (hits.size() > 0) {
                rapmap::utils::writeAlignmentsToStream(read, formatter,
                                                       hctr, hits, sstream);
              } else {
                rapmap::utils::writeUnalignedSingleToStream(read, sstream);
              }
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

    SACollector<RapMapIndexT> hitCollector(&rmi);
    if (mopts->sensitive) {
        hitCollector.disableNIP();
    }
    hitCollector.setStrictCheck(mopts->strictCheck);
    if (mopts->quasiCov > 0.0) {
        hitCollector.setCoverageRequirement(mopts->quasiCov);
    }

    auto logger = spdlog::get("stderrLog");

    fmt::MemoryWriter sstream;
    std::vector<QuasiAlignment> leftHits;
    std::vector<QuasiAlignment> rightHits;
    std::vector<QuasiAlignment> jointHits;
    std::string rc1; rc1.reserve(250);
    std::string rc2; rc2.reserve(250);

    size_t readLen{0};
    bool tooManyHits{false};

    rapmap::hit_manager::HitCollectorInfo<rapmap::utils::SAIntervalHit<OffsetT>> leftHCInfo;
    rapmap::hit_manager::HitCollectorInfo<rapmap::utils::SAIntervalHit<OffsetT>> rightHCInfo;
    rapmap::utils::MappingConfig mc;
    mc.consistentHits = false;//mopts->consistentHits;
    mc.doChaining = mopts->selAln;
    if (mc.doChaining) {
      auto consensusSlack = mopts->consensusSlack;
      mc.consensusFraction = (consensusSlack == 0.0) ? 1.0 : (1.0 - consensusSlack);
      mc.considerMultiPos = true;
      hitCollector.enableChainScoring();
      hitCollector.setMaxMMPExtension(mopts->maxMMPExtension);
    }
    bool hardFilter = mopts->hardFilter;

    // Start: Setup alignment engine
    using ksw2pp::KSW2Aligner;
    using ksw2pp::KSW2Config;
    using ksw2pp::EnumToType;
    using ksw2pp::KSW2AlignmentType;
    KSW2Config config;
    config.dropoff = -1;
    config.gapo = mopts->gapOpenPenalty;
    config.gape = mopts->gapExtendPenalty;
    config.bandwidth = mopts->dpBandwidth;
    config.flag = 0;
    config.flag |= KSW_EZ_SCORE_ONLY;
    int8_t a = static_cast<int8_t>(mopts->matchScore);
    int8_t b = static_cast<int8_t>(mopts->mismatchPenalty);
    KSW2Aligner aligner(static_cast<int8_t>(a), static_cast<int8_t>(b));
    aligner.config() = config;
    ksw_extz_t ez;
    memset(&ez, 0, sizeof(ksw_extz_t));
    bool noDovetail = mopts->noDovetail;
    bool noOrphans = mopts->noOrphans;
    size_t numOrphansRescued{0};
    auto ap = mopts->ap;

    size_t numMappingsDropped{0};
    size_t numFragsDropped{0};

    AlnCacheMap alnCacheLeft; alnCacheLeft.reserve(32);
    AlnCacheMap alnCacheRight; alnCacheRight.reserve(32);
    // Done: Setup alignment engine

    bool useSmartIntersect{mopts->fuzzy or mopts->selAln};

    // Create a formatter for alignments
    PairAlignmentFormatter<RapMapIndexT*> formatter(&rmi);

    SASearcher<RapMapIndexT> saSearcher(&rmi);

    // Get the read group by which this thread will
    // communicate with the parser (*once per-thread*)
    auto rg = parser->getReadGroup();

    while (parser->refill(rg)) {
      for (auto& rpair : rg) {
        tooManyHits = false;
        readLen = rpair.first.seq.length();
            ++hctr.numReads;
            leftHCInfo.clear();
            rightHCInfo.clear();
            jointHits.clear();
            leftHits.clear();
            rightHits.clear();

            bool lh = hitCollector(rpair.first.seq,
                                   saSearcher,
                                   leftHCInfo);

            bool rh = hitCollector(rpair.second.seq,
                                   saSearcher,
                                   rightHCInfo);

           rapmap::hit_manager::hitsToMappingsSimple(rmi, mc,
                                                     MateStatus::PAIRED_END_LEFT,
                                                     leftHCInfo, leftHits);

           rapmap::hit_manager::hitsToMappingsSimple(rmi, mc,
                                                     MateStatus::PAIRED_END_RIGHT,
                                                     rightHCInfo, rightHits);

           if (useSmartIntersect) {
             auto mergeRes = rapmap::utils::mergeLeftRightHitsFuzzy(
                                                                    lh, rh,
                                                                    leftHits, rightHits, jointHits,
                                                                    mc,
                                                                    readLen, mopts->maxNumHits, tooManyHits, hctr);
             bool mergeStatusOK = (mergeRes == rapmap::utils::MergeResult::HAD_EMPTY_INTERSECTION or
                                   mergeRes == rapmap::utils::MergeResult::HAD_ONLY_LEFT or
                                   mergeRes == rapmap::utils::MergeResult::HAD_ONLY_RIGHT);

             if ( mergeStatusOK and mopts->recoverOrphans and !tooManyHits) {
               if (leftHits.size() + rightHits.size() > 0) {
                 if (mergeRes == rapmap::utils::MergeResult::HAD_ONLY_LEFT) {
                   // In this case, we've "moved" the left hit's into joint, so put them back into
                   // left and make sure joint is clear.
                   leftHits.swap(jointHits);
                   jointHits.clear();
                 } else if (mergeRes == rapmap::utils::MergeResult::HAD_ONLY_RIGHT) {
                   // In this case, we've "moved" the right hit's into joint, so put them back into
                   // right and make sure joint is clear.
                   rightHits.swap(jointHits);
                   jointHits.clear();
                 }
                 selective_alignment::utils::recoverOrphans(rpair.first.seq,
                                                            rpair.second.seq,
                                                            rc1,
                                                            rc2,
                                                            //
                                                            rmi.seq,
                                                            rmi.txpNames,
                                                            rmi.txpOffsets,
                                                            rmi.txpLens,
                                                            //
                                                            leftHits,
                                                            rightHits,
                                                            jointHits);
                 if (!jointHits.empty()) { numOrphansRescued++; }
               }
             }
           } else {
             rapmap::utils::mergeLeftRightHits(
                                               leftHits, rightHits, jointHits,
                                               readLen, mopts->maxNumHits, tooManyHits, hctr);
           }

           // If the read mapped to > maxReadOccs places, discard it
           if (jointHits.size() > mopts->maxNumHits) {
             jointHits.clear();
           }

           // If we have mappings, then process them.
           if (!jointHits.empty()) {
             bool isPaired = jointHits.front().mateStatus ==
               rapmap::utils::MateStatus::PAIRED_END_PAIRED;
             // If we are ignoring orphans
             if (noOrphans) {
               // If the mappings for the current read are not properly-paired (i.e.
               // are orphans)
               // then just clear the group.
               if (!isPaired) {
                 jointHits.clear();
               }
             }
           }

           // Start: If requested, perform selective alignment
           if (mopts->selAln and !jointHits.empty()) {
             alnCacheLeft.clear();
             alnCacheRight.clear();
             auto* r1 = rpair.first.seq.data();
             auto* r2 = rpair.second.seq.data();
             auto l1 = static_cast<int32_t>(rpair.first.seq.length());
             auto l2 = static_cast<int32_t>(rpair.second.seq.length());
             // We compute the reverse complements below only if we
             // need them and don't have them.
             char* r1rc = nullptr;
             char* r2rc = nullptr;
             int32_t bestScore{std::numeric_limits<int32_t>::lowest()};
             std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
             size_t idx{0};
             double optFrac{mopts->minScoreFraction};
             int32_t maxLeftScore{a * static_cast<int32_t>(rpair.first.seq.length())};
             int32_t maxRightScore{a * static_cast<int32_t>(rpair.second.seq.length())};
             bool multiMapping{jointHits.size() > 1};

             for (auto& h : jointHits) {
               int32_t score{std::numeric_limits<int32_t>::min()};
               auto txpID = h.tid;
               char* tseq = const_cast<char*>(rmi.seq.data()) + rmi.txpOffsets[txpID];
               const int32_t tlen = static_cast<int32_t>(rmi.txpLens[txpID]);
               const uint32_t buf{20};

               if (h.mateStatus == rapmap::utils::MateStatus::PAIRED_END_PAIRED) {
                 if (!h.fwd and !r1rc) {
                   rapmap::utils::reverseRead(rpair.first.seq, rc1);
                   r1rc = const_cast<char*>(rc1.data());
                 }
                 if (!h.mateIsFwd and !r2rc) {
                   rapmap::utils::reverseRead(rpair.second.seq, rc2);
                   r2rc = const_cast<char*>(rc2.data());
                 }
                 auto* r1ptr = h.fwd ? r1 : r1rc;
                 auto* r2ptr = h.mateIsFwd ? r2 : r2rc;

                 int32_t s1 =
                   selective_alignment::utils::getAlnScore(aligner, ez, h.pos, r1ptr, l1, tseq, tlen, a, b, maxLeftScore, h.chainStatus.getLeft(),
                                                           multiMapping, ap, buf, alnCacheLeft);

                 int32_t s2 =
                   selective_alignment::utils::getAlnScore(aligner, ez, h.matePos, r2ptr, l2, tseq, tlen, a, b, maxRightScore, h.chainStatus.getRight(),
                                                           multiMapping, ap, buf, alnCacheRight);

                 // throw away dovetailed reads
                 if (h.fwd != h.mateIsFwd and noDovetail) {
                   if (h.fwd and (h.pos > h.matePos)) {
                     s1 = std::numeric_limits<int32_t>::min();
                     s2 = std::numeric_limits<int32_t>::min();
                   } else if (h.mateIsFwd and (h.matePos > h.pos)) {
                     s1 = std::numeric_limits<int32_t>::min();
                     s2 = std::numeric_limits<int32_t>::min();
                   }
                 }

                 // ends are scored separately
                 if ((s1 < (optFrac * maxLeftScore)) or (s2 < (optFrac * maxRightScore))) {
                   score = std::numeric_limits<decltype(score)>::min();
                 } else {
                   score = s1 + s2;
                 }
               } else if (h.mateStatus == rapmap::utils::MateStatus::PAIRED_END_LEFT) {
                 if (!h.fwd and !r1rc) {
                   rapmap::utils::reverseRead(rpair.first.seq, rc1);
                   r1rc = const_cast<char*>(rc1.data());
                 }
                 auto* rptr = h.fwd ? r1 : r1rc;

                 int32_t s =
                   selective_alignment::utils::getAlnScore(aligner, ez, h.pos, rptr, l1, tseq, tlen, a, b, maxLeftScore, h.chainStatus.getLeft(),
                                                           multiMapping, ap, buf, alnCacheLeft);
                 if (s < (optFrac * maxLeftScore)) {
                   score = std::numeric_limits<decltype(score)>::min();
                 } else {
                   score = s;
                 }
               } else if (h.mateStatus == rapmap::utils::MateStatus::PAIRED_END_RIGHT) {
                 if (!h.fwd and !r2rc) {
                   rapmap::utils::reverseRead(rpair.second.seq, rc2);
                   r2rc = const_cast<char*>(rc2.data());
                 }
                 auto* rptr = h.fwd ? r2 : r2rc;

                 int32_t s =
                   selective_alignment::utils::getAlnScore(aligner, ez, h.pos, rptr, l2, tseq, tlen, a, b, maxRightScore, h.chainStatus.getRight(),
                                                           multiMapping, ap, buf, alnCacheRight);
                 if (s < (optFrac * maxRightScore)) {
                   score = std::numeric_limits<decltype(score)>::min();
                 } else {
                   score = s;
                 }
               }

               bestScore = (score > bestScore) ? score : bestScore;
               scores[idx] = score;
               h.score(score);
               ++idx;
             }

             uint32_t ctr{0};
             if (bestScore > std::numeric_limits<int32_t>::min()) {
               // Note --- with soft filtering, only those hits that are given the minimum possible
               // score are filtered out.
               jointHits.erase(
                               std::remove_if(jointHits.begin(), jointHits.end(),
                                              [&ctr, &scores, &numMappingsDropped, bestScore, hardFilter] (const QuasiAlignment& qa) -> bool {
                                                // if soft filtering, we only drop things with an invalid score
                                                // if hard filtering, we drop everything with a sub-optimal score.
                                                bool rem = hardFilter ? (scores[ctr] < bestScore) :
                                                  (scores[ctr] == std::numeric_limits<int32_t>::min());
                                                ++ctr;
                                                numMappingsDropped += rem ? 1 : 0;
                                                return rem;
                                              }),
                               jointHits.end()
                               );

               double bestScoreD = static_cast<double>(bestScore);
               std::for_each(jointHits.begin(), jointHits.end(),
                             [bestScoreD, hardFilter](QuasiAlignment& qa) -> void {
                               qa.alnScore(static_cast<int32_t>(qa.score()));
                               double v = bestScoreD - qa.score();
                               qa.score( (hardFilter ? -1.0 : std::exp(-v)) );
                             });
             } else {
               ++numFragsDropped;
               jointHits.clear();
             }
           } else if (noDovetail) {
             jointHits.erase(
                             std::remove_if(jointHits.begin(), jointHits.end(),
                                            [](const QuasiAlignment& h) -> bool {
                                              if (h.fwd != h.mateIsFwd) {
                                                if (h.fwd and (h.pos > h.matePos)) {
                                                  return true;
                                                } else if (h.mateIsFwd and (h.matePos > h.pos)) {
                                                  return true;
                                                }
                                              }
                                              return false;
                                            }),
                             jointHits.end());
           }
           // Done: If requested, perform selective alignment

            hctr.totHits += jointHits.size();

            // If we have reads to output, and we're writing output.
            if (!mopts->noOutput) {
              if (jointHits.size() > 0 and jointHits.size() <= mopts->maxNumHits) {
                rapmap::utils::writeAlignmentsToStream(rpair, formatter,
                                                       hctr, jointHits, sstream);
              } else {
                rapmap::utils::writeUnalignedPairToStream(rpair, sstream);
              }
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
  const constexpr uint32_t buffSize{262144};
  std::vector<char> customStreamBuffer;
  customStreamBuffer.resize(buffSize);
	std::streambuf* outBuf;
	std::ofstream outFile;
  std::unique_ptr<std::ostream> outStream{nullptr};
	bool haveOutputFile{false};
  std::shared_ptr<spdlog::logger> outLog{nullptr};
  if (!mopts->noOutput) {
    if (mopts->outname == "") {
      outBuf = std::cout.rdbuf();
    } else {
      outBuf = outFile.rdbuf();
      outFile.open(mopts->outname);
      haveOutputFile = true;
    }
    // set the stream buffer size --- thanks for the suggestion @dnbaker!
    outBuf->pubsetbuf(customStreamBuffer.data(), customStreamBuffer.size());

    // Now set the output stream to the buffer, which is
    // either std::cout, or a file.

    if (mopts->compressedOutput) {
      outStream.reset(new zstr::ostream(outBuf));
    } else {
      outStream.reset(new std::ostream(outBuf));
    }
    //std::ostream outStream(outBuf);

    // Must be a power of 2
    size_t queueSize{268435456};
    spdlog::set_async_mode(queueSize);
    auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(*outStream);
    outLog = std::make_shared<spdlog::logger>("rapmap::outLog", outputSink);
    outLog->set_pattern("%v");

    rapmap::utils::writeSAMHeader(rmi, outLog);
  }

	uint32_t nthread = mopts->numThreads;
	std::unique_ptr<paired_parser> pairParserPtr{nullptr};
	std::unique_ptr<single_parser> singleParserPtr{nullptr};
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
      pairParserPtr->stop();
    } else {
      std::vector<std::string> unmatedReadVec = rapmap::utils::tokenize(mopts->unmatedReads, ',');


      uint32_t nprod = (unmatedReadVec.size() > 1) ? 2 : 1; 
      singleParserPtr.reset(new single_parser(unmatedReadVec, nthread, nprod, chunkSize));
      singleParserPtr->start();
      /** Create the threads depending on the collector type **/
      spawnProcessReadsThreads(nthread, singleParserPtr.get(), rmi, iomutex,
                               outLog, hctrs, mopts);
      singleParserPtr->stop();

    }
    if (!mopts->quiet) { std::cerr << "\n\n"; }


    consoleLog->info("Done mapping reads.");
    consoleLog->info("In total saw {:n} reads.", hctrs.numReads);
    consoleLog->info("Final # hits per read = {}", hctrs.totHits / static_cast<float>(hctrs.numReads));
    consoleLog->info("flushing output queue.");
    if (outLog) {
      outLog->flush();
    }
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

bool validateOpts(MappingOpts& mopts, spdlog::logger* log) {
  bool valid{true};

  if (!mopts.selAln) {
    log->warn("\n\n"
              "NOTE: It appears you are running rapmap without the selective-alignment (`--selAln`) option.\n"
              "Selective alignment can generally improve both the sensitivity and specificity of mapping,\n"
              "with only a moderate increase in use of computational resources. \n"
              "Unless there is a specific reason to do this (e.g. testing on clean simulated data),\n"
              "`--selAln` is generally recommended.\n");
  }

  if (mopts.maxMMPExtension < 1) {
    log->error("--maxMMPExtension must be at least 1, but {} was provided.", mopts.maxMMPExtension);
    valid = false;
  }
  if(mopts.selAln) {
    if (mopts.consensusSlack < 0 or mopts.consensusSlack > 1) {
      log->error("--consensusSlack must be between 0.0 and 1.0, you passed {}.", mopts.consensusSlack);
      valid = false;
    }
    if (mopts.minScoreFraction < 0 or mopts.minScoreFraction > 1) {
      log->error("--minScoreFrac slack must be between 0.0 and 1.0, you passed {}.", mopts.minScoreFraction);
      valid = false;
    }
    if (mopts.matchScore <= 0) {
      log->error("match score must be positive!");
      valid = false;
    }
    if (mopts.mismatchPenalty > 0) {
      log->error("mismatch penalty cannot be positve!");
      valid = false;
    }
    int32_t scoreVal = 2 * (static_cast<int32_t>(mopts.gapOpenPenalty) +
                            static_cast<int32_t>(mopts.gapExtendPenalty)) +
      static_cast<int32_t>(mopts.matchScore);
    int32_t maxAllowed = std::numeric_limits<int8_t>::max();
    if (scoreVal >= maxAllowed) {
      log->error(" [2* (gapOpen + gapExtend) + matchScore] cannot exceed {}!", maxAllowed);
      valid = false;
    }
  }
  return valid;
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
        optWriter.write("do chaining : {}\n", mopts.doChaining); 
        optWriter.write("selective alignment: {}\n", mopts.selAln); 
        optWriter.write("no orphans: {}\n", mopts.noOrphans); 
        optWriter.write("no dovetail: {}\n", mopts.noDovetail); 
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
  TCLAP::SwitchArg nosensitive("", "noSensitive", "Perform a less sensitive quasi-mapping by enabling NIP skipping", false);
  TCLAP::SwitchArg noStrict("", "noStrictCheck", "Don't perform extra checks to try and assure that only equally \"best\" mappings for a read are reported", false);
  TCLAP::SwitchArg fuzzy("f", "fuzzyIntersection", "Find paired-end mapping locations using fuzzy intersection", false);
  TCLAP::SwitchArg chain("c", "chaining", "Score the hits to find the best chain", false);
  TCLAP::SwitchArg compressedOutput("x", "compressed", "Compress the output SAM file using zlib", false);
  TCLAP::SwitchArg quiet("q", "quiet", "Disable all console output apart from warnings and errors", false);
  TCLAP::SwitchArg writeUnmapped("u", "writeUnmapped", "include unmapped reads in the output SAM records", false);

  TCLAP::SwitchArg recoverOrphans("", "recoverOrphans", "perform orphan recovery to try and recover endpoints of orphaned reads", false);
  TCLAP::SwitchArg noDovetail("", "noDovetail", "discard dovetailing mappings", false);
  TCLAP::SwitchArg noOrphans("", "noOrphans", "discard orphaned mappings", false);
  TCLAP::SwitchArg selAln("s", "selAln", "Perform selective alignment to validate mapping hits", false);
  TCLAP::ValueArg<int16_t> gapOpenPen("", "go", "[only with selAln]: gap open penalty", false, 4, "positive integer");
  TCLAP::ValueArg<int16_t> gapExtendPen("", "ge", "[only with selAln]: gap extend penalty", false, 2, "positive integer");
  TCLAP::ValueArg<int16_t> mismatchPen("", "mm", "[only with selAln]: mismatch penalty", false, -4, "negative integer");
  TCLAP::ValueArg<int16_t> matchScore("", "ma", "[only with selAln]: match score", false, 2, "positive integer");
  TCLAP::ValueArg<int32_t> dpBandwidth("", "dpBandwidth", "[only with selAln]: bandwidth to use in extension alignment", false, 15, "positive integer");
  TCLAP::ValueArg<double> minScoreFrac("", "minScoreFrac", "[only with selAln]: minimum score fraction for a valid alignment", false, 0.65, "ratio in (0,1]");
  TCLAP::ValueArg<double> consensusSlack("", "consensusSlack", "[only with selAln]: consensus slack to apply during mapping", false, 0.2, "ratio in (0,1]");
  TCLAP::SwitchArg hardFilter("", "hardFilter", "[only with selAln]: apply hard filter to only return best alignments for each read", false);
  TCLAP::SwitchArg mimicBT2("", "mimicBT2", "[only with selAln]: mimic Bowtie2-like default params but with --no-mixed and --no-discordant", false);
  TCLAP::SwitchArg mimicStrictBT2("", "mimicStrictBT2", "[only with selAln]: mimic strict Bowtie2-like default params (e.g. like those recommended for running RSEM)", false);
  TCLAP::ValueArg<int32_t> maxMMPExtension("", "maxMMPExtension", "[only with selAln or with chaining]: maximum allowable MMP extension", false, 7, "positive integer > 1");

  cmd.add(noout);
  cmd.add(maxNumHits);
  cmd.add(quasiCov);
  cmd.add(nosensitive);
  cmd.add(noStrict);
  cmd.add(fuzzy);
  cmd.add(chain);
  cmd.add(quiet);
  cmd.add(writeUnmapped);
  cmd.add(hardFilter);
  cmd.add(recoverOrphans);
  cmd.add(noDovetail);
  cmd.add(noOrphans);
  cmd.add(dpBandwidth);
  cmd.add(gapExtendPen);
  cmd.add(gapOpenPen);
  cmd.add(mismatchPen);
  cmd.add(matchScore);
  cmd.add(consensusSlack);
  cmd.add(maxMMPExtension);
  cmd.add(minScoreFrac);
  cmd.add(mimicStrictBT2);
  cmd.add(mimicBT2);
  cmd.add(selAln);
  cmd.add(numThreads);
  cmd.add(compressedOutput);
  cmd.add(outname);
  cmd.add(unmatedReads);
  cmd.add(read2);
  cmd.add(read1);
  cmd.add(index);

  auto consoleSink = std::make_shared<spdlog::sinks::ansicolor_stderr_sink_mt>();
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
    mopts.writeUnmapped = writeUnmapped.getValue();
    mopts.numThreads = numThreads.getValue();
    mopts.maxNumHits = maxNumHits.getValue();
    mopts.outname = (outname.isSet()) ? outname.getValue() : "";
    mopts.quasiCov = quasiCov.getValue();
    mopts.noOutput = noout.getValue();
    mopts.sensitive = !nosensitive.getValue();
    mopts.strictCheck = !noStrict.getValue();
    mopts.consistentHits = false;//consistent.getValue();
    mopts.doChaining = chain.getValue();
    mopts.compressedOutput = compressedOutput.getValue();
    mopts.fuzzy = fuzzy.getValue();
    mopts.quiet = quiet.getValue();
    mopts.noOrphans = noOrphans.getValue();
    mopts.noDovetail = noDovetail.getValue();
    mopts.maxMMPExtension = maxMMPExtension.getValue();

    mopts.selAln = selAln.getValue();
    if ( (mimicBT2.getValue() or mimicStrictBT2.getValue()) and !mopts.selAln ) {
      consoleLog->info("--mimicBT2 and --mimicStrictBT2 imply --selAln, turning selective alignment on.");
      mopts.selAln = true;
    }

    // selective alignment implies chaining
    if(mopts.selAln) {
      if (!mopts.doChaining) {
        consoleLog->info("The --selAln option implies --chaining.  Setting --chaining to true.");
        mopts.doChaining = true;
      }
      mopts.gapOpenPenalty = gapOpenPen.getValue();
      mopts.gapExtendPenalty = gapExtendPen.getValue();
      mopts.mismatchPenalty = mismatchPen.getValue();
      mopts.matchScore = matchScore.getValue();
      mopts.dpBandwidth = dpBandwidth.getValue();
      mopts.minScoreFraction = minScoreFrac.getValue();
      mopts.consensusSlack = consensusSlack.getValue();
      mopts.hardFilter = hardFilter.getValue();
      if (mimicBT2.getValue() and mimicStrictBT2.getValue()) {
        consoleLog->info("Cannot set --mimicBT2 and --mimicStrictBT2 simultaneously.  Please choose one");
        consoleLog->flush();
        spdlog::drop_all();
        std::exit(1);
      }
      if (mimicBT2.getValue()) {
        consoleLog->info("Setting parameters implied by --mimicBT2.");
        mopts.ap = selective_alignment::utils::AlignmentPolicy::BT2;
        mopts.noOrphans = true;
        mopts.noDovetail = true;
        mopts.consensusSlack = 0.35;
        mopts.maxNumHits = 1000;
      }
      if (mimicStrictBT2.getValue()) {
        consoleLog->info("Setting parameters implied by --mimicStrictBT2.");
        mopts.ap = selective_alignment::utils::AlignmentPolicy::BT2_STRICT;
        mopts.noOrphans = true;
        mopts.noDovetail = true;
        mopts.consensusSlack = 0.35;
        mopts.maxNumHits = 1000;
        mopts.minScoreFraction = 0.8;
        mopts.matchScore = 1;
        mopts.mismatchPenalty = 0;
        // NOTE: as a limitation of ksw2, we can't have
        // (gapOpenPenalty + gapExtendPenalty) * 2 + matchScore < numeric_limits<int8_t>::max()
        // these parameters below are sufficiently large penalties to
        // prohibit gaps, while not overflowing the above condition
        mopts.gapOpenPenalty = 25;
        mopts.gapExtendPenalty = 25;
      }
    }

    mopts.recoverOrphans = recoverOrphans.getValue();

    if (quasiCov.isSet() and nosensitive.isSet()) {
        consoleLog->info("The --quasiCoverage option is set to {}, but the --noSensitive flag was also set. The former forbids the later. Enabling sensitive mode.", quasiCov.getValue());
        mopts.sensitive = true;
    }

    if (!validateOpts(mopts, consoleLog.get())) {
        consoleLog->error("Failed to validate provided options!");
        consoleLog->flush();
        spdlog::drop_all();
        std::exit(1);
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
          RapMapSAIndex<int64_t, PerfectHashT<uint64_t, rapmap::utils::SAInterval<int64_t>>> rmi;
          rmi.load(indexPrefix);
          success = mapReads(rmi, consoleLog, &mopts);
      } else {
          RapMapSAIndex<int64_t,
                        RegHashT<uint64_t, rapmap::utils::SAInterval<int64_t>,
                                               rapmap::utils::KmerKeyHasher>> rmi;
          rmi.load(indexPrefix);
          success = mapReads(rmi, consoleLog, &mopts);
      }
    } else {
        //std::cerr << "Loading 32-bit suffix array index: \n";
      //SAIdxPtr.reset(new RapMapSAIndex<int32_t>);
      //SAIdxPtr->load(indexPrefix, h.kmerLen());
        if (h.perfectHash()) {
            RapMapSAIndex<int32_t, PerfectHashT<uint64_t, rapmap::utils::SAInterval<int32_t>>> rmi;
            rmi.load(indexPrefix);
            success = mapReads(rmi, consoleLog, &mopts);
        } else {
            RapMapSAIndex<int32_t,
                          RegHashT<uint64_t, rapmap::utils::SAInterval<int32_t>,
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
