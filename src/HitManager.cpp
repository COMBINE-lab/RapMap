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

#include "HitManager.hpp"
#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
#include "chobo/small_vector.hpp"
#include <type_traits>

namespace rapmap {
    namespace hit_manager {

      // fastlog2 and fasterlog2 from : https://raw.githubusercontent.com/romeric/fastapprox/master/fastapprox/src/fastlog.h
      // license information included at the bottom of this file.
      static inline float fastlog2 (float x) {
        union { float f; uint32_t i; } vx = { x };
        union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
        float y = vx.i;
        y *= 1.1920928955078125e-7f;

        return y - 124.22551499f
          - 1.498030302f * mx.f 
          - 1.72587999f / (0.3520887068f + mx.f);
      }

      static inline float fasterlog2 (float x) {
        union { float f; uint32_t i; } vx = { x };
        float y = vx.i;
        y *= 1.1920928955078125e-7f;
        return y - 126.94269504f;
      }

    	// Return hits from processedHits where position constraints
        // match maxDist
        bool collectHitsSimple(std::vector<ProcessedHit>& processedHits,
                uint32_t readLen,
                uint32_t maxDist,
                std::vector<QuasiAlignment>& hits,
                MateStatus mateStatus){
            //bool foundHit{false};
            // One processed hit per transcript
            for (auto& ph : processedHits) {
                auto tid = ph.tid;
                std::sort(ph.tqvec.begin(), ph.tqvec.end(),
                        [](const TxpQueryPos& x, const TxpQueryPos& y) -> bool {
                        return x.txpPosInfo.pos() < y.txpPosInfo.pos();
                        });
                auto& firstHit = ph.tqvec[0];
                bool hitRC = firstHit.queryRC;
                bool txpRC = ph.tqvec[0].txpPosInfo.isRC();
                bool isFwd = (hitRC == txpRC);
                int32_t hitPos = firstHit.txpPosInfo.pos() - firstHit.queryPos;

                // determine forward
                hits.emplace_back(tid, hitPos, isFwd, readLen);
                hits.back().mateStatus = mateStatus;
                hits.back().mateIsFwd = true;
            }

            return true;
        }


        // Return hits from processedHits where position constraints
        // match maxDist
        bool collectHitsSimpleSA(SAHitMap& processedHits,
                        uint32_t readLen,
                        int32_t maxDist,
                        std::vector<QuasiAlignment>& hits,
                        MateStatus mateStatus,
                        rapmap::utils::MappingConfig& mc) {
          //bool findBestChain){

          auto findBestChain = mc.doChaining;
          bool considerMultiPos = mc.considerMultiPos;

          //bool foundHit{false};
          // One processed hit per transcript
          //auto startOffset = hits.size();

          std::vector<double> f;
          std::vector<int32_t> p;
          chobo::small_vector<int32_t> bestChainEndInds;
          for (auto& ph : processedHits) {
            // If this is an *active* position list
            if (ph.second.active) {
              auto tid = ph.first;

              if (findBestChain) {
                // IF WE ARE GOING TO VALIDATE THE MAPPING LATER
                auto& hitVector = ph.second.tqvec;
                // sort the hits within this transcript, first by query position and then reference position.
                std::sort(hitVector.begin(), hitVector.end(),
                          [](const SATxpQueryPos& p1, const SATxpQueryPos& p2) -> bool {
                            auto r1 = p1.pos + p1.len;
                            auto r2 = p2.pos + p2.len;
                            auto q1 = p1.queryPos + p1.len;
                            auto q2 = p2.queryPos + p2.len;

                            return (r1 < r2) ? true :    // p1 ref < p2 ref
                                   ((r2 < r1 ) ? false : // p2 ref < p1 ref
                                    (q1 < q2));         // p1 ref == p2 ref (so sort by query pos)
                          });
                auto minPosIt = hitVector.begin();
                // find the valid chains
                // Use variant of minimap2 scoring (Li 2018) 
                // https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty191/4994778
                auto alpha = [](int32_t qdiff, int32_t rdiff, int32_t ilen) -> double {
                  double score = ilen;
                  double mindiff = (qdiff < rdiff) ? qdiff : rdiff;
                  return (score < mindiff) ? score : mindiff;
                };

                auto beta = [maxDist](int32_t qdiff, int32_t rdiff, double avgseed) -> double {
                  if (qdiff < 0 or (std::max(qdiff, rdiff) > maxDist)) {
                    return std::numeric_limits<double>::infinity();
                  }
                  double l = qdiff - rdiff;
                  int32_t al = std::abs(l);
                  return (l == 0) ? 0.0 : (0.01 * avgseed * al + 0.5 * fastlog2(static_cast<float>(al)));
                };

                double bottomScore = std::numeric_limits<double>::lowest();
                double bestScore = bottomScore;
                int32_t bestChainEnd = -1;
                double avgseed = 31.0;

                bestChainEndInds.clear();
                f.clear();
                p.clear();
                auto lastHitId = static_cast<int32_t>(hitVector.size() - 1);
                for (int32_t i = 0; i < static_cast<int32_t>(hitVector.size()); ++i) {
                  auto& hi = hitVector[i];

                  auto qposi = hi.queryPos + hi.len;
                  auto rposi = hi.pos + hi.len;

                  double baseScore = static_cast<double>(hi.len);
                  p.push_back(i);
                  f.push_back(baseScore);

                  // possible predecessors in the chain
                  int32_t numRounds{2};
                  (void)numRounds;
                  for (int32_t j = i-1; j >= 0; --j) {
                    auto& hj = hitVector[j];

                    auto qposj = hj.queryPos + hj.len;
                    auto rposj = hj.pos + hj.len;

                    auto qdiff = qposi - qposj;
                    auto rdiff = rposi - rposj;

                    auto extensionScore = f[j] + alpha(qdiff, rdiff, hi.len) - beta(qdiff, rdiff, avgseed);
                    bool extendWithJ = (extensionScore > f[i]);
                    p[i] = extendWithJ ? j : p[i];
                    f[i] = extendWithJ ? extensionScore : f[i];
                    // HEURISTIC : if we connected this match to an earlier one
                    // i.e. if we extended the chain.
                    // This implements Heng Li's heuristic ---
                    // "
                    // We note that if anchor i is chained to j, chaining i to a predecessor of j
                    // is likely to yield a lower score.
                    // "
                    // here we take this to the extreme, and stop at the first j to which we chain.
                    // we can add a parameter "h" as in the minimap paper.  But here we expect the
                    // chains of matches in short reads to be short enough that this may not be worth it.
                    if (p[i] < i) {
                      numRounds--;
                      if (numRounds <= 0) { break; }
                    }
                  }
                  if (f[i] > bestScore) {
                    bestScore = f[i];
                    bestChainEnd = i;
                    if (considerMultiPos) {
                      bestChainEndInds.clear();
                      bestChainEndInds.push_back(bestChainEnd);
                    }
                  } else if (considerMultiPos and f[i] == bestScore) {
                    bestChainEndInds.push_back(i);
                  }
                }

                // Since we use bestChainEndInds for backtracking,
                // regardless of if we are considering multiple positions
                // or not, make sure we populate it here if we have not above.
                if (!considerMultiPos) {
                  bestChainEndInds.push_back(bestChainEnd);
                }

                //  ===== Do backtracking ========

                // Multi-chain backtracking
                // =========================
                //
                // If considerMultiPos is true, then bestChainEndInds can have
                // >1 entry (the number of equally-optimal chains).  If
                // considerMultiPos is false, then bestChainEndInds will have
                // only a single (best-scoring) chain end, and the below will
                // behave as single-chain backtracking.
                auto numEnds = f.size();
                size_t numDistinctOpt{0};
                chobo::small_vector<int8_t> seen(numEnds, 0);
                chobo::small_vector<decltype(minPosIt)> startPositions;
                auto lastChainHit = bestChainEnd;
                for (auto bestChainEndInd : bestChainEndInds) {
                  bool validChain{true};
                  if (bestChainEndInd >= 0) {
                    auto lastPtr = p[bestChainEndInd];
                    while (lastPtr < bestChainEndInd) {
                      if (seen[bestChainEndInd]) {
                        validChain = false;
                        break;
                      }
                      // It is always OK to mark this, because even if we end up
                      // not using this chain (i.e. it later encounters a seen mem)
                      // then _all_ chains passing through this mem would encounter the
                      // same seen mem and hence be invalid.
                      seen[bestChainEndInd] = 1;
                      bestChainEndInd = lastPtr;
                      lastPtr = p[bestChainEndInd];
                    }
                    if (seen[bestChainEndInd]) {
                      validChain = false;
                    }
                    if (validChain) {
                      ++numDistinctOpt;
                      startPositions.push_back(minPosIt + lastPtr);
                    }
                  } else {
                    // should not happen
                    std::cerr << "[FATAL] : Cannot find any valid chain for quasi-mapping\n";
                    std::cerr << "num hits = " << hitVector.size() << "\n";
                    std::cerr << "bestChainEnd = " << bestChainEnd << "\n";
                    std::cerr << "bestChainScore = " << bestScore << "\n";
                    std::exit(1);
                  }
                }

                {
                  auto posIt = startPositions.begin();
                  bool hitRC = (*posIt)->queryRC;
                  bool isFwd = !hitRC;
                  int32_t hitPos = (*posIt)->pos - (*posIt)->queryPos;
                  hits.emplace_back(tid, hitPos, isFwd, readLen);
                  auto& currHit = hits.back();
                  currHit.mateIsFwd = true; // address issue #46
                  currHit.setChainScore(bestScore);
                  currHit.mateStatus = mateStatus;
                  currHit.allPositions.push_back(hitPos);
                  if (startPositions.size() > 1) {
                    currHit.hasMultiPos = true;
                    while (++posIt != startPositions.end()) {
                      currHit.allPositions.push_back((*posIt)->pos - (*posIt)->queryPos);
                    }
                    // Sort the equally-best hit positions in this transcript
                    std::sort(currHit.allPositions.begin(), currHit.allPositions.end());
                  } else {
                    currHit.hasMultiPos = false;
                  }
                }

                // See if we have a gapless chain
                if (hitVector.size() > 1 and numDistinctOpt == 1 and lastChainHit == lastHitId) {
                  // For this to be the case we must have that
                  // the length of the chain on the query (L) is *exactly*
                  // the same as the length of the chain on the reference,
                  // AND that L is the read length.
                  auto& lastHit = hitVector[lastHitId];
                  int64_t queryRange = static_cast<int64_t>(lastHit.queryPos + lastHit.len) - minPosIt->queryPos;
                  int64_t refRange = static_cast<int64_t>(lastHit.pos + lastHit.len) - minPosIt->pos;
                  int64_t signedReadLen = static_cast<int64_t>(readLen);
                  if (queryRange == refRange and queryRange == signedReadLen) {
                    rapmap::utils::ChainStatus s = rapmap::utils::ChainStatus::UNGAPPED;
                    switch (mateStatus) {
                    case rapmap::utils::MateStatus::SINGLE_END:
                    case rapmap::utils::MateStatus::PAIRED_END_LEFT:
                      hits.back().chainStatus.setLeft(s);
                      break;
                    case rapmap::utils::MateStatus::PAIRED_END_RIGHT:
                      hits.back().chainStatus.setRight(s);
                      break;
                    default:
                      break;
                    }
                  }
                }
                // END OF VALIDATE MAPPING VALIDATION
              } else {
                auto minPosIt = std::min_element(ph.second.tqvec.begin(),
                                            ph.second.tqvec.end(),
                                            [](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
                                              return a.pos < b.pos;
                                            });
                bool hitRC = minPosIt->queryRC;
                int32_t hitPos = minPosIt->pos - minPosIt->queryPos;
                bool isFwd = !hitRC;
                hits.emplace_back(tid, hitPos, isFwd, readLen);
                auto& currHit = hits.back();
                currHit.mateIsFwd = true; // address issue #46
                currHit.mateStatus = mateStatus;
                currHit.allPositions.push_back(hitPos);
              }
            }
          }
          return true;
        }

        /** from http://en.cppreference.com/w/cpp/algorithm/lower_bound **/
        template <typename ForwardIt>
        ForwardIt binarySearch(
                ForwardIt first,
                ForwardIt last,
                uint32_t value) {
            ForwardIt it;
            typename std::iterator_traits<ForwardIt>::difference_type count, step;
            count = std::distance(first, last);

            while (count > 0) {
                it = first;
                step = count / 2;
                std::advance(it, step);
                if (*it < value) {
                    first = ++it;
                    count -= step + 1;
                }
                else {
                    count = step;
                }
            }
            return first;
        }

        /** from http://en.cppreference.com/w/cpp/algorithm/find **/
        template<class InputIt>
        InputIt linearSearch(InputIt first, InputIt last, uint32_t value) {
            for (; first != last; ++first) {
                if (*first == value) {
                    return first;
                }
            }
            return last;
        }

        /** adapted from https://schani.wordpress.com/2010/04/30/linear-vs-binary-search/ **/
        uint32_t binarySearchFast(const std::vector<uint32_t>& arr, size_t n, uint32_t key) {
            uint32_t min = 0, max = n;
            while (min < max) {
                int middle = (min + max) >> 1;
                min = (key > arr[middle]) ? middle+1 : min;
                max = (key <= arr[middle]) ? middle : max;
            }
            return (arr[min] == key) ? min : std::numeric_limits<uint32_t>::max();
        }

        /** adapted from https://schani.wordpress.com/2010/04/30/linear-vs-binary-search/ **/
        // ASSUMES SENTINEL VALUE (value in array >= key *MUST* exist)
        uint32_t linearSearchUnrolled16(const std::vector<uint32_t>& arr, size_t n, uint32_t key) {
            uint32_t i{0};
                for (;;) {
                    if ( arr[i + 0] >= key) return  i + 0;
                    if ( arr[i + 1] >= key) return  i + 1;
                    if ( arr[i + 2] >= key) return  i + 2;
                    if ( arr[i + 3] >= key) return  i + 3;
                    if ( arr[i + 4] >= key) return  i + 4;
                    if ( arr[i + 5] >= key) return  i + 5;
                    if ( arr[i + 6] >= key) return  i + 6;
                    if ( arr[i + 7] >= key) return  i + 7;
                    if ( arr[i + 8] >= key) return  i + 8;
                    if ( arr[i + 9] >= key) return  i + 9;
                    if ( arr[i + 10] >= key) return i + 10;
                    if ( arr[i + 11] >= key) return i + 11;
                    if ( arr[i + 12] >= key) return i + 12;
                    if ( arr[i + 13] >= key) return i + 13;
                    if ( arr[i + 14] >= key) return i + 14;
                    if ( arr[i + 15] >= key) return i + 15;
                    i += 16;
                }
            }

          template <typename RapMapIndexT>
        void intersectSAIntervalWithOutput2(SAIntervalHit<typename RapMapIndexT::IndexType>& h,
                RapMapIndexT& rmi,
                //fbs::eytzinger_array_bfp<uint32_t, uint32_t, true>& outTxps,
                //std::vector<uint32_t>& outTxps,
                SAProcessedHitVec& processedHits) {
            // Convenient bindings for variables we'll use
            auto& SA = rmi.SA;
            auto& txpIDs = rmi.positionIDs;
            auto& txpStarts = rmi.txpOffsets;

            auto& outStructs = processedHits.hits;
            auto& outTxps = processedHits.txps;

            // Iterator to the beginning and end of the output hits
            auto txpIt = processedHits.txps.begin();
            auto txpEnd = processedHits.txps.end();

            uint32_t arraySize = processedHits.txps.size();

            uint32_t rightTxp;
            uint32_t pos;
            //decltype(processedHits.txps)::iterator searchIt = txpEnd;
            uint32_t searchInd{0};
            for (auto i = h.begin; i < h.end; ++i) {
                rightTxp = txpIDs[SA[i]];
                if (arraySize > 64) {
                    searchInd = binarySearchFast(outTxps, arraySize, rightTxp);
                } else {
                    searchInd = linearSearchUnrolled16(outTxps, arraySize, rightTxp);
                }
                // If we found this transcript (make sure it's not the sentinel) then
                // add it to the list.
                if ( searchInd < arraySize - 1 ) {
                    //auto offset = std::distance(txpIt, searchIt);
                    pos = static_cast<uint32_t>(SA[i]) - txpStarts[rightTxp];
                    outStructs[searchInd].tqvec.emplace_back(pos, h.queryPos, h.queryRC);
                }
                /*
                auto searchIdx = outTxps.search(rightTxp);
                if (searchIdx < arraySize) {
                    pos = static_cast<uint32_t>(SA[i]) - txpStarts[rightTxp];
                    outStructs[searchIdx].tqvec.emplace_back(pos, h.queryPos, h.queryRC);
                }
                */
            }
        }


        template <typename RapMapIndexT>
        void intersectSAIntervalWithOutput(SAIntervalHit<typename RapMapIndexT::IndexType>& h,
                                           RapMapIndexT& rmi,
                                           uint32_t intervalCounter,
                                           int32_t maxSlack,
                                           SAHitMap& outHits) {
            using OffsetT = typename RapMapIndexT::IndexType;
            // Convenient bindings for variables we'll use
            auto& SA = rmi.SA;
            //auto& txpIDs = rmi.positionIDs;
            auto& txpStarts = rmi.txpOffsets;
            bool nonStrictIntersection = maxSlack > 0;

            // Walk through every hit in the new interval 'h'
            for (OffsetT i = h.begin; i != h.end; ++i) {
              auto txpID = rmi.transcriptAtPosition(SA[i]);
              auto txpListIt = outHits.find(txpID);
              // True if this transcript is already in the output set, and false otherwise.
              bool inOutputSet = (txpListIt != outHits.end());
              // The number of intervals in which this transcript has occurred.
              int32_t txpOccCount = inOutputSet ? txpListIt->second.numActive : 0;

              // Difference between the number of hits for this transcript and
              // the total number of intervals examined so far.
              int32_t slack = ((static_cast<int32_t>(intervalCounter)- 1) - txpOccCount);

              // If we found this transcript
              // Add this position to the list
              if (nonStrictIntersection or (slack <= maxSlack)) {
                auto globalPos = SA[i];
                auto localPos = globalPos - txpStarts[txpID];
                // We already have records for this transcript
                if (inOutputSet) {
                  txpListIt->second.numActive += (txpListIt->second.lastActiveInterval == intervalCounter) ? 0 : 1;
                  txpListIt->second.lastActiveInterval = intervalCounter;
                  txpListIt->second.tqvec.emplace_back(localPos, h.queryPos, h.queryRC, h.len);
                } else { // We need to add this transcript
                  // The constructor in emplace-back will set numActive = 1.
                  auto& oh = outHits[txpID];
                  oh.tqvec.emplace_back(localPos, h.queryPos, h.queryRC, h.len);
                  oh.lastActiveInterval = intervalCounter;
                }
              }
            }
          }


        template <typename RapMapIndexT>
        std::vector<ProcessedSAHit> intersectSAHits2(
                std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& inHits,
                RapMapIndexT& rmi
                ) {
            using OffsetT = typename RapMapIndexT::IndexType;

            // Each inHit is a SAIntervalHit structure that contains
            // an SA interval with all hits for a particuar query location
            // on the read.
            //
            // We want to find the transcripts that appear in *every*
            // interavl.  Further, for each transcript, we want to
            // know the positions within this txp.

            // Check this --- we should never call this function
            // with less than 2 hits.
            SAProcessedHitVec outHits;
            if (inHits.size() < 2) {
                std::cerr << "intersectHitsSA() called with < 2 k-mer "
                    " hits; this shouldn't happen\n";
                return outHits.hits;
            }

            auto& SA = rmi.SA;
            auto& txpStarts = rmi.txpOffsets;
            auto& txpIDs = rmi.positionIDs;

            // Start with the smallest interval
            // i.e. interval with the fewest hits.
            SAIntervalHit<OffsetT>* minHit = &inHits[0];
            for (auto& h : inHits) {
                if (h.span() < minHit->span()) {
                    minHit = &h;
                }
            }

            auto& outStructs = outHits.hits;
            auto& outTxps = outHits.txps;
            outStructs.reserve(minHit->span());
            outTxps.reserve(minHit->span());
            std::map<int, uint32_t> posMap;
            // =========
            //{ // Add the info from minHit to outHits
                for (int i = minHit->begin; i < minHit->end; ++i) {
                    auto globalPos = SA[i];
                    auto tid = txpIDs[globalPos];
                    auto txpPos = globalPos - txpStarts[tid];
                    auto posIt = posMap.find(tid);
                    if (posIt == posMap.end()) {
                        posMap[tid] = outStructs.size();
                        outStructs.emplace_back(tid, txpPos, minHit->queryPos, minHit->queryRC);
                    } else {
                        outStructs[posIt->second].tqvec.emplace_back(txpPos, minHit->queryPos, minHit->queryRC);
                    }
                }
                std::sort(outStructs.begin(), outStructs.end(),
                          [] (const ProcessedSAHit& a, const ProcessedSAHit& b) -> bool {
                            return a.tid < b.tid;
                          });
                for (auto it = outStructs.begin(); it != outStructs.end(); ++it) {
                    outTxps.emplace_back(it->tid);
                }
                // Sentinel value for search
                outTxps.emplace_back(std::numeric_limits<uint32_t>::max());
                /*
                fbs::eytzinger_array_bfp<uint32_t, uint32_t, true> searchArray(
                        txpIndices.begin(), txpIndices.size()
                        );
                        */
            //}
            // =========

            // Now intersect everything in inHits (apart from minHits)
            // to get the final set of mapping info.
            for (auto& h : inHits) {
                if (&h != minHit) { // don't intersect minHit with itself
                    intersectSAIntervalWithOutput2(h, rmi, outHits);
                }
            }

            size_t requiredNumHits = inHits.size();
            // Mark as active any transcripts with the required number of hits.
            for (auto it = outStructs.begin(); it != outStructs.end(); ++it) {
                if (it->tqvec.size() >= requiredNumHits) {
                    it->active = true;
                }
            }
            return outStructs;
        }

        template <typename RapMapIndexT>
        SAHitMap intersectSAHits(
                SAIntervalVector<SAIntervalHit<typename RapMapIndexT::IndexType>>& inHits,
                RapMapIndexT& rmi,
                size_t readLen,
                float consensusFraction, // The fraction of intervals a transcript must appear in to still be considered a valid mapping
                bool strictFilter) {
          using OffsetT = typename RapMapIndexT::IndexType;
          // Each inHit is a SAIntervalHit structure that contains
          // an SA interval with all hits for a particuar query location
          // on the read.
          //
          // We want to find the transcripts that appear in *every*
          // interavl.  Further, for each transcript, we want to
          // know the positions within this txp.

          // Check this --- we should never call this function
          // with less than 2 hits.
          SAHitMap outHits;
          if (inHits.size() < 2) {
            std::cerr << "intersectHitsSA() called with < 2 hits "
              " hits; this shouldn't happen\n";
            return outHits;
          }

          int32_t sInHitsSize = static_cast<int32_t>(inHits.size());
          float requiredFrac = sInHitsSize * consensusFraction;
          // The number of SA intervals a target must appear in to generate a valid mapping
          int32_t requiredNumHits = sInHitsSize;
          // The maximum slack (maximum number of SA intervals a target can miss and still generate
          // a valid mapping).
          int32_t maxSlack = 0;

          // If the consensusFraction is less than one, then we will allow some slack.
          // We compute the integer number of required intervals (and hence the maximum allowable slack)
          // in this case below.
          if(consensusFraction < 1.0) {
            // always require at least one hit
            requiredNumHits = std::max(static_cast<int32_t>(1), static_cast<int32_t>(std::floor(requiredFrac)));
            // the maximum slack is simply the number of intervals minus the required num hits
            maxSlack = sInHitsSize - requiredNumHits;
          }

          auto& SA = rmi.SA;
          auto& txpStarts = rmi.txpOffsets;
          //auto& txpIDs = rmi.positionIDs;

          // Start with the smallest interval
          // i.e. interval with the fewest hits.
          SAIntervalHit<OffsetT>* minHit = &inHits[0];
          for (auto& h : inHits) {
            if (h.span() < minHit->span()) {
              minHit = &h;
            }
          }

          //outHits.reserve(minHit->span());
          // =========
          { // Add the info from minHit to outHits
            for (OffsetT i = minHit->begin; i < minHit->end; ++i) {
              auto globalPos = SA[i];
              auto tid = rmi.transcriptAtPosition(globalPos);
              auto txpPos = globalPos - txpStarts[tid];
              auto& oh = outHits[tid];
              oh.tqvec.emplace_back(txpPos, minHit->queryPos, minHit->queryRC, minHit->len);
              oh.lastActiveInterval = 1;
            }
          }
          // =========

          // Now intersect everything in inHits (apart from minHits)
          // to get the final set of mapping info.
          size_t intervalCounter{2};
          for (auto& h : inHits) {
            if (&h != minHit) { // don't intersect minHit with itself
              intersectSAIntervalWithOutput(h, rmi, intervalCounter, maxSlack, outHits);
              ++intervalCounter;
            }
          }

          size_t numActive{0};
          // Mark as active any transcripts with the required number of hits.
          for (auto it = outHits.begin(); it != outHits.end(); ++it) {
            int32_t intervalsForHit = it->second.numActive;
            //maxNumIntervals = (intervalsForHit > maxNumIntervals) ? intervalsForHit : maxNumIntervals;
            bool enoughHits = (intervalsForHit >= requiredNumHits);
            bool setActive = (strictFilter) ?
              (enoughHits and it->second.checkConsistent(readLen, requiredNumHits)) :
              (enoughHits);
            it->second.active = setActive;
            numActive += setActive ? 1 : 0;
          }

          // If we had no valid hits according to the interval criteria, then,
          // unless we are using a strict intersection, allow all hits to be potentially considered.
          if (maxSlack > 0 and numActive == 0) {
            for (auto it = outHits.begin(); it != outHits.end(); ++it) {
              it->second.active = true;
            }
          }

          return outHits;
        }

      template <typename RapMapIndexT>
      void hitsToMappingsSimple(RapMapIndexT& rmi,
                                rapmap::utils::MappingConfig& mc,
                                rapmap::utils::MateStatus mateStatus,
                                HitCollectorInfo<SAIntervalHit<typename RapMapIndexT::IndexType>>& hcinfo,
                                std::vector<rapmap::utils::QuasiAlignment>& hits) {

        using OffsetT = typename RapMapIndexT::IndexType;
        auto& txpStarts = rmi.txpOffsets;
        auto& SA = rmi.SA;
        auto& fwdSAInts = hcinfo.fwdSAInts;
        auto& rcSAInts = hcinfo.rcSAInts;
        auto readLen = hcinfo.readLen;
        auto maxDist = hcinfo.maxDist;

        auto consistentHits = mc.consistentHits;
        //auto doChaining = mc.doChaining;
        bool considerMultiPos = mc.considerMultiPos;
        float consensusFraction = mc.consensusFraction;

        auto fwdHitsStart = hits.size();

        // If the hit we have for a read is contained within a single suffix array interval
        // (i.e. there is only one MMP), then this function is called to collect the relevant
        // hits into `outHits`.
        auto collectFromSingleInterval = [&](decltype(hcinfo.fwdSAInts)& saInts,
                                             bool isFw,
                                             std::vector<rapmap::utils::QuasiAlignment>& outHits) {
            auto& saIntervalHit = saInts.front();
            auto initialSize = outHits.size();

            // First, we put all (txp, position) pairs into their own
            // QuasiAlignment
            for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
              auto globalPos = SA[i];
              auto txpID = rmi.transcriptAtPosition(globalPos);
              // the offset into this transcript
              auto pos = globalPos - txpStarts[txpID];
              int32_t hitPos = pos - saIntervalHit.queryPos;
              outHits.emplace_back(txpID, hitPos, isFw, readLen);
              auto& lastHit = outHits.back();
              lastHit.mateIsFwd = true; // address issue #46

              //lastHit.queryOffset = saIntervalHit.queryPos;
              lastHit.mateStatus = mateStatus;
              lastHit.allPositions.push_back(hitPos);
              lastHit.hasMultiPos = false;
              switch (mateStatus) {
              case rapmap::utils::MateStatus::PAIRED_END_LEFT:
              case rapmap::utils::MateStatus::SINGLE_END:
                lastHit.chainStatus.setLeft( (saIntervalHit.len == readLen) ? rapmap::utils::ChainStatus::PERFECT :
                                             rapmap::utils::ChainStatus::REGULAR );
                break;
              case rapmap::utils::MateStatus::PAIRED_END_RIGHT:
                lastHit.chainStatus.setRight( (saIntervalHit.len == readLen) ? rapmap::utils::ChainStatus::PERFECT :
                                              rapmap::utils::ChainStatus::REGULAR );
                break;
              default:
                break;
              }
              //lastHit.completeMatchType = (saIntervalHit.len == readLen) ? mateStatus : MateStatus::NOTHING;
            }

            // Now sort these QuasiAlignments by transcript ID (then position)
            auto sortStartIt = outHits.begin() + initialSize;
            auto sortEndIt = outHits.end();
            std::sort(sortStartIt, sortEndIt,
                      [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                        if (a.tid == b.tid) {
                          return a.pos < b.pos;
                        } else {
                          return a.tid < b.tid;
                        }
                      });


            // Custom variant of unique algorithm from : https://en.cppreference.com/w/cpp/algorithm/unique
            // that will merge the position lists of adjacent hits from the same transcript.
            auto mergeUnique = [](decltype(outHits.begin()) first, decltype(outHits.end()) last) -> decltype(outHits.end()) {
              if (first == last)
                return last;

              decltype(first) result = first;
              while (++first != last) {
                bool distinct = !(result->tid == first->tid);
                // If we just completed a run of transcript-identical hits
                // then move the new transcript's first hit into the next
                // position.
                if (distinct && ++result != first) {
                  *result = std::move(*first);
                } else if(!distinct) {
                  // If we are looking at a pair of transcript-identical elements
                  // then push the position of the current hit onto the position
                  // list for the element that will remain after uniquifying.
                  result->hasMultiPos = true;
                  result->allPositions.push_back(first->pos);
                }
              }

              // NOTE: positions within transcripts should already be sorted, so we don't need
              // to do it again here.
              return ++result;
            };

            // If we are considering multiple hits, then merge the position lists
            // otherwise, just make the hit lists unique.
            auto newEnd = considerMultiPos ?
               mergeUnique(outHits.begin() + initialSize, outHits.end()) :
               std::unique(
                 outHits.begin() + initialSize, outHits.end(),
                 [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                 return a.tid == b.tid;
               });

            // Get rid of redundant QuasiAlignemnts.
            outHits.resize(std::distance(outHits.begin(), newEnd));
        };


        // If we had > 1 forward hit
        if (fwdSAInts.size() > 1) {
          auto processedHits = rapmap::hit_manager::intersectSAHits(
                                                                    fwdSAInts, rmi, readLen, consensusFraction, consistentHits);
          rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                                   hits, mateStatus, mc);
        } else if (fwdSAInts.size() == 1) { // only 1 hit!
          collectFromSingleInterval(fwdSAInts, true, hits);
        }

        auto fwdHitsEnd = hits.size();
        auto rcHitsStart = fwdHitsEnd;
        // If we had > 1 rc hit
        if (rcSAInts.size() > 1) {
          auto processedHits = rapmap::hit_manager::intersectSAHits(
                                                                    rcSAInts, rmi, readLen, consensusFraction, consistentHits);
          rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                                   hits, mateStatus, mc);
        } else if (rcSAInts.size() == 1) { // only 1 hit!
          collectFromSingleInterval(rcSAInts, false, hits);
        }
        auto rcHitsEnd = hits.size();

        // If we had both forward and RC hits, then merge them
        if ((fwdHitsEnd > fwdHitsStart) and (rcHitsEnd > rcHitsStart)) {
          // Merge the forward and reverse hits
          std::inplace_merge(
                             hits.begin() + fwdHitsStart, hits.begin() + fwdHitsEnd,
                             hits.begin() + rcHitsEnd,
                             [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                               // If we have a tie, then the one with the better chain score
                               // comes *before* the one with the worse chain score.
                               return (a.tid == b.tid) ? a.chainScore() > b.chainScore() : a.tid < b.tid;
                             });

          ////
          auto mergeOrientationUnique = [](decltype(hits.begin()) first, decltype(hits.end()) last) -> decltype(hits.end()) {
            if (first == last)
              return last;

            decltype(first) result = first;
            while (++first != last) {
              bool distinct = !(result->tid == first->tid);
              // If we just completed a run of transcript-identical hits
              // then move the new transcript's first hit into the next
              // position.
              if (distinct && ++result != first) {
                *result = std::move(*first);
              } else if(!distinct) {
                // We copy the hits from the second hit to the opposite strand
                result->oppositeStrandPositions = first->allPositions;
              }
            }

            // NOTE: positions within transcripts should already be sorted, so we don't need
            // to do it again here.
            return ++result;
          };
          auto newEnd = mergeOrientationUnique(hits.begin() + fwdHitsStart, hits.begin() + rcHitsEnd);
          ////


          // And get rid of duplicate transcript IDs
          // TODO: We generally don't want to get rid of duplicate transcripts if we are allowing multiple
          // positions, since we may want to consider both fw and rev on the same transcript.
          //auto newEnd = std::unique(
          //                          hits.begin() + fwdHitsStart, hits.begin() + rcHitsEnd,
          //                          [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
          //                            return a.tid == b.tid;
          //                          });
          hits.resize(std::distance(hits.begin(), newEnd));
        }
      }

        /**
        * Need to explicitly instantiate the versions we use
        */
      using SAIndex32BitDense = RapMapSAIndex<int32_t, RegHashT<uint64_t, rapmap::utils::SAInterval<int32_t>,
									     rapmap::utils::KmerKeyHasher>>;
      using SAIndex64BitDense = RapMapSAIndex<int64_t, RegHashT<uint64_t, rapmap::utils::SAInterval<int64_t>,
									     rapmap::utils::KmerKeyHasher>>;
      using SAIndex32BitPerfect = RapMapSAIndex<int32_t, PerfectHashT<uint64_t, rapmap::utils::SAInterval<int32_t>>>;
      using SAIndex64BitPerfect = RapMapSAIndex<int64_t, PerfectHashT<uint64_t, rapmap::utils::SAInterval<int64_t>>>;

        template
        void intersectSAIntervalWithOutput<SAIndex32BitDense>(SAIntervalHit<int32_t>& h,
                                                              SAIndex32BitDense& rmi, 
                                                              uint32_t intervalCounter, 
                                                              int32_t maxSlack,
                                                              SAHitMap& outHits);

        template
        void intersectSAIntervalWithOutput<SAIndex64BitDense>(SAIntervalHit<int64_t>& h,
                                                              SAIndex64BitDense& rmi, 
                                                              uint32_t intervalCounter, 
                                                              int32_t maxSlack,
                                                              SAHitMap& outHits); 

        template
        SAHitMap intersectSAHits<SAIndex32BitDense>(SAIntervalVector<SAIntervalHit<int32_t>>& inHits,
                                                    SAIndex32BitDense& rmi, size_t readLen, float consensusFraction, bool strictFilter);

        template
        SAHitMap intersectSAHits<SAIndex64BitDense>(SAIntervalVector<SAIntervalHit<int64_t>>& inHits,
                                                    SAIndex64BitDense& rmi, size_t readLen, float consensusFraction, bool strictFilter);

        template
        void intersectSAIntervalWithOutput<SAIndex32BitPerfect>(SAIntervalHit<int32_t>& h,
                                                                SAIndex32BitPerfect& rmi, 
                                                                uint32_t intervalCounter, 
                                                                int32_t maxSlack,
                                                                SAHitMap& outHits);

        template
        void intersectSAIntervalWithOutput<SAIndex64BitPerfect>(SAIntervalHit<int64_t>& h,
                                                                SAIndex64BitPerfect& rmi, 
                                                                uint32_t intervalCounter, 
                                                                int32_t maxSlack,
                                                                SAHitMap& outHits);

        template
        SAHitMap intersectSAHits<SAIndex32BitPerfect>(SAIntervalVector<SAIntervalHit<int32_t>>& inHits,
                                                      SAIndex32BitPerfect& rmi, size_t readLen, float consensusFraction, bool strictFilter);

        template
        SAHitMap intersectSAHits<SAIndex64BitPerfect>(SAIntervalVector<SAIntervalHit<int64_t>>& inHits,
                                                      SAIndex64BitPerfect& rmi, size_t readLen, float consensusFraction, bool strictFilter);
      template
      void hitsToMappingsSimple<SAIndex32BitDense>(SAIndex32BitDense& rmi,
                                                   rapmap::utils::MappingConfig& mc,
                                                   rapmap::utils::MateStatus mateStatus,
                                                   HitCollectorInfo<SAIntervalHit<typename SAIndex32BitDense::IndexType>>& hcinfo, std::vector<rapmap::utils::QuasiAlignment>& hits);
      template
      void hitsToMappingsSimple<SAIndex64BitDense>(SAIndex64BitDense& rmi,
                                                   rapmap::utils::MappingConfig& mc,
                                                   rapmap::utils::MateStatus mateStatus,
                                                   HitCollectorInfo<SAIntervalHit<typename SAIndex64BitDense::IndexType>>& hcinfo, std::vector<rapmap::utils::QuasiAlignment>& hits);

      template
      void hitsToMappingsSimple<SAIndex32BitPerfect>(SAIndex32BitPerfect& rmi,
                                                   rapmap::utils::MappingConfig& mc,
                                                   rapmap::utils::MateStatus mateStatus,
                                                   HitCollectorInfo<SAIntervalHit<typename SAIndex32BitPerfect::IndexType>>& hcinfo, std::vector<rapmap::utils::QuasiAlignment>& hits);

      template
      void hitsToMappingsSimple<SAIndex64BitPerfect>(SAIndex64BitPerfect& rmi,
                                                   rapmap::utils::MappingConfig& mc,
                                                   rapmap::utils::MateStatus mateStatus,
                                                   HitCollectorInfo<SAIntervalHit<typename SAIndex64BitPerfect::IndexType>>& hcinfo, std::vector<rapmap::utils::QuasiAlignment>& hits);
}
}


// LICENSE FOR fastlog code :
/*=====================================================================*
 *                   Copyright (C) 2011 Paul Mineiro                   *
 * All rights reserved.                                                *
 *                                                                     *
 * Redistribution and use in source and binary forms, with             *
 * or without modification, are permitted provided that the            *
 * following conditions are met:                                       *
 *                                                                     *
 *     * Redistributions of source code must retain the                *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer.                                       *
 *                                                                     *
 *     * Redistributions in binary form must reproduce the             *
 *     above copyright notice, this list of conditions and             *
 *     the following disclaimer in the documentation and/or            *
 *     other materials provided with the distribution.                 *
 *                                                                     *
 *     * Neither the name of Paul Mineiro nor the names                *
 *     of other contributors may be used to endorse or promote         *
 *     products derived from this software without specific            *
 *     prior written permission.                                       *
 *                                                                     *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND              *
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,         *
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES               *
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE             *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER               *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,                 *
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES            *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE           *
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                *
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF          *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY              *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE             *
 * POSSIBILITY OF SUCH DAMAGE.                                         *
 *                                                                     *
 * Contact: Paul Mineiro <paul@mineiro.com>                            *
 *=====================================================================*/

