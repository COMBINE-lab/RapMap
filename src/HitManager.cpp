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
#include <type_traits>

namespace rapmap {
    namespace hit_manager {
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
                        bool findBestChain){

          //bool foundHit{false};
          // One processed hit per transcript
          //auto startOffset = hits.size();

          std::vector<double> f;
          std::vector<int32_t> p;
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
                  return (l == 0) ? 0.0 : (0.01 * avgseed * al + 0.5 * std::log2(al));
                };

                double bottomScore = std::numeric_limits<double>::lowest();
                double bestScore = bottomScore;
                int32_t bestChainEnd = -1;
                double avgseed = 31.0;
                f.clear();
                p.clear();
                auto lastHitId = static_cast<int32_t>(hitVector.size() - 1);
                for (int32_t i = 0; i < static_cast<int32_t>(hitVector.size()); ++i) {
                  auto& hi = hitVector[i];

                  auto qposi = hi.queryPos;
                  auto rposi = hi.pos;

                  double baseScore = static_cast<double>(hi.len);
                  p.push_back(i);
                  f.push_back(baseScore);

                  // possible predecessors in the chain
                  int32_t numRounds{1};
                  (void)numRounds;
                  for (int32_t j = i-1; j >= 0; --j) {
                    auto& hj = hitVector[j];

                    auto qposj = hj.queryPos;
                    auto rposj = hj.pos;

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
                    if (p[i] < i) { break; }
                  }
                  if (f[i] > bestScore) {
                    bestScore = f[i];
                    bestChainEnd = i;
                  }
                }

                // Do backtracking
                auto lastChainHit = bestChainEnd;
                if (bestChainEnd >= 0) {
                  auto lastPtr = p[bestChainEnd];
                  while (lastPtr < bestChainEnd) {
                    lastPtr = bestChainEnd;
                    bestChainEnd = p[bestChainEnd];
                  }
                  minPosIt += lastPtr;
                } else {
                  // should not happen
                  std::cerr << "[FATAL] : Cannot find any valid chain for quasi-mapping\n";
                  std::cerr << "num hits = " << hitVector.size() << "\n";
                  std::cerr << "bestChainEnd = " << bestChainEnd << "\n";
                  std::cerr << "bestChainScore = " << bestScore << "\n";
                  std::exit(1);
                }

                bool hitRC = minPosIt->queryRC;
                int32_t hitPos = minPosIt->pos - minPosIt->queryPos;
                bool isFwd = !hitRC;
                hits.emplace_back(tid, hitPos, isFwd, readLen);
                hits.back().mateStatus = mateStatus;

                // See if we have a gapless chain
                if (hitVector.size() > 1 and lastChainHit == lastHitId) {
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
                hits.back().mateStatus = mateStatus;
              }
            }
          }
          // if SAHitMap is sorted, no need to sort here
          /*
            std::sort(hits.begin() + startOffset, hits.end(),
            [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
            return a.tid < b.tid;
            });
          */
          return true;
        }


        // Return hits from processedHits where position constraints
        // match maxDist
        bool collectHitsSimpleSA2(std::vector<ProcessedSAHit>& processedHits,
                        uint32_t readLen,
                        uint32_t maxDist,
                        std::vector<QuasiAlignment>& hits,
                        MateStatus mateStatus){
                //bool foundHit{false};

                // One processed hit per transcript
                for (auto& ph : processedHits) {
                        // If this is an *active* position list
                        if (ph.active) {
                                auto tid = ph.tid;
                                auto minPosIt =
                                    std::min_element(ph.tqvec.begin(),
                                                     ph.tqvec.end(),
                                                     [](const SATxpQueryPos& a, const SATxpQueryPos& b) -> bool {
                                                        return a.pos < b.pos;
                                                        });

                                bool hitRC = minPosIt->queryRC;
                                int32_t hitPos = minPosIt->pos - minPosIt->queryPos;
                                bool isFwd = !hitRC;
                                hits.emplace_back(tid, hitPos, isFwd, readLen);
                                hits.back().mateStatus = mateStatus;
                        }
                }
                return true;
        }




        // Intersects the hit h2 with outHits.
        // This will modify outHits so that the tqvec field of the
        // entries in outHits that are labeled by the transcripts in
        // which h2 appears will have an iterator to the beginning of
        // the position list for h2.
        /*void intersectWithOutput(HitInfo& h2, RapMapIndex& rmi,
                std::vector<ProcessedHit>& outHits) {

            // Convenient bindings for variables we'll use
            auto& eqClasses = rmi.eqClassList;
            auto& eqClassLabels = rmi.eqLabelList;
            auto& posList = rmi.posList;

            // Iterator to the beginning and end of the output hits
            auto outHitIt = outHits.begin();
            auto outHitEnd = outHits.end();

            // Equiv. class for h2
            auto& eqClassRight = eqClasses[h2.kinfo->eqId];

            // Iterator into, length of and end of the positon list for h2
            auto rightPosIt = posList.begin() + h2.kinfo->offset;
            auto rightPosLen = h2.kinfo->count;
            // auto rightPosEnd = rightPosIt + rightPosLen;
            // Iterator into, length of and end of the transcript list for h2
            auto rightTxpIt = eqClassLabels.begin() + eqClassRight.txpListStart;
            auto rightTxpListLen = eqClassRight.txpListLen;
            auto rightTxpEnd = rightTxpIt + rightTxpListLen;

            auto rightQueryPos = h2.queryPos;
            auto rightQueryRC = h2.queryRC;
            PositionListHelper rightPosHelper(rightPosIt, posList.end());

            uint32_t leftTxp, rightTxp;
            while (outHitIt != outHitEnd and rightTxpIt != rightTxpEnd) {
                // Get the current transcript ID for the left and right eq class
                leftTxp = outHitIt->tid;
                rightTxp = *rightTxpIt;
                // If we need to advance the left txp, do it
                if (leftTxp < rightTxp) {
                    // Advance to the next transcript in the
                    // equivalence class label
                    ++outHitIt;
                } else {
                    // If the transcripts are equal (i.e. leftTxp >= rightTxp and !(rightTxp < leftTxp))
                    // Then see if there are any hits here.
                    if (!(rightTxp < leftTxp)) {
                        // Add the position list iterator and query pos for the
                        // hit from h2 to the back of outHits' tqvec.
                        outHitIt->tqvec.emplace_back(rightPosHelper, rightQueryPos, rightQueryRC);
                        ++outHitIt;
                    }
                    // advance the hit we're intersecting to the next transcript
                    rightPosHelper.advanceToNextTranscript();
                    // Advance the right transcript id regardless of whether
                    // we found a hit or not.
                    ++rightTxpIt;
                }
            }

        }
        */

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


        /*
        void intersectSAIntervalWithOutput3(SAIntervalHit& h,
                RapMapSAIndex& rmi,
                SAProcessedHitVec& outHits) {
            // Convenient bindings for variables we'll use
            auto& SA = rmi.SA;
            auto& txpIDs = rmi.positionIDs;
            auto& txpStarts = rmi.txpOffsets;

            // Iterator to the beginning and end of the output hits
            auto outHitIt = outHits.begin();
            auto outHitEnd = outHits.end();

            // Make a vector of iterators into the right interval
            std::vector<int*> rightHitIterators;
            rightHitIterators.reserve(h.span());
            for (auto i = h.begin; i < h.end; ++i) {
                rightHitIterators.emplace_back(&SA[i]);
            }
            // Sort the iterators by their transcript ID
            std::sort(rightHitIterators.begin(), rightHitIterators.end(),
                    [&txpIDs](const int* a, const int* b) -> bool {
                    return txpIDs[*a] < txpIDs[*b];
                    });
            auto rightIntHit = rightHitIterators.begin();
            auto rightIntHitEnd = rightHitIterators.end();

            uint32_t leftTxp, rightTxp;
            uint32_t pos;
            while (outHitIt != outHitEnd and rightIntHit != rightIntHitEnd) {
                // Get the current transcript ID for the left and right eq class
                leftTxp = outHitIt->tid;
                rightTxp = txpIDs[(*(*rightIntHit))];
                // If we need to advance the left txp, do it
                if (leftTxp < rightTxp) {
                    // Advance to the next transcript in the
                    // equivalence class label
                    ++outHitIt;
                } else {
                    // If the transcripts are equal (i.e. leftTxp >= rightTxp and !(rightTxp < leftTxp))
                    // Then see if there are any hits here.
                    if (!(rightTxp < leftTxp)) {
                        // Add the position list iterator and query pos for the
                        // hit from h2 to the back of outHits' tqvec.
                        pos = static_cast<uint32_t>(*(*rightIntHit)) - txpStarts[rightTxp];
                        outHitIt->tqvec.emplace_back(pos, h.queryPos, h.queryRC);
                        //++outHitIt;
                    }
                    ++rightIntHit;
                }
            }
        }
        */



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

            // Walk through every hit in the new interval 'h'
            for (OffsetT i = h.begin; i != h.end; ++i) {
              auto txpID = rmi.transcriptAtPosition(SA[i]);
              auto txpListIt = outHits.find(txpID);
              // True if this transcript is already in the output set, and false otherwise.
              bool inOutputSet = (txpListIt != outHits.end());
              // The number of intervals in which this transcript has occurred.
              int32_t txpOccCount = inOutputSet ? 0 : txpListIt->second.numActive;

              // Difference between the number of hits for this transcript and
              // the total number of intervals examined so far.
              int32_t slack = ((static_cast<int32_t>(intervalCounter)- 1) - txpOccCount);

              // If we found this transcript
              // Add this position to the list
              if (slack <= maxSlack) {
                auto globalPos = SA[i];
                auto localPos = globalPos - txpStarts[txpID];
                // We already have records for this transcript
                if (inOutputSet) {
                  txpListIt->second.numActive += (txpListIt->second.numActive == intervalCounter - 1) ? 1 : 0;
                  txpListIt->second.tqvec.emplace_back(localPos, h.queryPos, h.queryRC, h.len);
                } else { // We need to add this transcript
                  // The constructor in emplace-back 
                  auto& oh = outHits[txpID];
                  oh.tqvec.emplace_back(localPos, h.queryPos, h.queryRC, h.len);
                }
              }
            }
          }


          /*
        std::vector<ProcessedHit> intersectHits(
                std::vector<HitInfo>& inHits,
                RapMapIndex& rmi
                ) {
            // Each inHit is a HitInfo structure that contains
            // an iterator to the KmerInfo for this k-mer, the k-mer ID,
            // and the query position where this k-mer appeared.
            // We want to find the transcripts that appear in *every*
            // hit.  Further, for each transcript, we want to
            // know the k-mers that appear in this txp.

            // Check this --- we should never call this function
            // with less than 2 hits.
            if (inHits.size() < 2) {
                std::cerr << "intersectHits() called with < 2 k-mer "
                    " hits; this shouldn't happen\n";
                return {};
            }

            auto& eqClasses = rmi.eqClassList;
            auto& eqClassLabels = rmi.eqLabelList;
            auto& posList = rmi.posList;

            // The HitInfo with the smallest equivalence class
            // i.e. label with the fewest transcripts.
            HitInfo* minHit = &inHits[0];
            for (auto& h : inHits) {
                if (h.kinfo->count < minHit->kinfo->count) {
                    minHit = &h;
                }
            }

            std::vector<ProcessedHit> outHits;
            outHits.reserve(minHit->kinfo->count);
            // =========
            { // Add the info from minHit to outHits
                // Equiv. class for minHit
                auto& eqClass = eqClasses[minHit->kinfo->eqId];
                // Iterator into, length of and end of the positon list
                auto posIt = posList.begin() + minHit->kinfo->offset;
                auto posLen = minHit->kinfo->count;
                // auto posEnd = posIt + posLen;
                // Iterator into, length of and end of the transcript list
                auto txpIt = eqClassLabels.begin() + eqClass.txpListStart;
                auto txpListLen = eqClass.txpListLen;
                auto txpEnd = txpIt + txpListLen;
                PositionListHelper posHelper(posIt, posList.end());

                while (txpIt != txpEnd) {
                    auto tid = *txpIt;
                    outHits.emplace_back(tid, posHelper, minHit->queryPos, minHit->queryRC);
                    posHelper.advanceToNextTranscript();
                    ++txpIt;
                }
            }
            // =========

            // Now intersect everything in inHits (apart from minHits)
            // to get the final set of mapping info.
            for (auto& h : inHits) {
                if (&h != minHit) { // don't intersect minHit with itself
                    intersectWithOutput(h, rmi, outHits);
                }
            }

            size_t requiredNumHits = inHits.size();
            // do we need stable_partition? --- don't think so.
            auto newEnd = std::stable_partition(outHits.begin(), outHits.end(),
                    [requiredNumHits] (const ProcessedHit& ph) -> bool {
                    // should never really be greater.
                    return (ph.tqvec.size() >= requiredNumHits);
                    });
d
            // return only the valid hits
            outHits.resize(std::distance(outHits.begin(), newEnd));
            return outHits;
        }
        */

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
                std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& inHits,
                RapMapIndexT& rmi,
                size_t readLen,
                int32_t maxSlack, // The number of intervals a transcript is allowed to miss and still be considered a valid mapping
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

          size_t requiredNumHits = inHits.size() - maxSlack;
          // Mark as active any transcripts with the required number of hits.
          for (auto it = outHits.begin(); it != outHits.end(); ++it) {
            bool enoughHits = (it->second.numActive >= requiredNumHits);
            it->second.active = (strictFilter) ? 
              (enoughHits and it->second.checkConsistent(readLen, requiredNumHits)) :
              (enoughHits);
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
        auto& rankDict = rmi.rankDict;
        auto& txpStarts = rmi.txpOffsets;
        auto& SA = rmi.SA;
        auto& fwdSAInts = hcinfo.fwdSAInts;
        auto& rcSAInts = hcinfo.rcSAInts;
        auto readLen = hcinfo.readLen;
        auto maxDist = hcinfo.maxDist;

        auto consistentHits = mc.consistentHits;
        auto doChaining = mc.doChaining;

        auto fwdHitsStart = hits.size();
        int32_t maxSlack = (doChaining) ? 1 : 0;

        // If we had > 1 forward hit
        if (fwdSAInts.size() > 1) {
          auto processedHits = rapmap::hit_manager::intersectSAHits(
                                                                    fwdSAInts, rmi, readLen, maxSlack, consistentHits);
          rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                                   hits, mateStatus, doChaining);
        } else if (fwdSAInts.size() == 1) { // only 1 hit!
          auto& saIntervalHit = fwdSAInts.front();
          auto initialSize = hits.size();
          for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
            auto globalPos = SA[i];
            auto txpID = rmi.transcriptAtPosition(globalPos);
            // the offset into this transcript
            auto pos = globalPos - txpStarts[txpID];
            int32_t hitPos = pos - saIntervalHit.queryPos;
            hits.emplace_back(txpID, hitPos, true, readLen);
            auto& lastHit = hits.back();
            lastHit.mateStatus = mateStatus;
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
          // Now sort by transcript ID (then position) and eliminate
          // duplicates
          auto sortStartIt = hits.begin() + initialSize;
          auto sortEndIt = hits.end();
          std::sort(sortStartIt, sortEndIt,
                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                      if (a.tid == b.tid) {
                        return a.pos < b.pos;
                      } else {
                        return a.tid < b.tid;
                      }
                    });
          auto newEnd = std::unique(
                                    hits.begin() + initialSize, hits.end(),
                                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                      return a.tid == b.tid;
                                    });
          hits.resize(std::distance(hits.begin(), newEnd));
        }
        auto fwdHitsEnd = hits.size();

        auto rcHitsStart = fwdHitsEnd;
        // If we had > 1 rc hit
        if (rcSAInts.size() > 1) {
          auto processedHits = rapmap::hit_manager::intersectSAHits(
                                                                    rcSAInts, rmi, readLen, maxSlack, consistentHits);
          rapmap::hit_manager::collectHitsSimpleSA(processedHits, readLen, maxDist,
                                                   hits, mateStatus, doChaining);
        } else if (rcSAInts.size() == 1) { // only 1 hit!
          auto& saIntervalHit = rcSAInts.front();
          auto initialSize = hits.size();
          for (OffsetT i = saIntervalHit.begin; i != saIntervalHit.end; ++i) {
            auto globalPos = SA[i];
            auto txpID = rmi.transcriptAtPosition(globalPos);
            // the offset into this transcript
            auto pos = globalPos - txpStarts[txpID];
            int32_t hitPos = pos - saIntervalHit.queryPos;
            hits.emplace_back(txpID, hitPos, false, readLen);
            auto& lastHit = hits.back();
            lastHit.mateStatus = mateStatus;
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
          // Now sort by transcript ID (then position) and eliminate
          // duplicates
          auto sortStartIt = hits.begin() + rcHitsStart;
          auto sortEndIt = hits.end();
          std::sort(sortStartIt, sortEndIt,
                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                      if (a.tid == b.tid) {
                        return a.pos < b.pos;
                      } else {
                        return a.tid < b.tid;
                      }
                    });
          auto newEnd = std::unique(
                                    sortStartIt, sortEndIt,
                                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                      return a.tid == b.tid;
                                    });
          hits.resize(std::distance(hits.begin(), newEnd));
        }
        auto rcHitsEnd = hits.size();

        // If we had both forward and RC hits, then merge them
        if ((fwdHitsEnd > fwdHitsStart) and (rcHitsEnd > rcHitsStart)) {
          // Merge the forward and reverse hits
          std::inplace_merge(
                             hits.begin() + fwdHitsStart, hits.begin() + fwdHitsEnd,
                             hits.begin() + rcHitsEnd,
                             [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                               return a.tid < b.tid;
                             });
          // And get rid of duplicate transcript IDs
          auto newEnd = std::unique(
                                    hits.begin() + fwdHitsStart, hits.begin() + rcHitsEnd,
                                    [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                      return a.tid == b.tid;
                                    });
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
        SAHitMap intersectSAHits<SAIndex32BitDense>(std::vector<SAIntervalHit<int32_t>>& inHits,
                                                    SAIndex32BitDense& rmi, size_t readLen, int32_t maxSlack, bool strictFilter);

        template
        SAHitMap intersectSAHits<SAIndex64BitDense>(std::vector<SAIntervalHit<int64_t>>& inHits,
                                                    SAIndex64BitDense& rmi, size_t readLen, int32_t maxSlack, bool strictFilter);

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
        SAHitMap intersectSAHits<SAIndex32BitPerfect>(std::vector<SAIntervalHit<int32_t>>& inHits,
                                                      SAIndex32BitPerfect& rmi, size_t readLen, int32_t maxSlack, bool strictFilter);

        template
        SAHitMap intersectSAHits<SAIndex64BitPerfect>(std::vector<SAIntervalHit<int64_t>>& inHits,
                                                      SAIndex64BitPerfect& rmi, size_t readLen, int32_t maxSlack, bool strictFilter);
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
