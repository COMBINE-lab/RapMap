#include "HitManager.hpp"
#include "MPHMap.hpp"
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
            bool foundHit{false};
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
                        uint32_t maxDist,
                        std::vector<QuasiAlignment>& hits,
                        MateStatus mateStatus){
                bool foundHit{false};
                // One processed hit per transcript
	            auto startOffset = hits.size();
                for (auto& ph : processedHits) {
                        // If this is an *active* position list
                        if (ph.second.active) {
                                auto tid = ph.first;
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
                bool foundHit{false};

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
        void intersectWithOutput(HitInfo& h2, RapMapIndex& rmi,
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
            auto rightPosEnd = rightPosIt + rightPosLen;
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
          SAHitMap& outHits) {
            using OffsetT = typename RapMapIndexT::IndexType;
            // Convenient bindings for variables we'll use
            auto& SA = rmi.SA;
            //auto& txpIDs = rmi.positionIDs;
            auto& rankDict = rmi.rankDict;
            auto& txpStarts = rmi.txpOffsets;

            // Walk through every hit in the new interval 'h'
            for (OffsetT i = h.begin; i != h.end; ++i) {
              //auto txpID = txpIDs[SA[i]];
              // auto txpID = rankDict.Rank(SA[i], 1);
              auto txpID = rmi.transcriptAtPosition(SA[i]);
              auto txpListIt = outHits.find(txpID);
              // If we found this transcript
              // Add this position to the list
              if (txpListIt != outHits.end()) {
                txpListIt->second.numActive += (txpListIt->second.numActive == intervalCounter - 1) ? 1 : 0;
                if (txpListIt->second.numActive == intervalCounter) {
                  auto globalPos = SA[i];
                  auto localPos = globalPos - txpStarts[txpID];
                  txpListIt->second.tqvec.emplace_back(localPos, h.queryPos, h.queryRC);
                }
              }
            }
          }



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
                auto posEnd = posIt + posLen;
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
            /*
               bool didDrop = false;
               for (auto it = newEnd; it != outHits.end(); ++it) {
               std::cerr << "Dropped hit for txp " << it->tid << "\n";
               didDrop = true;
               }
               if (didDrop) {
               auto& eqClass = eqClasses[inHits[0].kinfo->eqId];
               auto txpIt = eqClassLabels.begin() + eqClass.txpListStart;
               auto txpListLen = eqClass.txpListLen;
               auto txpEnd = txpIt + txpListLen;
               std::cerr << "hits1: {";
               while (txpIt != txpEnd) {
               std::cerr << *txpIt << ", ";
               ++txpIt;
               }
               std::cerr << "}\n";
               auto& eqClass2 = eqClasses[inHits[1].kinfo->eqId];
               txpIt = eqClassLabels.begin() + eqClass2.txpListStart;
               txpListLen = eqClass2.txpListLen;
               txpEnd = txpIt + txpListLen;
               std::cerr << "hits2: {";
               while (txpIt != txpEnd) {
               std::cerr << *txpIt << ", ";
               ++txpIt;
               }
               std::cerr << "}\n";
               }
               */
            // return only the valid hits
            outHits.resize(std::distance(outHits.begin(), newEnd));
            return outHits;
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
            SAHitMap outHits;
            if (inHits.size() < 2) {
                std::cerr << "intersectHitsSA() called with < 2 k-mer "
                    " hits; this shouldn't happen\n";
                return outHits;
            }

            auto& SA = rmi.SA;
            auto& txpStarts = rmi.txpOffsets;
            //auto& txpIDs = rmi.positionIDs;
	    auto& rankDict = rmi.rankDict;

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
                for (int i = minHit->begin; i < minHit->end; ++i) {
                    auto globalPos = SA[i];
                    //auto tid = txpIDs[globalPos];
		    auto tid = rmi.transcriptAtPosition(globalPos);
                    auto txpPos = globalPos - txpStarts[tid];
                    outHits[tid].tqvec.emplace_back(txpPos, minHit->queryPos, minHit->queryRC);
                }
            }
            // =========

            // Now intersect everything in inHits (apart from minHits)
            // to get the final set of mapping info.
	    size_t intervalCounter{2};
            for (auto& h : inHits) {
                if (&h != minHit) { // don't intersect minHit with itself
                    intersectSAIntervalWithOutput(h, rmi, intervalCounter, outHits);
		    ++intervalCounter;
                }
            }

            size_t requiredNumHits = inHits.size();
            // Mark as active any transcripts with the required number of hits.
            for (auto it = outHits.begin(); it != outHits.end(); ++it) {
                if (it->second.numActive >= requiredNumHits) {
                    it->second.active = true;
                }
            }
            return outHits;
        }


        /**
        * Need to explicitly instantiate the versions we use
        */
      using SAIndex32BitDense = RapMapSAIndex<int32_t,google::dense_hash_map<uint64_t, rapmap::utils::SAInterval<int32_t>,
									     rapmap::utils::KmerKeyHasher>>;
      using SAIndex64BitDense = RapMapSAIndex<int64_t,google::dense_hash_map<uint64_t, rapmap::utils::SAInterval<int64_t>,
									     rapmap::utils::KmerKeyHasher>>;
      using SAIndex32BitPerfect = RapMapSAIndex<int32_t, MPHMap<uint64_t, std::pair<uint64_t, rapmap::utils::SAInterval<int32_t>>>>;
      using SAIndex64BitPerfect = RapMapSAIndex<int64_t, MPHMap<uint64_t, std::pair<uint64_t, rapmap::utils::SAInterval<int64_t>>>>;

        template
        void intersectSAIntervalWithOutput<SAIndex32BitDense>(SAIntervalHit<int32_t>& h,
          SAIndex32BitDense& rmi, uint32_t intervalCounter, SAHitMap& outHits);

        template
        void intersectSAIntervalWithOutput<SAIndex64BitDense>(SAIntervalHit<int64_t>& h,
          SAIndex64BitDense& rmi, uint32_t intervalCounter, SAHitMap& outHits);

        template
        SAHitMap intersectSAHits<SAIndex32BitDense>(std::vector<SAIntervalHit<int32_t>>& inHits,
          SAIndex32BitDense& rmi);

        template
        SAHitMap intersectSAHits<SAIndex64BitDense>(std::vector<SAIntervalHit<int64_t>>& inHits,
          SAIndex64BitDense& rmi);

        template
        void intersectSAIntervalWithOutput<SAIndex32BitPerfect>(SAIntervalHit<int32_t>& h,
          SAIndex32BitPerfect& rmi, uint32_t intervalCounter, SAHitMap& outHits);

        template
        void intersectSAIntervalWithOutput<SAIndex64BitPerfect>(SAIntervalHit<int64_t>& h,
          SAIndex64BitPerfect& rmi, uint32_t intervalCounter, SAHitMap& outHits);

        template
        SAHitMap intersectSAHits<SAIndex32BitPerfect>(std::vector<SAIntervalHit<int32_t>>& inHits,
          SAIndex32BitPerfect& rmi);

        template
        SAHitMap intersectSAHits<SAIndex64BitPerfect>(std::vector<SAIntervalHit<int64_t>>& inHits,
          SAIndex64BitPerfect& rmi);


    }
}
