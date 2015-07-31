#include "HitManager.hpp"


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
	    auto startOffset = hits.size();
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
	    std::sort(hits.begin() + startOffset, hits.end(),
		      [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
		      	return a.tid < b.tid;
		      });
            return true;
        }


	// Return hits from processedHits where position constraints
        // match maxDist
        bool collectHitsSimpleSA(std::unordered_map<int, std::vector<SATxpQueryPos>>& processedHits,
                uint32_t readLen,
                uint32_t maxDist,
                std::vector<QuasiAlignment>& hits,
                MateStatus mateStatus){
            bool foundHit{false};
            // One processed hit per transcript
            for (auto& ph : processedHits) {
		// If this is an *active* position list
		if (ph.second.front().active) {
		    auto tid = ph.first;
		    std::sort(ph.second.begin(), ph.second.end(),
				    [](const SATxpQueryPos& x, const SATxpQueryPos& y) -> bool {
				    	return x.pos < y.pos;
				    });
		    auto& firstHit = ph.second.front();
		    bool hitRC = firstHit.queryRC;
		    int32_t hitPos = firstHit.pos - firstHit.queryPos;
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

	void intersectSAIntervalWithOutput(SAIntervalHit& h, 
   				           RapMapSAIndex& rmi, 
				           std::unordered_map<int, std::vector<SATxpQueryPos>>& outHits) {
            // Convenient bindings for variables we'll use
            auto& SA = rmi.SA;
            auto& txpIDs = rmi.positionIDs;
            auto& txpStarts = rmi.txpOffsets;

	    // Walk through every hit in the new interval 'h'
	    for (int i = h.begin; i != h.end; ++i) {
	        auto txpID = txpIDs[SA[i]];
	    	auto txpListIt = outHits.find(txpID);
		// If we found this transcript
		// Add this position to the list
		if (txpListIt != outHits.end()) {
		    auto globalPos = SA[i];
		    auto localPos = globalPos - txpStarts[txpID];
		    txpListIt->second.emplace_back(localPos, h.queryPos, h.queryRC);
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



	std::unordered_map<int, std::vector<SATxpQueryPos>> intersectSAHits(
                std::vector<SAIntervalHit>& inHits,
                RapMapSAIndex& rmi
                ) {

            // Each inHit is a SAIntervalHit structure that contains
            // an SA interval with all hits for a particuar query location 
	    // on the read.
	    // 
	    // We want to find the transcripts that appear in *every*
            // interavl.  Further, for each transcript, we want to
            // know the positions within this txp.

            // Check this --- we should never call this function
            // with less than 2 hits.
            if (inHits.size() < 2) {
                std::cerr << "intersectHitsSA() called with < 2 k-mer "
                    " hits; this shouldn't happen\n";
                return {};
            }

	    auto& SA = rmi.SA;
	    auto& txpStarts = rmi.txpOffsets;
	    auto& txpIDs = rmi.positionIDs;

            // Start with the smallest interval 
            // i.e. interval with the fewest hits. 
            SAIntervalHit* minHit = &inHits[0];
            for (auto& h : inHits) {
                if (h.span() < minHit->span()) {
                    minHit = &h;
                }
            }

	    std::unordered_map<int, std::vector<SATxpQueryPos>> outHits;
            outHits.reserve(minHit->span());
            // =========
            { // Add the info from minHit to outHits
		for (int i = minHit->begin; i < minHit->end; ++i) {
		    auto globalPos = SA[i];
		    auto tid = txpIDs[globalPos];
		    auto txpPos = globalPos - txpStarts[tid];
		    outHits[tid].emplace_back(txpPos, minHit->queryPos, minHit->queryRC);
		}
            }
            // =========

            // Now intersect everything in inHits (apart from minHits)
            // to get the final set of mapping info.
            for (auto& h : inHits) {
                if (&h != minHit) { // don't intersect minHit with itself
                    intersectSAIntervalWithOutput(h, rmi, outHits);
                }
            }

            size_t requiredNumHits = inHits.size();
	    // Mark as active any transcripts with the required number of hits. 
	    for (auto it = outHits.begin(); it != outHits.end(); ++it) {
	        if (it->second.size() >= requiredNumHits) {
		    it->second.front().active = true;
		}
	    }
            return outHits;
        }

    }
}
