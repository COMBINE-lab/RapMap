#ifndef __HIT_MANAGER_HPP__
#define __HIT_MANAGER_HPP__

#include "RapMapUtils.hpp"
#include "RapMapIndex.hpp"

#include <vector>
#include <algorithm>

namespace rapmap {
    namespace hit_manager {
        using HitInfo = rapmap::utils::HitInfo;
        using ProcessedHit = rapmap::utils::ProcessedHit;
        using MateStatus = rapmap::utils::MateStatus;
        using PositionListHelper = rapmap::utils::PositionListHelper;
        using QuasiAlignment = rapmap::utils::QuasiAlignment;
        using TxpQueryPos = rapmap::utils::TxpQueryPos;


        // Return hits from processedHits where position constraints
        // match maxDist
        bool collectHitsSimple(std::vector<ProcessedHit>& processedHits,
                uint32_t readLen,
                uint32_t maxDist,
                std::vector<QuasiAlignment>& hits,
                MateStatus mateStatus);

        // Intersects the hit h2 with outHits.
        // This will modify outHits so that the tqvec field of the
        // entries in outHits that are labeled by the transcripts in
        // which h2 appears will have an iterator to the beginning of
        // the position list for h2.
        void intersectWithOutput(HitInfo& h2, RapMapIndex& rmi,
                std::vector<ProcessedHit>& outHits); 

        std::vector<ProcessedHit> intersectHits(
                std::vector<HitInfo>& inHits,
                RapMapIndex& rmi); 
    }
}


#endif // __HIT_MANAGER_HPP__
