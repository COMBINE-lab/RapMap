#ifndef __HIT_MANAGER_HPP__
#define __HIT_MANAGER_HPP__

#include "RapMapUtils.hpp"
#include "RapMapIndex.hpp"
#include "RapMapSAIndex.hpp"

#include "eytzinger_array.h"

#include <tuple>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>

namespace rapmap {
    namespace hit_manager {
        using HitInfo = rapmap::utils::HitInfo;
        using ProcessedHit = rapmap::utils::ProcessedHit;
        using MateStatus = rapmap::utils::MateStatus;
        using PositionListHelper = rapmap::utils::PositionListHelper;
        using QuasiAlignment = rapmap::utils::QuasiAlignment;
        using TxpQueryPos = rapmap::utils::TxpQueryPos;
        using SATxpQueryPos = rapmap::utils::SATxpQueryPos;
        using SAIntervalHit = rapmap::utils::SAIntervalHit;
        using SAHitMap = std::map<int, std::vector<SATxpQueryPos>>;
        using ProcessedSAHit = rapmap::utils::ProcessedSAHit;

        class SAProcessedHitVec {
            public:
                std::vector<ProcessedSAHit> hits;
                std::vector<uint32_t> txps;
        };
        /*
        using SAProcessedHitVec = std::tuple<std::vector<ProcessedSAHit>, std::vector<uint32_t>>;
        */

        // Return hits from processedHits where position constraints
        // match maxDist
        bool collectHitsSimple(std::vector<ProcessedHit>& processedHits,
                uint32_t readLen,
                uint32_t maxDist,
                std::vector<QuasiAlignment>& hits,
                MateStatus mateStatus);

        // Return hits from processedHits where position constraints
        // match maxDist
        bool collectHitsSimpleSA(SAHitMap& processedHits,
                uint32_t readLen,
                uint32_t maxDist,
                std::vector<QuasiAlignment>& hits,
                MateStatus mateStatus);

        // Return hits from processedHits where position constraints
        // match maxDist
        bool collectHitsSimpleSA2(std::vector<ProcessedSAHit>& processedHits,
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

        void intersectSAIntervalWithOutput(SAIntervalHit& h,
                RapMapSAIndex& rmi,
                SAHitMap& outHits);

        void intersectSAIntervalWithOutput2(SAIntervalHit& h,
                RapMapSAIndex& rmi,
                SAProcessedHitVec& outStructs);

        /*
        void intersectSAIntervalWithOutput3(SAIntervalHit& h,
                RapMapSAIndex& rmi,
                SAProcessedHitVec& outHits);
                */

        std::vector<ProcessedHit> intersectHits(
                std::vector<HitInfo>& inHits,
                RapMapIndex& rmi);

        SAHitMap intersectSAHits(
                std::vector<SAIntervalHit>& inHits,
                RapMapSAIndex& rmi);

        std::vector<ProcessedSAHit> intersectSAHits2(
                std::vector<SAIntervalHit>& inHits,
                RapMapSAIndex& rmi);
    }
}


#endif // __HIT_MANAGER_HPP__
