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

#ifndef __HIT_MANAGER_SA_HPP__
#define __HIT_MANAGER_SA_HPP__

#include "RapMapUtils.hpp"
#include "RapMapIndex.hpp"
#include "RapMapSAIndex.hpp"
#include "FastxParser.hpp"
//#include "eytzinger_array.h"


#include <tuple>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>

namespace rapmap {
    namespace hit_manager_sa {
        using HitInfo = rapmap::utils::HitInfo;
        using ProcessedHit = rapmap::utils::ProcessedHit;
        using MateStatus = rapmap::utils::MateStatus;
        using PositionListHelper = rapmap::utils::PositionListHelper;
        using QuasiAlignment = rapmap::utils::QuasiAlignment;
        using TxpQueryPos = rapmap::utils::TxpQueryPos;
        using SATxpQueryPos = rapmap::utils::SATxpQueryPos;

        template <typename T>
        using SAIntervalHit = rapmap::utils::SAIntervalHit<T>;
        using SAHitMap = std::map<int, rapmap::utils::ProcessedSAHit>;
        using ProcessedSAHit = rapmap::utils::ProcessedSAHit;

        //previously used methods

        // This function is called from
        // SalmonQuantify
        // So we have to include this
        // hpp file there

        //template <typename readT>
	    bool mergeLeftRightMap(fastx_parser::ReadPair& rpair,
							 SAHitMap& leftMap,
							 SAHitMap& rightMap,
                             std::vector<QuasiAlignment>& jointHits
							  );

        template <typename RapMapIndexT>
        bool mergeLeftRightSAInts(
        					    fastx_parser::ReadPair& rpair,
								bool lhp,
								bool rhp,
                                std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftFwdSAInts,
                                std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& leftRcSAInts,
                                std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& rightFwdSAInts,
                                std::vector<SAIntervalHit<typename RapMapIndexT::IndexType>>& rightRcSAInts,
								std::vector<QuasiAlignment>& jointHits,
                                RapMapIndexT& rmi,
								bool maxNumHits,
								bool consistentHits,
                                rapmap::utils::HitCounters& hctr);

    }
}


#endif // __HIT_MANAGER_HPP__
