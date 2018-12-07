//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015-2018 Rob Patro, Avi Srivastava, Hirak Sarkar
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

#ifndef __RAPMAP_SEGMENT_MAPPING_INFO_HPP__
#define __RAPMAP_SEGMENT_MAPPING_INFO_HPP__

#include <vector>
#include <utility>
#include "nonstd/span.hpp"
#include "cuckoohash_map.hh"
#include "spdlog/spdlog.h"

using SegmentIDType = uint32_t;

struct SegmentCountValue {
  std::array<uint32_t, 8> typeCounts{{0, 0, 0, 0, 0, 0, 0, 0}};
};


// Parse segment information from yanagi (https://github.com/mgunady/yanagi)
class SegmentMappingInfo {
  using TranscriptIDType = int32_t;
public:
  SegmentMappingInfo();
  bool loadFromFile(const std::string& fname,
                    const std::vector<std::string>& segmentNames,
                    std::shared_ptr<spdlog::logger> log);

  void serialize(const std::string& outDir);

  void load(const std::string& outDir);

  size_t numSegments() const { return txpListRanges_.size(); }

  nonstd::span<TranscriptIDType> transcriptsForSegment(int64_t segmentID) const {
    auto& p = txpListRanges_[segmentID];
    return nonstd::span<TranscriptIDType>(static_cast<int32_t*>(const_cast<int32_t*>(txpList_.data()) + p.first), static_cast<size_t>(p.second));
  }

  size_t tableSize() const { return countMap_.size(); }

  void addHit(SegmentIDType s1, SegmentIDType s2, uint8_t mappingType) {
    auto p = std::make_pair(s1, s2);
    auto upfn = [mappingType](SegmentCountValue& x) -> void {
      x.typeCounts[mappingType] += 1;
    };
    SegmentCountValue v;
    countMap_.upsert(p, upfn, v);
  }

private:
    struct pairhash {
    public:
      template <typename T, typename U>
      std::size_t operator()(const std::pair<T, U> &x) const
      {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
      }
    };

    // The segment id provides an
    // interval into the transcript list vector of
    // all transcripts corresponding to this segment.
    std::vector<std::pair<uint32_t, uint32_t>> txpListRanges_;
    std::vector<TranscriptIDType> txpList_;
    std::vector<std::string> txpNames_;
  cuckoohash_map<std::pair<SegmentIDType, SegmentIDType>, SegmentCountValue, pairhash> countMap_;
};

#endif // __RAPMAP_SEGMENT_MAPPING_INFO_HPP__
