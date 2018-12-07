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
#include "nonstd/span.hpp"
#include "spdlog/spdlog.h"

// Parse segment information from yanagi (https://github.com/mgunady/yanagi)
class SegmentMappingInfo {
  using TranscriptIDType = int32_t;
public:
  SegmentMappingInfo();
  bool loadFromFile(const std::string& fname,
                    const std::vector<std::string>& segmentNames,
                    std::shared_ptr<spdlog::logger> log);
  nonstd::span<TranscriptIDType> transcriptsForSegment(int64_t segmentID) const;
 
private:
  // The segment id provides an
  // interval into the transcript list vector of
  // all transcripts corresponding to this segment.
  std::vector<nonstd::span<TranscriptIDType>> txpListRanges_;
  std::vector<TranscriptIDType> txpList_;
  std::vector<std::string> txpNames_;
};

#endif // __RAPMAP_SEGMENT_MAPPING_INFO_HPP__
