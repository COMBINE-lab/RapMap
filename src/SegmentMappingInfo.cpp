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

#include <fstream>
#include <map>
#include <utility>
#include <sstream>
#include "sparsepp/spp.h"
#include "SegmentMappingInfo.hpp"
#include "RapMapUtils.hpp"

SegmentMappingInfo::SegmentMappingInfo() {}

bool SegmentMappingInfo::loadFromFile(const std::string& fname,
                                      const std::vector<std::string>& segmentNames,
                                      std::shared_ptr<spdlog::logger> log
                                      ) {
  /*
    segID   chrom   geneID  txAnnIDs        binIDs  st      end     strand
    SEG0000001      10      ENSG00000012779 ENST00000542434 57010,57011     45869661        45869774        +
    SEG0000002      10      ENSG00000012779 ENST00000374391,ENST00000542434 57011,57012,57013,57014,57015,57016,57017       45869675        45920450        +
    SEG0000003      10      ENSG00000012779 ENST00000374391,ENST00000483623,ENST00000542434 57016,57017     45919539        45920580        +
    SEG0000004      10      ENSG00000012779 ENST00000483623 57017,57018     45920481        45923934        +
  */

  std::ifstream ifile(fname);
  std::string line;
  // The first line is the header
  std::getline(ifile, line);
  auto header = rapmap::utils::tokenize(line, '\t');
  std::map<std::string, int32_t> hmap;
  for (size_t i = 0; i < header.size(); ++i) {
    hmap[header[i]] = i;
  }

  spp::sparse_hash_map<std::string, uint32_t> segmentID;
  segmentID.reserve(segmentNames.size());

  spp::sparse_hash_map<std::string, uint32_t> txpID;

  std::vector<bool> seenSegment(segmentNames.size(), false);

  uint32_t ctr{0};
  for (const auto& sn : segmentNames) {
    segmentID[sn] = ctr;
    ++ctr;
  }

  // number of segments we care about
  //txpListRanges_.resize(segmentNames.size());

  int32_t segID = hmap["segID"];
  int32_t txAnnIDs = hmap["txAnnIDs"];

  std::vector<std::pair<uint32_t, uint32_t>> segmentIntervals;

  uint64_t lastIdx{0};
  while (std::getline(ifile, line)) {
    auto tokens = rapmap::utils::tokenize(line, '\t');
    std::string sid = tokens[segID];
    std::string txlist = tokens[txAnnIDs];

    auto sidIt = segmentID.find(sid);
    uint32_t segIdx{std::numeric_limits<uint32_t>::max()};
    if (sidIt == segmentID.end()) {
      //log->error("Segment mapping {} referred to a segment {} not in the sequence file.", fname, sid);
      continue;
      //return false;
    } else {
      segIdx = sidIt->second;
    }

    if (seenSegment[segIdx]) {
      log->error("Segment {} was seen more than once in the segment mapping {}.", sidIt->first, fname);
      return false;
    }
    seenSegment[segIdx] = true;

    auto txtokens = rapmap::utils::tokenize(txlist, ',');
    auto firstIdx = lastIdx;
    uint32_t len{0};
    for (auto& tx : txtokens) {
      // index for this transcript
      auto txIt = txpID.find(tx);
      uint32_t txpIdx{std::numeric_limits<uint32_t>::max()};
      if (txIt == txpID.end()) {
        txpIdx = txpID.size();
        txpID[tx] = txpIdx;
        txpNames_.push_back(tx);
      } else {
        txpIdx = txIt->second;
      }
      txpList_.push_back(txpIdx);
      ++lastIdx;
      ++len;
    }
    segmentIntervals.push_back(std::make_pair(firstIdx, len));
    //std::cerr << "sid : " << sid << " :: txlist : " << txlist << "\n";
  }

  txpListRanges_.resize(segmentIntervals.size());
  size_t siIdx{0};
  for (auto& si : segmentIntervals) {
    txpListRanges_[siIdx] = nonstd::span<int32_t>(txpList_.data() + si.first, si.second);
    ++siIdx;
  }

  log->info("Total txps {} for {} segments.", txpNames_.size(), txpListRanges_.size());
  log->info("Total txps list length = {}", txpList_.size());

  siIdx = 0;
  for (auto r : txpListRanges_) {
    std::stringstream s;
    s << "transcripts for " << siIdx << " : {";
    for (auto tid : r) {
      s << tid << ',';
    }
    s << "}\n";
    log->info(s.str());
    siIdx++;
  }
  ifile.close();
  return true;
}

