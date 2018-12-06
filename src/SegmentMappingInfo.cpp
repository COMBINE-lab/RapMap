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
#include "SegmentMappingInfo.hpp"
#include "RapMapUtils.hpp"

SegmentMappingInfo::SegmentMappingInfo() {}

bool SegmentMappingInfo::loadFromFile(const std::string& fname) {
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

  int32_t segID = hmap["segID"];
  int32_t txAnnIDs = hmap["txAnnIDs"];
  while (std::getline(ifile, line)) {
    auto tokens = rapmap::utils::tokenize(line, '\t');
    std::string sid = tokens[segID];
    std::string txlist = tokens[txAnnIDs];
    std::cerr << "sid : " << sid << " :: txlist : " << txlist << "\n";
  }

  ifile.close();
  return true;
}

