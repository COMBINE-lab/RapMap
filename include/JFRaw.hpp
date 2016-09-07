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

#ifndef __JF_RAW_H__
#define __JF_RAW_H__

#include "jellyfish/file_header.hpp"
// Type for values
/*
struct value_type {
  char foo;
  int  bar;
  bool baz;
};
*/

// Special header type. Just like the jellyfish header type, but save
// one extra piece of information about the hash array.
class SpecialHeader : public jellyfish::file_header {
public:
  SpecialHeader() = default;
  SpecialHeader(std::istream& is) : jellyfish::file_header(is) { }

  template<typename storage>
  void update_from_ary(const storage& ary) {
    jellyfish::file_header::update_from_ary(ary);
    root_["size_bytes"] = (Json::UInt64)ary.size_bytes();
  }

  size_t size_bytes() const { return root_["size_bytes"].asLargestUInt(); }
};

#endif /* __JF_RAW_H__ */
