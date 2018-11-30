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

#ifndef SPARSEPP_HASH_SERALIZER_HPP
#define SPARSEPP_HASH_SERALIZER_HPP

#include "sparsepp/spp.h"

namespace spp_utils {

// Can not be used with datatypes containing internal pointers
// OK with POD datatypes
template <typename key_type, typename value_type> struct pod_hash_serializer {
  using KeySerializer = spp::sparsehash_internal::pod_serializer<key_type>;
  using ValueSerializer = spp::sparsehash_internal::pod_serializer<value_type>;

  KeySerializer ks_;
  ValueSerializer vs_;

  template <typename OUTPUT>
  bool operator()(OUTPUT* fp, const std::pair<const key_type, value_type>& value) const {
    return ks_(fp, value.first) && vs_(fp, value.second);
  }

  template <typename INPUT>
  bool operator()(INPUT* fp, std::pair<const key_type, value_type>* value) const {
    return ks_(fp, (key_type*)&value->first) && vs_(fp, (value_type*)&value->second);
  }
};

}

#endif // SPARSEPP_HASH_SERALIZER_HPP
