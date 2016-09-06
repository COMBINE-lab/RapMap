#ifndef SPARSEPP_HASH_SERALIZER_HPP
#define SPARSEPP_HASH_SERALIZER_HPP

#include "sparsepp.h"

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
