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
