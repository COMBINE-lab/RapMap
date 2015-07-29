/*
 *  Copyright (c) 2012 Daisuke Okanohara
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#ifndef RSDIC_RSDIC_HPP_
#define RSDIC_RSDIC_HPP_

#include <vector>
#include <iostream>
#include <stdint.h>
#include "Type.hpp"

namespace rsdic{

class RSDic{
public:
  RSDic();
  ~RSDic();

  void Clear();
  bool GetBit(uint64_t pos) const;
  uint64_t Rank(uint64_t pos, bool bit) const;
  uint64_t Select(uint64_t ind, bool bit) const;
  std::pair<uint64_t, uint64_t> GetBitAndRank(uint64_t pos) const;
  void Save(std::ostream& os) const;
  void Load(std::istream& is);
  uint64_t GetUsageBytes() const;

  uint64_t num() const{
    return num_;
  }

  uint64_t one_num() const{
    return one_num_;
  }

  uint64_t zero_num() const{
    return num_ - one_num_;
  }

  bool operator == (const RSDic& bv) const;


private:
  friend class RSDicBuilder;

  uint64_t Select1(uint64_t ind) const;
  uint64_t Select0(uint64_t ind) const;

  template <class T>
  void Save(std::ostream& os, const std::vector<T>& vs) const{
    uint64_t size = vs.size();
    os.write((const char*)&size, sizeof(size));
    os.write((const char*)&vs[0], sizeof(vs[0]) * size);
  }

  template <class T>
  void Load(std::istream& is, std::vector<T>& vs){
    uint64_t size = 0;
    is.read((char*)&size, sizeof(size));
    vs.resize(size);
    is.read((char*)&vs[0], sizeof(vs[0]) * size);
  }

  std::vector<uint64_t> bits_;
  std::vector<rsdic_uint> pointer_blocks_;
  std::vector<rsdic_uint> rank_blocks_;
  std::vector<rsdic_uint> select_one_inds_;
  std::vector<rsdic_uint> select_zero_inds_;
  std::vector<uint8_t> rank_small_blocks_;
  uint64_t num_;
  uint64_t one_num_;
};

} // rsdic

#endif // RSDIC_BITVEC_HPP_
