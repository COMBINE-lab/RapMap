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

#ifndef RSDIC_RSDIC_BUILDER_HPP_
#define RSDIC_RSDIC_BUILDER_HPP_

#include <vector>
#include <stdint.h>
#include "Type.hpp"

namespace rsdic{

class RSDic;

class RSDicBuilder{

public:
  RSDicBuilder();
  void Clear();
  void PushBack(bool bit);
  void Build(RSDic& bitvec);

private:
  void WriteBlock();
  
  std::vector<uint64_t> bits_;
  std::vector<rsdic_uint> pointer_blocks_;
  std::vector<rsdic_uint> select_one_inds_;
  std::vector<rsdic_uint> select_zero_inds_;
  std::vector<rsdic_uint> rank_blocks_;
  std::vector<uint8_t> rank_small_blocks_;
  uint64_t buf_;
  rsdic_uint offset_;
  rsdic_uint bit_num_;
  rsdic_uint one_num_;
  rsdic_uint prev_one_num_;
  rsdic_uint zero_num_;
};

}

#endif // RSDIC_RSDIC_BUILDER_HPP_
