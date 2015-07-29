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

#include "RSDic.hpp"
#include "Const.hpp"
#include "Util.hpp"
#include "EnumCoder.hpp"
#include "RSDicBuilder.hpp"

using namespace std;

namespace rsdic{

RSDicBuilder::RSDicBuilder() : buf_(0), offset_(0), bit_num_(0), one_num_(0), prev_one_num_(0), zero_num_(0) {
}

void RSDicBuilder::Clear(){
  bits_.clear();
  pointer_blocks_.clear();
  select_one_inds_.clear();
  select_zero_inds_.clear();
  rank_blocks_.clear();
  rank_small_blocks_.clear();
  buf_ = 0;
  offset_ = 0;
  bit_num_ = 0;
  one_num_ = 0;
  prev_one_num_ = 0;
  zero_num_ = 0;
}

void RSDicBuilder::PushBack(bool bit) {
  if (bit_num_ % kSmallBlockSize == 0){
    WriteBlock();
  }
  if (bit){
    buf_ |= (1LLU << (bit_num_ % kSmallBlockSize));
    if ((one_num_ % kSelectBlockSize) == 0){
      select_one_inds_.push_back(bit_num_ / kLargeBlockSize);
    }
    ++one_num_;
  } else {
    if ((zero_num_ % kSelectBlockSize) == 0){
      select_zero_inds_.push_back(bit_num_ / kLargeBlockSize);
    }
    ++zero_num_;
  }
  ++bit_num_;
}

void RSDicBuilder::WriteBlock(){
  if (bit_num_ > 0) {
    uint64_t rank_sb = one_num_ - prev_one_num_;
    rank_small_blocks_.push_back(rank_sb);
    prev_one_num_ = one_num_;


    uint64_t len = EnumCoder::Len(rank_sb);
    uint64_t code = 0;
    if (len == kSmallBlockSize){
      code = buf_; // use raw 
    } else {
      code = EnumCoder::Encode(buf_, rank_sb);
    }
    uint64_t new_size =  Util::Floor(offset_ + len, kSmallBlockSize);
    if (new_size > bits_.size()) {
      bits_.push_back(0);
    }
    Util::SetSlice(bits_, offset_, len, code);
    buf_ = 0;
    offset_ += len;
  }
  if ((bit_num_ % kLargeBlockSize) == 0){
    rank_blocks_.push_back(one_num_);
    pointer_blocks_.push_back(offset_);
  }
}


void RSDicBuilder::Build(RSDic& bv){
  bv.Clear();
  if (bit_num_ == 0) return;
  WriteBlock();
  bv.num_ = bit_num_;
  bv.one_num_ = one_num_;
  // use copy instead of swap to allocate adequate working space
  bv.bits_ = bits_; 
  bv.select_one_inds_ = select_one_inds_;
  bv.select_zero_inds_ = select_zero_inds_;
  bv.pointer_blocks_ = pointer_blocks_;
  bv.rank_blocks_ = rank_blocks_;
  bv.rank_small_blocks_ = rank_small_blocks_;
}

} // rsdic
