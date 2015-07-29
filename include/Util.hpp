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

#ifndef RSDIC_UTIL_HPP_
#define RSDIC_UTIL_HPP_

#include <vector>
#include <stdint.h>
#include "Const.hpp"

namespace rsdic{

class Util{
public:
  static uint64_t GetSlice(const std::vector<uint64_t>& bits,
                           uint64_t pos, uint64_t len) {
    if (len == 0) return 0;
    uint64_t block = pos / kSmallBlockSize;
    uint64_t offset = pos % kSmallBlockSize;
    uint64_t ret = bits[block] >> offset;
    if (offset + len > kSmallBlockSize){
      ret |= (bits[block+1] << (kSmallBlockSize - offset));
    }
    if (len == 64) return ret;
    return ret & ((1LLU << len) - 1);
  }

  static void SetSlice(std::vector<uint64_t>& bits,
                       uint64_t pos, uint64_t len, uint64_t val) {
    if (len == 0) return;
    uint64_t block = pos / kSmallBlockSize;
    uint64_t offset = pos % kSmallBlockSize;
    bits[block] |= val << offset;
    if (offset + len > kSmallBlockSize){
      bits[block+1] |= val >> (kSmallBlockSize - offset);
    }
  }

  static uint64_t Floor(uint64_t num, uint64_t div){
    return (num + div - 1) / div;
  }

  static uint64_t GetNum(bool bit, uint64_t num, uint64_t total) {
    if (bit) return num;
    else return total - num;
  }


};

}

#endif // RSDIC_UTIL_HPP_
