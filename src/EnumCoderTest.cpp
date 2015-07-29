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

#include <gtest/gtest.h>
#include "EnumCoder.hpp"

using namespace std;
using namespace rsdic;

uint64_t PopCount(uint64_t x){
  uint64_t count = 0;
  for (uint64_t i = 0; i < 64; ++i){
    if ((x >> i) & 1LLU) ++count;
  }
  return count;
}

TEST(EnumCoder, small){
  uint64_t code = EnumCoder::Encode(0, PopCount(0));
  ASSERT_EQ(0, EnumCoder::Decode(0, code));
}

TEST(EnumCoder, random){
  for (uint64_t i = 0; i < 10000; ++i){
    uint64_t x = rand();
    uint64_t rank_sb = PopCount(x);
    uint64_t code = EnumCoder::Encode(x, rank_sb);
    ASSERT_EQ(x, EnumCoder::Decode(code, rank_sb));
  }
}
