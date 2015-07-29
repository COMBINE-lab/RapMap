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
#include "RSDicBuilder.hpp"

#define private public
#include "RSDic.hpp"
#include "EnumCoder.hpp"


using namespace std;
using namespace rsdic;

TEST(BitVec, combination){
  for (uint64_t i = 0; i <= 64; ++i){
    ASSERT_EQ(1, EnumCoder::kCombinationTable64_[i][i]);
  }
  for (uint64_t i = 1; i <= 64; ++i){
    for (uint64_t j = 1; j <= i; ++j){
      ASSERT_EQ(EnumCoder::kCombinationTable64_[i][j],
                EnumCoder::kCombinationTable64_[i-1][j-1] +
                EnumCoder::kCombinationTable64_[i-1][j]);
    }
  }
}


TEST(RSDic, small){
  RSDicBuilder bvb;
  const uint64_t n = 65;
  for (int i = 0; i < n; ++i){
    bvb.PushBack(1);
  }

  RSDic bv;
  bvb.Build(bv);
  ASSERT_EQ(n, bv.num());
  ASSERT_EQ(n, bv.one_num());
  for (size_t i = 0; i < bv.num(); ++i){
    ASSERT_EQ(1, bv.GetBit(i));
    ASSERT_EQ(i, bv.Rank(i, 1));
    ASSERT_EQ(i, bv.Select(i, 1));
  }
}


TEST(RSDic, trivial_zero){
  RSDicBuilder bvb;
  const uint64_t n = 10000;
  for (int i = 0; i < n; ++i){
    bvb.PushBack(0);
  }

  RSDic bv;
  bvb.Build(bv);
  ASSERT_EQ(n, bv.num());
  ASSERT_EQ(0, bv.one_num());
  for (size_t i = 0; i < bv.num(); ++i){
    ASSERT_EQ(0, bv.GetBit(i));
    ASSERT_EQ(i, bv.Rank(i, 0));
    ASSERT_EQ(i, bv.Select(i, 0)) << " i=" << i;
  }
}

TEST(RSDic, trivial_one){
  RSDicBuilder bvb;
  const uint64_t n = 10000;
  for (int i = 0; i < n; ++i){
    bvb.PushBack(1);
  }

  RSDic bv;
  bvb.Build(bv);
  ASSERT_EQ(n, bv.num());
  ASSERT_EQ(n, bv.one_num());
  for (size_t i = 0; i < bv.num(); ++i){
    ASSERT_EQ(1, bv.GetBit(i));
    ASSERT_EQ(i, bv.Rank(i, 1)) << " i=" << i;
    ASSERT_EQ(i, bv.Select(i, 1)) << " i=" << i;
  }
}


TEST(RSDic, random){
  RSDicBuilder bvb;
  vector<int> B;
  const uint64_t n = 100000;
  for (int i = 0; i < n; ++i){
    int b = rand() % 2;
    bvb.PushBack(b);
    B.push_back(b);
  }

  RSDic bv;
  bvb.Build(bv);
  ASSERT_EQ(n, bv.num());
  int sum = 0;
  for (size_t i = 0; i < bv.num(); ++i){
    ASSERT_EQ(B[i]  , bv.GetBit(i)) << " i=" << i;
    pair<uint64_t, uint64_t> bit_rank = bv.GetBitAndRank(i);
    ASSERT_EQ(B[i], bit_rank.first);
    if (B[i]){
      ASSERT_EQ(sum, bv.Rank(i, 1));
      ASSERT_EQ(sum, bit_rank.second);
      ASSERT_EQ(i,bv.Select(sum, 1));
    } else {
      ASSERT_EQ(i - sum, bv.Rank(i, 0));
      ASSERT_EQ(i - sum, bit_rank.second);
      ASSERT_EQ(i, bv.Select(i-sum, 0));
    }

    sum += B[i];
  }

  ostringstream os;
  bv.Save(os);
  istringstream is(os.str());

  RSDic bv_load;
  bv_load.Load(is);

  ASSERT_EQ(bv, bv_load);
}

TEST(RSDic, large){
  RSDicBuilder rsdb;
  const uint64_t n = 26843545;
  vector<uint64_t> poses;
  for (uint64_t i = 0; i < n; ++i){
    float r = (float)rand() / RAND_MAX;
    if (r < 0.001) {
      rsdb.PushBack(1);
      poses.push_back(i);
    }
    else rsdb.PushBack(0);
  }

  RSDic bv;
  rsdb.Build(bv);
  uint64_t one_num = bv.one_num();
  for (uint64_t i = 0; i < one_num; ++i){
    ASSERT_EQ(poses[i], bv.Select(i, 1));
  }
}

