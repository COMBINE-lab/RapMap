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

#ifndef RSDIC_CONST_HPP_
#define RSDIC_CONST_HPP_

#include <stdint.h>

namespace rsdic{

static const uint64_t kLargeBlockSize = 1024;
static const uint64_t kSmallBlockSize = 64;
static const uint64_t kSelectBlockSize = 2048;
static const uint64_t kUseRawLen = 48;
static const uint64_t kSmallBlockPerLargeBlock = kLargeBlockSize / kSmallBlockSize;


}

#endif // RSDIC_CONST_HPP_
