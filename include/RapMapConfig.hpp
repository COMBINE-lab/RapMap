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

#ifndef __RAPMAP_CONFIG_HPP__
#define __RAPMAP_CONFIG_HPP__

#include <string>

namespace rapmap {
    constexpr char majorVersion[] = "0";
    constexpr char minorVersion[] = "6";
    constexpr char patchVersion[] = "0";
    constexpr char version [] = "0.6.0";
    constexpr uint32_t indexVersion = 3;
}

#endif //__RAPMAP_CONFIG_HPP__
