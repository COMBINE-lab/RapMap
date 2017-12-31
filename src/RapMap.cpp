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

#include <iostream>
#include <fstream>
#include <vector>

#include <cereal/archives/json.hpp>

#include "RapMapConfig.hpp"
#include "IndexHeader.hpp"

int rapMapIndex(int argc, char* argv[]);
int rapMapSAIndex(int argc, char* argv[]);
int rapMapMap(int argc, char* argv[]);
int rapMapSAMap(int argc, char* argv[]);

void printUsage() {
    std::string versionString = rapmap::version;
    std::cerr << "RapMap Transcriptome Aligner v"
              << versionString << '\n';
    std::cerr << "=====================================\n";
    auto usage =
        R"(
There are currently 4 RapMap subcommands
    quasiindex --- builds a suffix array-based (SA) index
    quasimap   --- map reads using the SA-based index

Run a corresponding command "rapmap <cmd> -h" for
more information on each of the possible RapMap
commands.)";
    std::cerr << usage << '\n';
}

bool isIndexArg(char* arg) {
    std::string argStr(arg);
    return (argStr == "-i") or (argStr == "--index");
}


int main(int argc, char* argv[]) {

    std::vector<char*> args;
    args.push_back(argv[0]);

    if (argc < 2) {
        printUsage();
        std::exit(0);
    }

    for (int i = 2; i < argc; ++i) {
        args.push_back(argv[i]);
    }

    if (std::string(argv[1]) == "-h" or
        std::string(argv[1]) == "--help") {
        printUsage();
        std::exit(0);
    }

    if (std::string(argv[1]) == "pseudoindex") {
	std::cerr << "This mode not enabled in this branch\n";
        return 1;// return rapMapIndex(argc - 1, args.data());
    } else if (std::string(argv[1]) == "quasiindex") {
        return rapMapSAIndex(argc - 1, args.data());
    } else if (std::string(argv[1]) == "pseudomap") {
	std::cerr << "This mode not enabled in this branch\n";
        return 1;//rapMapIndex(argc - 1, args.data());
    } else if (std::string(argv[1]) == "quasimap") {
        return rapMapSAMap(argc - 1, args.data());
    } else {
        std::cerr << "the command " << argv[1]
                  << " is not yet implemented\n";
        return 1;
    }
    return 0;
}
