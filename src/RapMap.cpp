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
#include <clocale>

#include <cereal/archives/json.hpp>

#include "RapMapConfig.hpp"
#include "IndexHeader.hpp"

int rapMapIndex(int argc, char* argv[]);
int rapMapSAIndex(int argc, char* argv[]);
int rapMapMap(int argc, char* argv[]);
int rapMapSAMap(int argc, char* argv[]);

void printUsage() {
    std::string versionString = rapmap::version;
    std::cerr << "rapmap v"
              << versionString << '\n';
    std::cerr << "=====================================\n";
    auto usage =
        R"(
There are currently 2 rapmap subcommands
    quasiindex --- builds a quasi-mapping index on the reference
    quasimap   --- maps reads to a quasi-mapping index

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
    std::setlocale(LC_NUMERIC, "en_US.UTF-8");
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

    if (std::string(argv[1]) == "-v" or
        std::string(argv[1]) == "--version") {
      std::cout << "rapmap " << rapmap::version << "\n";
      std::exit(0);
    }

    if (std::string(argv[1]) == "quasiindex") {
        return rapMapSAIndex(argc - 1, args.data());
    } else if (std::string(argv[1]) == "quasimap") {
        return rapMapSAMap(argc - 1, args.data());
    } else {
        std::cerr << "the command " << argv[1]
                  << " is not yet implemented\n";
        return 1;
    }
    return 0;
}
