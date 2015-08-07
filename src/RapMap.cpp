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
    pseudoindex   --- builds a k-mer-based index
    pseudomap     --- map reads using a k-mer-based index
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
        return rapMapIndex(argc - 1, args.data());
    } else if (std::string(argv[1]) == "quasiindex") {
        return rapMapSAIndex(argc - 1, args.data());
    } else if (std::string(argv[1]) == "pseudomap") {
        return rapMapMap(argc - 1, args.data());
    } else if (std::string(argv[1]) == "quasimap") {
        return rapMapSAMap(argc - 1, args.data());
    } else {
        std::cerr << "the command " << argv[1]
                  << " is not yet implemented\n";
        return 1;
    }
    return 0;
}
