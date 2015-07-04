#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>

#include "xxhash.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include "PairSequenceParser.hpp"

using paired_parser = pair_sequence_parser<char**>;
using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

class KmerKeyHasher {
    public:
        size_t operator()(const my_mer& m) const {
            auto k = my_mer::k();
            auto v = m.get_bits(0, 2*k);
            return XXH64(static_cast<void*>(&v), 8, 0);
        }
};

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT>//, typename CoverageCalculator>
void processTranscripts(ParserT* parser,
                        std::mutex& iomutex) {
    // Seed with a real random value, if available
    std::random_device rd;

    // Create a random uniform distribution
    std::default_random_engine eng(rd());

    std::uniform_int_distribution<> dis(0, 3);

    uint32_t n{0};
    uint32_t k = my_mer::k();
    std::vector<std::string> transcriptNames;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};

    using TranscriptList = std::vector<uint32_t>;
    std::unordered_map<my_mer, TranscriptList, KmerKeyHasher> hmap;
  //size_t locRead{0};
  //uint64_t localUpperBoundHits{0};
  while(true) {
    typename ParserT::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;           // If got nothing, quit

    for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
        std::string& readStr = j->data[i].seq;
        uint32_t readLen  = readStr.size();
        uint32_t txpIndex = n++;
        transcriptNames.push_back(j->data[i].header);
        my_mer mer;
        mer.polyT();
        for (size_t b = 0; b < readLen; ++b) {
            int c = jellyfish::mer_dna::code(readStr[b]);
            if (jellyfish::mer_dna::not_dna(c)) {
                c = jellyfish::mer_dna::code(bases[dis(eng)]);
            }
            mer.shift_left(c);

            if (b >= k) {
                //auto key = mer.get_bits(0, 2*k);
                hmap[mer].emplace_back(txpIndex);
            }
        }
        if (n % 10000 == 0) {
            std::cerr << "hashed k-mers for " << n << " transcripts\n";
        }
    }
  }
  std::cerr << "parsed " << transcriptNames.size() << " transcripts\n";

}



int rapMapIndex(int argc, char* argv[]) {
    std::cerr << "indexing\n";
    // stupid parsing for now
    std::string reads1(argv[0]);
    std::vector<std::string> readFiles({ reads1 });

    int k = std::atoi(argv[1]);
    my_mer::k(static_cast<uint32_t>(k));

    size_t maxReadGroup{1000}; // Number of reads in each "job"
    size_t concurrentFile{2}; // Number of files to read simultaneously
    size_t numThreads{2};
    stream_manager streams(readFiles.begin(), readFiles.end(), concurrentFile);
    std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
    transcriptParserPtr.reset(new single_parser(4 * numThreads, maxReadGroup,
                              concurrentFile, streams));
    std::mutex iomutex;
    processTranscripts(transcriptParserPtr.get(), iomutex);
    return 0;
}


