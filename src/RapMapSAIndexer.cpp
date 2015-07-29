#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <type_traits>
#include <fstream>
#include <cctype>

#include "tclap/CmdLine.h"

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "xxhash.h"
#include "btree/btree_map.h"

#include "spdlog/spdlog.h"

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

#include "sais.h"
#include "RSDic.hpp"
#include "RSDicBuilder.hpp"
#include "RapMapUtils.hpp"
#include "RapMapFileSystem.hpp"
#include "ScopedTimer.hpp"

#include "jellyfish/file_header.hpp"
#include "jellyfish/binary_dumper.hpp"
#include "jellyfish/thread_exec.hpp"
#include "jellyfish/hash_counter.hpp"
#include "jellyfish/mer_overlap_sequence_parser.hpp"
#include "jellyfish/mer_iterator.hpp"
#include "JFRaw.hpp"

#include <chrono>

using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using MerMapT = jellyfish::cooperative::hash_counter<rapmap::utils::my_mer>;

//uint64_t encode(uint32_t tid, uint32_t pos) {
//    uint64_t res = tid;
//    res = res << 32;
//    res |= pos;
//    return res;
//}
//
//constexpr uint32_t uint32HighestBitMask = 0x80000000;
//constexpr uint32_t uint31LowBitMask = 0x7FFFFFFF;
//
//constexpr uint32_t uint30LowBitMask = 0x3FFFFFFF;
//
//constexpr uint32_t rcSetMask = 0x40000000;
//
//// marks the second highest bit
//void markRCBit(uint32_t& i) { i |= rcSetMask; }
//
//// marks the highest bit
//void markNewTxpBit(uint32_t& i) { i |= uint32HighestBitMask; }
//
//bool wasSeen(uint32_t i) { return ((i & uint32HighestBitMask) >> 31) == 1; }
//
//void markSeen(uint32_t& i) { i |= uint32HighestBitMask; }
//
//uint32_t unmarked(uint32_t i) { return (i & uint31LowBitMask); }
//
//class VectorHasher {
//    public:
//    size_t operator()(const std::vector<uint32_t>& vec) const {
//        size_t hashVal{0};
//        for (auto v : vec) {
//            rapmap::utils::hashCombine(hashVal, v);
//        }
//        return hashVal;
//    }
//};
//
//struct PosInfo {
//    PosInfo(uint64_t merIDIn, bool isRCIn, uint32_t posIn) :
//        merID(merIDIn), isRC(isRCIn), pos(posIn) {}
//
//    uint64_t merID;
//    bool isRC;
//    uint32_t pos;
//};
//
//// maybe unify with the above?
//struct JumpCell {
//    JumpCell(uint32_t merIdxIn, int32_t posIn, bool isRCIn) :
//        merIdx(merIdxIn), pos(posIn), isRC(isRCIn) {}
//    uint32_t merIdx;
//    int32_t pos;
//    bool isRC;
//};
//
//template <typename T>
//void printVec(std::vector<T>& vec) {
//    std::cerr << "{ ";
//    for (auto& e : vec) {
//        std::cerr << e << ", ";
//    }
//    std::cerr << "}";
//}
//
//// There may be a better way to do this, but here we just check the
//// possible neighbors
//bool isBreakpoint(MerMapT& merIntMap, rapmap::utils::my_mer canonicalMer) {
//    uint32_t inDegree{0};
//    uint32_t outDegree{0};
//    const auto& ary = merIntMap.ary();
//    // extend forward
//    for (char b : {'A', 'C', 'G', 'T'}) {
//        auto newMer = canonicalMer;
//        newMer.shift_left(b);
//        newMer.canonicalize();
//        if (outDegree > 1) { return true; }
//    }
//    // extend backward
//    for (char b : {'A', 'C', 'G', 'T'}) {
//        auto newMer = canonicalMer;
//        newMer.shift_right(b);
//        newMer.canonicalize();
//        inDegree += ary->has_key(newMer) ? 1 : 0;
//        if (inDegree > 1) { return true; }
//    }
//    return false;
//}
//
//void emptyJumpQueue(std::vector<JumpCell>& jumpQueue, int32_t lastBreak,
//                    int32_t pos,
//                    MerMapT& merIntMap,
//                    std::vector<uint8_t>& fwdJump,
//                    std::vector<uint8_t>& revJump) {
//
//    // The maximum representable jump
//    constexpr auto maxJump =
//        std::numeric_limits<std::remove_reference<decltype(fwdJump[0])>::type>::max();
//
//    while (!jumpQueue.empty()) {
//        auto& jumpCell = jumpQueue.back();
//        uint8_t revJumpDist = static_cast<uint8_t>(
//                std::min(jumpCell.pos - lastBreak + 1,
//                static_cast<int32_t>(maxJump)));
//        uint8_t fwdJumpDist = static_cast<uint8_t>(
//                std::min(pos - jumpCell.pos + 1,
//                static_cast<int32_t>(maxJump)));
//        fwdJump[jumpCell.merIdx] =
//            std::min(fwdJumpDist, fwdJump[jumpCell.merIdx]);
//        revJump[jumpCell.merIdx] =
//            std::min(revJumpDist, revJump[jumpCell.merIdx]);
//        // Now that we marked the jump, no need to keep this around.
//        jumpQueue.pop_back();
//    }
//}
//


//
// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT>//, typename CoverageCalculator>
void indexTranscriptsSA(ParserT* parser,
            		std::string& outputDir,
                        std::mutex& iomutex) {
    // Seed with a real random value, if available
    std::random_device rd;

    // Create a random uniform distribution
    std::default_random_engine eng(rd());

    std::uniform_int_distribution<> dis(0, 3);

    uint32_t n{0};
    uint32_t k = rapmap::utils::my_mer::k();
    std::vector<std::string> transcriptNames;
    std::vector<uint32_t> transcriptStarts;
    std::vector<uint32_t> positionIDs;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};
    uint32_t polyAClipLength{10};
    uint32_t numPolyAsClipped{0};
    std::string polyA(polyAClipLength, 'A');

    using TranscriptList = std::vector<uint32_t>;
    using eager_iterator = MerMapT::array::eager_iterator;
    using KmerBinT = uint64_t;

    size_t numDistinctKmers{0};
    size_t numKmers{0};
    size_t currIndex{1};
    std::cerr << "\n[Step 1 of 4] : counting k-mers\n";

    fmt::MemoryWriter txpSeqStream;
    {
        ScopedTimer timer;
        while(true) {
            typename ParserT::job j(*parser);
            if(j.is_empty()) break;
            for(size_t i = 0; i < j->nb_filled; ++i) { // For each sequence
                std::string& readStr = j->data[i].seq;
                readStr.erase(std::remove_if(readStr.begin(), readStr.end(),
                               [](const char a) -> bool {
                                    return !(isprint(a));
                                }), readStr.end());
                // Do Kallisto-esque clipping of polyA tails
                if (readStr.size() > polyAClipLength and
                        readStr.substr(readStr.length() - polyAClipLength) == polyA) {

                    auto newEndPos = readStr.find_last_not_of("Aa");
                    // If it was all As
                    if (newEndPos == std::string::npos) {
                        readStr.resize(0);
                    } else {
                        readStr.resize(newEndPos + 1);
                    }
                    ++numPolyAsClipped;
                }

                uint32_t readLen  = readStr.size();
                uint32_t txpIndex = n++;
                transcriptNames.push_back(j->data[i].header);
                transcriptStarts.push_back(currIndex);

                bool firstBase{true};
                rapmap::utils::my_mer mer;
                mer.polyT();
                for (size_t b = 0; b < readLen; ++b) {
                    readStr[b] = ::toupper(readStr[b]);
                    int c = jellyfish::mer_dna::code(readStr[b]);
                    // Replace non-ACGT bases with pseudo-random bases
                    if (jellyfish::mer_dna::not_dna(c)) {
                        char rbase = bases[dis(eng)];
                        c = jellyfish::mer_dna::code(rbase);
                        readStr[b] = rbase;
                    }
                    positionIDs.push_back(txpIndex);
                }
                txpSeqStream << readStr;
                txpSeqStream << '$';
                positionIDs.push_back(txpIndex);
                currIndex += readLen + 1;
            }
            if (n % 10000 == 0) {
                std::cerr << "\r\rcounted k-mers for " << n << " transcripts";
            }
        }
    }
    std::cerr << "\n";

    std::cerr << "Clipped poly-A tails from " << numPolyAsClipped << " transcripts\n";

    /*
    std::ofstream rsStream(outputDir + "rsd.bin", std::ios::binary);
    {
        ScopedTimer timer;
        std::cerr << "Building rank-select dictionary and saving to disk ";
        rsdic::RSDic rsd;
        rsdb.Build(rsd);
        rsd.Save(rsStream);
        std::cerr << "done\n";
    }
    rsStream.close();
    */

    std::string concatText = txpSeqStream.str();
    txpSeqStream.clear();
    std::ofstream seqStream(outputDir + "txpInfo.bin", std::ios::binary);
    {
        ScopedTimer timer;
        std::cerr << "Writing sequence data to file . . . ";
        cereal::BinaryOutputArchive seqArchive(seqStream);
        seqArchive(transcriptNames);
        seqArchive(transcriptStarts);
        seqArchive(positionIDs);
        seqArchive(concatText);
        std::cerr << "done\n";
    }
    seqStream.close();

    // clear stuff we no longer need
    positionIDs.clear();
    positionIDs.shrink_to_fit();
    transcriptStarts.clear();
    transcriptStarts.shrink_to_fit();
    transcriptNames.clear();
    transcriptNames.shrink_to_fit();
    // done clearing


    // Build the suffix array
    size_t tlen = concatText.length();
    std::vector<int> SA(tlen, 0);
    std::ofstream saStream(outputDir + "sa.bin", std::ios::binary);
    {
        ScopedTimer timer;
        std::cerr << "Building suffix array ";
        auto ret = sais(reinterpret_cast<unsigned char*>(
                        const_cast<char*>(concatText.c_str())),
                        SA.data(), tlen + 1);
        if (ret == 0) {
            std::cerr << "SUCCESS!\n";
            {
                ScopedTimer timer2;
                std::cerr << "saving to disk . . . ";
                cereal::BinaryOutputArchive saArchive(saStream);
                saArchive(SA);
		// don't actually need the LCP right now
                // saArchive(LCP);
                std::cerr << "done\n";
            }
        } else {
            std::cerr << "FAILURE: return code was " << ret << "\n";
        }
        std::cerr << "done\n";
    }
    saStream.close();

    // clear things we don't need
    //LCP.clear();
    // LCP.shrink_to_fit();
    // done clearing

    // Now, build the k-mer lookup table
    std::unordered_map<uint64_t,
                       rapmap::utils::SAInterval,
                       rapmap::utils::KmerKeyHasher> khash;
    /*
    concatText.erase(std::remove_if(concatText.begin(),
                                    concatText.end(),
                                    [] (const char a) -> bool { return !isprint(a); }),
                     concatText.end());
                     */

    // The start and stop of the current interval
    uint32_t start = 0, stop = 0;
    // An iterator to the beginning of the text
    auto textB = concatText.begin();
    auto textE = concatText.end();
    // The current k-mer as a string
    rapmap::utils::my_mer mer;
    bool currentValid{false};
    std::string currentKmer;
    std::string nextKmer;
    while (stop < tlen) {
        // Check if the string starting at the
        // current position is valid (i.e. doesn't contain $)
        // and is <= k bases from the end of the string
        nextKmer = concatText.substr(SA[stop], k);
        if (nextKmer.length() == k and
            nextKmer.find_first_of('$') == std::string::npos) {
            // If this is a new k-mer, then hash the current k-mer
            if (nextKmer != currentKmer) {
                if (currentKmer.length() == k and
                    currentKmer.find_first_of('$') == std::string::npos) {
                    mer = rapmap::utils::my_mer(currentKmer);
                    auto bits = mer.get_bits(0, 2*k);
                    auto hashIt = khash.find(bits);
                    if (hashIt == khash.end()) {
                        if (start > 1) {
                            if (concatText.substr(SA[start-1], k) ==
                                concatText.substr(SA[start], k)) {
                                std::cerr << "T[SA["
                                          << start-1 << "]:" << k << "] = "
                                          << concatText.substr(SA[start-1], k)
                                          << " = T[SA[" << start << "]:" << k << "]\n";
                                std::cerr << "start = " << start << ", stop = " << stop << "\n";
                                std::cerr << "(1) THIS SHOULD NOT HAPPEN\n";
                                std::exit(1);
                            }
                        }
                        if (start == stop) {
                            std::cerr << "AHH (1) : Interval is empty! (start = " << start
                                      << ") = (stop =  " << stop << ")\n";
                        }
                        if (start == stop) {
                            std::cerr << "AHH (2) : Interval is empty! (start = " << start
                                << ") = (stop =  " << stop << ")\n";
                        }

                        khash[bits] = {start, stop};
                    } else {
                        std::cerr << "\nERROR (1): trying to add same suffix "
                                  << currentKmer << " (len = "
                                  << currentKmer.length() << ") multiple times!\n";
                        auto prevInt = hashIt->second;
                        std::cerr << "existing interval is ["
                                  << prevInt.begin << ", " << prevInt.end << ")\n";
                        for (auto x = prevInt.begin; x < prevInt.end; ++x) {
                            auto suff = concatText.substr(SA[x], k);
                            for (auto c : suff) {
                                std::cerr << "*" << c << "*";
                            }
                            std::cerr << " (len = " << suff.length() <<")\n";
                        }
                        std::cerr << "new interval is ["
                                  << start << ", " << stop << ")\n";
                        for (auto x = start; x < stop; ++x) {
                            auto suff = concatText.substr(SA[x], k);
                            for (auto c : suff) {
                                std::cerr << "*" << c << "*";
                            }
                            std::cerr << "\n";
                        }
                    }
                }
                currentKmer = nextKmer;
                start = stop;
            }
        } else {
            // If this isn't a valid suffix (contains a $)

            // If the previous interval was valid, put it
            // in the hash.
            if (currentKmer.length() == k and
                currentKmer.find_first_of('$') == std::string::npos) {
                mer = rapmap::utils::my_mer(currentKmer);
                auto bits = mer.get_bits(0, 2*k);
                auto hashIt = khash.find(bits);
                if (hashIt == khash.end()) {
                    if (start > 2) {
                        if (concatText.substr(SA[start-1], k) ==
                            concatText.substr(SA[start], k)) {
                            std::cerr << "T[SA["
                                << start-1 << "]:" << k << "] = "
                                << concatText.substr(SA[start-1], k)
                                << " = T[SA[" << start << "]:" << k << "]\n";
                            std::cerr << "start = " << start << ", stop = " << stop << "\n";
                            std::cerr << "(2) THIS SHOULD NOT HAPPEN\n";
                            std::exit(1);
                        }
                    }
                    khash[bits] = {start, stop};
                } else {
                    std::cerr << "\nERROR (2): trying to add same suffix "
                              << currentKmer << "multiple times!\n";
                    auto prevInt = hashIt->second;
                    std::cerr << "existing interval is ["
                        << prevInt.begin << ", " << prevInt.end << ")\n";
                    for (auto x = prevInt.begin; x < prevInt.end; ++x) {
                        std::cerr << concatText.substr(SA[x], k) << "\n";
                    }
                    std::cerr << "new interval is ["
                        << start << ", " << stop << ")\n";
                    for (auto x = start; x < stop; ++x) {
                        std::cerr << concatText.substr(SA[x], k) << "\n";
                    }
                }

            }
            // The current interval is invalid and empty
            currentKmer = nextKmer;
            start = stop;
        }
        if (stop % 1000000 == 0) {
            std::cerr << "\r\rprocessed " << stop << " positions";
        }
        // We always update the end position
        ++stop;
    }
    if (start < tlen) {
        if (currentKmer.length() == k and
            currentKmer.find_first_of('$') != std::string::npos) {
            mer = rapmap::utils::my_mer(currentKmer);
            khash[mer.get_bits(0, 2*k)] = {start, stop};
        }
    }
    std::cerr << "\nkhash had " << khash.size() << " keys\n";
    std::ofstream hashStream(outputDir + "hash.bin", std::ios::binary);
    {
        ScopedTimer timer;
        std::cerr << "saving hash to disk . . . ";
        cereal::BinaryOutputArchive hashArchive(hashStream);
        hashArchive(k);
        hashArchive(khash);
        std::cerr << "done\n";
    }
    hashStream.close();
    std::exit(1);
//    std::ofstream txpLenStream(outputDir + "txplens.bin", std::ios::binary);
//    {
//        cereal::BinaryOutputArchive txpLenArchive(txpLenStream);
//        txpLenArchive(transcriptLengths);
//    }
//    txpLenStream.close();
//    transcriptLengths.clear();
//    transcriptLengths.shrink_to_fit();
//
//    constexpr uint32_t uint32Invalid = std::numeric_limits<uint32_t>::max();
//    std::vector<uint32_t> transcriptIDs(numKmers, uint32Invalid);
//
//    std::cerr << "\n[Step 2 of 4] : marking k-mers\n";
//    // Mark the transcript in which each occurence oc a k-mer appears
//    // in the transcriptIDs vector.
//
//
//    bool isRC{false};
//    int32_t pos{0};
//    uint32_t offset{0};
//    uint32_t transcriptID{0};
//    {
//    ScopedTimer timer;
//
//
//    for (auto& transcriptSeq : transcriptSeqs) {
//        auto readLen = transcriptSeq.length();
//        rapmap::utils::my_mer mer;
//        mer.polyT();
//        std::vector<KmerBinT> kmers;
//        for (size_t b = 0; b < readLen; ++b) {
//            int c = jellyfish::mer_dna::code(transcriptSeq[b]);
//            mer.shift_left(c);
//            if (b >= k) {
//                auto canonicalMer = mer.get_canonical();
//                uint64_t kmerIndex;
//                auto found = merIntMap.ary()->get_val_for_key(canonicalMer, &kmerIndex);
//                // Should ALWAYS find the key
//                assert(found);
//
//
//                auto& v = kmerInfos[kmerIndex];
//                // use the highest bit to mark if we've seen this k-mer yet or not
//                // If we haven't seen this k-mer yet
//                if (!wasSeen(v.count)) {
//                   // Where we start looking for transcripts for this k-mer
//                   v.offset = offset;
//                   offset += v.count;
//                   // The number of transcripts we've currently added
//                   v.count = 0;
//                   markSeen(v.count);
//                }
//
//                // Note: We allow duplicate transcripts here --- they will always be adjacent
//                auto lastOffset = v.offset + unmarked(v.count);
//                transcriptIDs[lastOffset] = transcriptID;
//                v.count++;
//            }
//
//
//        }
//        if (transcriptID % 10000 == 0) {
//            std::cerr << "\r\rmarked kmers for " << transcriptID << " transcripts";
//        }
//        ++transcriptID;
//    }
//    	std::cerr << "\n";
//    }
//
//    //printVec(transcriptIDs);
//    // A hash to quickly and easily determine the equivalence classes
//    std::unordered_map<std::vector<uint32_t>, uint32_t, VectorHasher> eqClassMap;
//    // Holds the members of each equivalence class in the order in which classes
//    // are assigned.  The final size should be \sum_{c \in eqclasses} |c|.
//    std::vector<uint32_t> eqClassLabelVec;
//    // Holds pointer information about a k-mer's equivalence class.
//    // Specifically, where the label for the eq can be found
//    // as an offset and length into eqClassTxpVec, and where the
//    std::vector<rapmap::utils::EqClass> eqClasses;
//
//    uint32_t eqClassVecSize{0};
//
//    std::cerr << "\n[Step 3 of 4] : building k-mers equivalence classes\n";
//    // Compute the equivalence classes for the k-mers
//    {
//        ScopedTimer timer;
//        const auto ary = merIntMap.ary();
//        auto hashIt = ary->iterator_all<eager_iterator>();
//        std::vector<uint32_t> tlist;
//        while (hashIt.next()) {
//            auto& key = hashIt.key();
//            auto& val = kmerInfos[hashIt.val()];
//            //auto& val = kv.second;
//            auto offset = val.offset;
//            auto num = unmarked(val.count);
//
//            tlist.clear();
//            tlist.reserve(num);
//
//            for (size_t idx = offset; idx < offset + num; ++idx) {
//                auto tid = transcriptIDs[idx];
//                // We won't consider duplicate transcript IDs when building the
//                // equivalence classes
//		//
//		if (tlist.size() > 0 and tlist.back() > tid) {
//		    std::cerr << "Non monotnoically increasing transcript id!\n";
//		}
//                if (tlist.size() == 0 or tlist.back() != tid) {
//                    tlist.push_back(tid);
//                }
//            }
//
//            auto eqIt = eqClassMap.find(tlist);
//            uint32_t eqId{0};
//            // If there is no such equivalence class yet, then add it
//            if (eqIt == eqClassMap.end()) {
//                eqId = eqClassMap.size();
//                eqClassMap[tlist] = eqId;
//                // The label of this eq-class starts at eqClassVecSize and is
//                // tlist.size() elements long.
//                eqClasses.emplace_back(eqClassVecSize, tlist.size());
//                // Insert the label information into eqClassTxpVec
//                eqClassLabelVec.insert(eqClassLabelVec.end(), tlist.begin(), tlist.end());
//                eqClassVecSize += tlist.size();
//            } else {
//                eqId = eqIt->second;
//            }
//            // Set the equivalence class ID here for this transcript
//            val.eqId = eqId;
//            val.count = 0;
//        }
//    	std::cerr << "done! There were " << eqClassMap.size() << " classes\n";
//    }
//    // reuse the transcript IDs vector
//    auto& posVec = transcriptIDs;
//
//
//    constexpr uint8_t maxJump = std::numeric_limits<uint8_t>::max();
//    // Also, attempt to build *jump* tables here!
//    // How far we can move "forward" before hitting a new eq. class
//    std::vector<uint8_t> fwdJump(numDistinctKmers, maxJump);
//    // How far we can move "forward" backward hitting a new eq. class
//    std::vector<uint8_t> revJump(numDistinctKmers, maxJump);
//    int32_t lastBreak{0};
//
//    std::cerr << "\n[Step 4 of 4] : finalizing index\n";
//    transcriptID = 0;
//    {
//        ScopedTimer finalizeTimer;
//        // Local vector to hold k-mers per transcript
//        btree::btree_map<rapmap::utils::my_mer,
//                         std::vector<PosInfo>> posHash;//std::vector<PosInfo> posInfos;
//
//        // k-mers in the forward orientation w.r.t the reference txp
//        std::vector<JumpCell> fwdJumpQueue;
//        // k-mers in the reverse orientation w.r.t the reference txp
//        std::vector<JumpCell> revJumpQueue;
//
//        for (auto& transcriptSeq : transcriptSeqs) {
//	    // We can always jump to the beginning of a
//	    // new transcript
//	    lastBreak = 0;
//
//            fwdJumpQueue.clear();
//            revJumpQueue.clear();
//            posHash.clear();
//            auto readLen = transcriptSeq.length();
//            uint32_t currEqClass;
//            uint32_t prevEqClass = std::numeric_limits<uint32_t>::max();
//            rapmap::utils::my_mer mer;
//            mer.polyT();
//            for (size_t b = 0; b < readLen; ++b) {
//                int c = jellyfish::mer_dna::code(transcriptSeq[b]);
//                mer.shift_left(c);
//                if (b >= k) {
//                    auto canonicalMer = mer.get_canonical();
//                    bool isRC = (mer != canonicalMer);
//
//                    uint64_t kmerIndex;
//                    auto found = merIntMap.ary()->get_val_for_key(
//                                            canonicalMer, &kmerIndex);
//
//                    auto& val = kmerInfos[kmerIndex];
//                    currEqClass = val.eqId;
//
//                    // Record the position of this k-mer in the transcript
//                    uint32_t pos = b - k;
//                    if (pos > readLen) {
//                        std::cerr << "Pos is " << pos << ", but transcript length is " << readLen << "\n";
//                    }
//
//                    // === Jumping
//                    // if we hit a node with in-degree > 1 or out-degree > 1
//                    // then this defines the new breakpoint.  At this time, we
//                    // clear out the queues and mark the appropriate skips for
//                    // each k-mer we encountered.
//                    if ( currEqClass != prevEqClass ) {
//                        // For each k-mer in the forward direction, it can
//                        // skip forward to this breakpoint, which is at
//                        // position pos.
//                        emptyJumpQueue(fwdJumpQueue, lastBreak, pos, merIntMap,
//                                fwdJump, revJump);
//                        // The only difference here is that we reverse the
//                        // revJump and fwdJump arguments, since these are RC
//                        // mers.
//                        emptyJumpQueue(revJumpQueue, lastBreak, pos, merIntMap,
//                                revJump, fwdJump);
//                        lastBreak = pos;
//                    }
//                    // Does this k-mer exists in the table in the forward
//                    // or reverse complement direction.
//                    if (isRC) {
//                        revJumpQueue.emplace_back(kmerIndex, pos, isRC);
//                    } else {
//                        fwdJumpQueue.emplace_back(kmerIndex, pos, isRC);
//                    }
//                    prevEqClass = currEqClass;
//                    // === Jumping
//
//                    //posInfos.emplace_back(canonicalMer, isRC, pos);
//                    if (posHash[canonicalMer].size() > 0) {
//                        if (pos < posHash[canonicalMer].back().pos) {
//                            std::cerr << "NON-MONOTONIC POS\n";
//                        }
//                    }
//                    posHash[canonicalMer].emplace_back(kmerIndex, isRC, pos);
//                }
//            }
//
//            // === Jumping
//            // Empty anything remaining out of the jump queues
//            //
//            // The last k-mer in the transcript is, by definition a breakpoint.
//            // So, empty the queues.
//            emptyJumpQueue(fwdJumpQueue, lastBreak, pos, merIntMap,
//                    fwdJump, revJump);
//            // The only difference here is that we reverse the
//            // revJump and fwdJump arguments, since these are RC
//            // mers.
//            emptyJumpQueue(revJumpQueue, lastBreak, pos, merIntMap,
//                    revJump, fwdJump);
//            // === Jumping
//
//            for (auto kv = posHash.begin(); kv != posHash.end(); ++kv) {
//                    auto mer = kv->first;
//                    auto& list = kv->second;
//                    // Should ALWAYS find the key
//                    assert(found);
//                    auto& val = kmerInfos[list.front().merID];
//                    uint32_t offset;
//                    markNewTxpBit(list.front().pos);
//                    for (auto& pi : list) {
//                        offset = val.offset + val.count;
//                        if (pi.isRC) {
//                            markRCBit(pi.pos);
//                        }
//                        transcriptIDs[offset] = pi.pos;
//                        ++val.count;
//                    }
//                }
//            /*
//            std::sort(posInfos.begin(), posInfos.end(),
//                      [] (const PosInfo& a, const PosInfo& b) -> bool {
//                        if (a.mer < b.mer) {
//                            return true;
//                        } else if (a.mer == b.mer) {
//                            return a.pos < b.pos;
//                        } else {
//                            return false;
//                        }
//                      });
//            PosInfo* prev{nullptr};
//            for (auto& posInfo : posInfos) {
//                    uint64_t kmerIndex;
//                    auto found = merIntMap.ary()->get_val_for_key(
//                                    posInfo.mer, &kmerIndex);
//                    // Should ALWAYS find the key
//                    assert(found);
//                    auto& val = kmerInfos[kmerIndex];
//                    auto offset = val.offset + val.count;
//                    // Check if this offset is for a new transcript in the
//                    // position list and, if so, set the proper flag bit.
//                    if ( prev == nullptr or prev->mer != posInfo.mer) {
//                        markNewTxpBit(posInfo.pos);
//                    } else {
//                        // These are the same k-mer so the pos better
//                        // be increasing!
//                        if (prev != nullptr) {
//                            if ((prev->pos & uint30LowBitMask) >= posInfo.pos) {
//                                std::cerr << "prev pos = " << (prev->pos & uint30LowBitMask)
//                                          << ", but curr pos = " << posInfo.pos
//                                          << "\n";
//                            }
//                        }
//                    }
//                    if ( posInfo.isRC ) {
//                        markRCBit(posInfo.pos);
//                    }
//        		    prev = &posInfo;
//                    transcriptIDs[offset] = posInfo.pos;
//                    ++val.count;
//                }
//            */
//            if (transcriptID % 10000 == 0) {
//                std::cerr << "\r\rfinalized kmers for " << transcriptID << " transcripts";
//            }
//            ++transcriptID;
//        }
//	std::cerr << "\n";
//    }
//
//    /*
//    std::cerr << "[DEBUG]: Verifying Index\n";
//    {
//        ScopedTimer st;
//    auto hashIt = merIntMap.ary()->iterator_all<eager_iterator>();
//    std::vector<uint32_t> tlist;
//    while (hashIt.next()) {
//        auto& key = hashIt.key();
//        auto& val = kmerInfos[hashIt.val()];
//        if (!(*(transcriptIDs.begin() + val.offset) & uint32HighestBitMask)) {
//            std::cerr << "found un-marked k-mer position at beginning of list!\n";
//        }
//    }
//    }
//    */
//    // merIntMap --- jf hash
//    // kmerInfos --- info for each k-mer
//    // std::vector<uint32_t> eqClassLabelVec --- transcripts for each eq class
//    // std::vector<EqClass> eqClasses --- where each eq class starts and ends
//    // transcriptIDs --- position vec
//
//    std::cerr << "Writing the index to " << outputDir << "\n";
//
//    using JFFileHeader = jellyfish::file_header;
//    using JFDumper = jellyfish::binary_dumper<MerMapT::array>;
//
//    SpecialHeader fh;
//    fh.update_from_ary(*merIntMap.ary());
//    fh.canonical(true);
//    fh.format("gus/special"); // Thanks, Guillaume
//    fh.counter_len(8*sizeof(uint32_t)); // size of counter in bits
//    fh.fill_standard();
//    //fh.set_cmdline(argc, argv);
//
//    std::ofstream jfos(outputDir + "rapidx.jfhash");
//    fh.write(jfos);
//    merIntMap.ary()->write(jfos);
//    jfos.close();
//
//    // === Dump the jump tables to disk and reclaim the space
//    std::ofstream fwdJumpStream(outputDir + "fwdjump.bin", std::ios::binary);
//    {
//        cereal::BinaryOutputArchive fwdJumpArchive(fwdJumpStream);
//        fwdJumpArchive(fwdJump);
//    }
//    fwdJumpStream.close();
//    fwdJump.clear();
//    fwdJump.shrink_to_fit();
//
//    std::ofstream revJumpStream(outputDir + "revjump.bin", std::ios::binary);
//    {
//        cereal::BinaryOutputArchive revJumpArchive(revJumpStream);
//        revJumpArchive(revJump);
//    }
//    revJumpStream.close();
//    revJump.clear();
//    revJump.shrink_to_fit();
//    // === Done dumping the jump tables
//
//
//
//
//
//    std::ofstream kinfoStream(outputDir + "kinfo.bin", std::ios::binary);
//    {
//        cereal::BinaryOutputArchive kinfoArchive(kinfoStream);
//        kinfoArchive(kmerInfos);
//    }
//    kinfoStream.close();
//
//  std::ofstream eqLabelStream(outputDir + "eqlab.bin", std::ios::binary);
//  {
//      cereal::BinaryOutputArchive eqLabelArchive(eqLabelStream);
//      eqLabelArchive(eqClassLabelVec);
//  }
//  eqLabelStream.close();
//
//  std::ofstream eqClassStream(outputDir + "eqclass.bin", std::ios::binary);
//  {
//      cereal::BinaryOutputArchive eqClassArchive(eqClassStream);
//      eqClassArchive(eqClasses);
//  }
//  eqClassStream.close();
//
//  std::ofstream posStream(outputDir + "pos.bin", std::ios::binary);
//  {
//      cereal::BinaryOutputArchive posArchive(posStream);
//      posArchive(transcriptIDs);
//  }
//  posStream.close();
//  std::ofstream txpStream(outputDir + "txpnames.bin", std::ios::binary);
//  {
//    cereal::BinaryOutputArchive txpArchive(txpStream);
//    txpArchive(transcriptNames);
//  }
//  txpStream.close();
//
//  std::cerr << "transcriptIDs.size() = " << transcriptIDs.size() << "\n";
//  std::cerr << "parsed " << transcriptNames.size() << " transcripts\n";
//  std::cerr << "There were " << numDistinctKmers << " distinct k-mers (canonicalized)\n";
//
  // Data structure idea:
  // k-mer => equivalence class, position array offset
  // equivalence class = ordered (unique) list of transcripts
  // position array = *self-delimited* list of positions with the same order as txps in the equivalence class
}



int rapMapSAIndex(int argc, char* argv[]) {
    std::cerr << "RapMap Indexer\n";

    TCLAP::CmdLine cmd("RapMap Indexer");
    TCLAP::ValueArg<std::string> transcripts("t", "transcripts", "The transcript file to be indexed", true, "", "path");
    TCLAP::ValueArg<std::string> index("i", "index", "The location where the index should be written", true, "", "path");
    TCLAP::ValueArg<uint32_t> kval("k", "klen", "The length of k-mer to index", false, 31, "positive integer less than 32");
    cmd.add(transcripts);
    cmd.add(index);
    cmd.add(kval);

    cmd.parse(argc, argv);

    // stupid parsing for now
    std::string transcriptFile(transcripts.getValue());
    std::vector<std::string> transcriptFiles({ transcriptFile });

    uint32_t k = kval.getValue();
    rapmap::utils::my_mer::k(k);

    std::string indexDir = index.getValue();
    if (indexDir.back() != '/') {
	indexDir += '/';
    }
    bool dirExists = rapmap::fs::DirExists(indexDir.c_str());
    bool dirIsFile = rapmap::fs::FileExists(indexDir.c_str());
    if (dirIsFile) {
        std::cerr << "The requested index directory already exists as a file.";
        std::exit(1);
    }
    if (!dirExists) {
        rapmap::fs::MakeDir(indexDir.c_str());
    }

    size_t maxReadGroup{1000}; // Number of reads in each "job"
    size_t concurrentFile{2}; // Number of files to read simultaneously
    size_t numThreads{2};
    stream_manager streams(transcriptFiles.begin(), transcriptFiles.end(), concurrentFile);
    std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
    transcriptParserPtr.reset(new single_parser(4 * numThreads, maxReadGroup,
                              concurrentFile, streams));
    std::mutex iomutex;
    indexTranscriptsSA(transcriptParserPtr.get(), indexDir, iomutex);
    return 0;
}

