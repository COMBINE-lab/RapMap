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
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <type_traits>
#include <fstream>

#include "tclap/CmdLine.h"

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>

#include "xxhash.h"
#include "btree/btree_map.h"

#include "FastxParser.hpp"
// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"

#include "RapMapUtils.hpp"
#include "RapMapFileSystem.hpp"
#include "ScopedTimer.hpp"
#include "IndexHeader.hpp"

#include "jellyfish/file_header.hpp"
#include "jellyfish/binary_dumper.hpp"
#include "jellyfish/thread_exec.hpp"
#include "jellyfish/hash_counter.hpp"
#include "jellyfish/mer_overlap_sequence_parser.hpp"
#include "jellyfish/mer_iterator.hpp"
#include "JFRaw.hpp"

#include <chrono>

using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;
using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using MerMapT = jellyfish::cooperative::hash_counter<rapmap::utils::my_mer>;

uint64_t encode(uint32_t tid, uint32_t pos) {
    uint64_t res = tid;
    res = res << 32;
    res |= pos;
    return res;
}

constexpr uint32_t uint32HighestBitMask = 0x80000000;
constexpr uint32_t uint31LowBitMask = 0x7FFFFFFF;

constexpr uint32_t uint30LowBitMask = 0x3FFFFFFF;

constexpr uint32_t rcSetMask = 0x40000000;

// marks the second highest bit
void markRCBit(uint32_t& i) { i |= rcSetMask; }

// marks the highest bit
void markNewTxpBit(uint32_t& i) { i |= uint32HighestBitMask; }

bool wasSeen(uint32_t i) { return ((i & uint32HighestBitMask) >> 31) == 1; }

void markSeen(uint32_t& i) { i |= uint32HighestBitMask; }

uint32_t unmarked(uint32_t i) { return (i & uint31LowBitMask); }

class VectorHasher {
    public:
    size_t operator()(const std::vector<uint32_t>& vec) const {
        size_t hashVal{0};
        for (auto v : vec) {
            rapmap::utils::hashCombine(hashVal, v);
        }
        return hashVal;
    }
};

struct PosInfo {
    PosInfo(uint64_t merIDIn, bool isRCIn, uint32_t posIn) :
        merID(merIDIn), isRC(isRCIn), pos(posIn) {}

    uint64_t merID;
    bool isRC;
    uint32_t pos;
};

// maybe unify with the above?
struct JumpCell {
    JumpCell(uint32_t merIdxIn, int32_t posIn, bool isRCIn) :
        merIdx(merIdxIn), pos(posIn), isRC(isRCIn) {}
    uint32_t merIdx;
    int32_t pos;
    bool isRC;
};

template <typename T>
void printVec(std::vector<T>& vec) {
    std::cerr << "{ ";
    for (auto& e : vec) {
        std::cerr << e << ", ";
    }
    std::cerr << "}";
}

// There may be a better way to do this, but here we just check the
// possible neighbors
bool isBreakpoint(MerMapT& merIntMap, rapmap::utils::my_mer canonicalMer) {
    uint32_t inDegree{0};
    uint32_t outDegree{0};
    const auto& ary = merIntMap.ary();
    // extend forward
    for (char b : {'A', 'C', 'G', 'T'}) {
        auto newMer = canonicalMer;
        newMer.shift_left(b);
        newMer.canonicalize();
        outDegree += ary->has_key(newMer) ? 1 : 0;
        if (outDegree > 1) { return true; }
    }
    // extend backward
    for (char b : {'A', 'C', 'G', 'T'}) {
        auto newMer = canonicalMer;
        newMer.shift_right(b);
        newMer.canonicalize();
        inDegree += ary->has_key(newMer) ? 1 : 0;
        if (inDegree > 1) { return true; }
    }
    return false;
}

void emptyJumpQueue(std::vector<JumpCell>& jumpQueue, int32_t lastBreak,
                    int32_t pos,
                    MerMapT& merIntMap,
                    std::vector<uint8_t>& fwdJump,
                    std::vector<uint8_t>& revJump) {

    // The maximum representable jump
    constexpr auto maxJump =
        std::numeric_limits<std::remove_reference<decltype(fwdJump[0])>::type>::max();

    while (!jumpQueue.empty()) {
        auto& jumpCell = jumpQueue.back();
        uint8_t revJumpDist = static_cast<uint8_t>(
                std::min(jumpCell.pos - lastBreak + 1,
                static_cast<int32_t>(maxJump)));
        uint8_t fwdJumpDist = static_cast<uint8_t>(
                std::min(pos - jumpCell.pos + 1,
                static_cast<int32_t>(maxJump)));
        fwdJump[jumpCell.merIdx] =
            std::min(fwdJumpDist, fwdJump[jumpCell.merIdx]);
        revJump[jumpCell.merIdx] =
            std::min(revJumpDist, revJump[jumpCell.merIdx]);
        // Now that we marked the jump, no need to keep this around.
        jumpQueue.pop_back();
    }
}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT>//, typename CoverageCalculator>
void processTranscripts(ParserT* parser,
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
    std::vector<uint32_t> transcriptLengths;
    constexpr char bases[] = {'A', 'C', 'G', 'T'};
    uint32_t polyAClipLength{10};
    uint32_t numPolyAsClipped{0};
    std::string polyA(polyAClipLength, 'A');

    using TranscriptList = std::vector<uint32_t>;
    using eager_iterator = MerMapT::array::eager_iterator;
    using KmerBinT = uint64_t;
    //create the hash
    size_t hashSize = 100000000;
    MerMapT merIntMap(hashSize, rapmap::utils::my_mer::k()*2, 32, 1, 126);
    std::vector<rapmap::utils::KmerInfo> kmerInfos;

    std::vector<std::string> transcriptSeqs;
    size_t numDistinctKmers{0};
    size_t numKmers{0};

    std::cerr << "\n[Step 1 of 4] : counting k-mers\n";

    {
        ScopedTimer timer;
	// Get the read group by which this thread will
	// communicate with the parser (*once per-thread*)
	auto rg = parser->getReadGroup();

	while (parser->refill(rg)) {
	  for (auto& read : rg) { // for each sequence
	    std::string& readStr = read.seq; 

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
                transcriptLengths.push_back(readLen);
                auto& recHeader = read.name;
                transcriptNames.emplace_back(recHeader.substr(0, recHeader.find_first_of(" \t")));

                rapmap::utils::my_mer mer;
                mer.polyT();
                for (size_t b = 0; b < readLen; ++b) {
                    int c = jellyfish::mer_dna::code(readStr[b]);
                    if (jellyfish::mer_dna::not_dna(c)) {
                        char rbase = bases[dis(eng)];
                        c = jellyfish::mer_dna::code(rbase);
                        readStr[b] = rbase;
                    }
                    mer.shift_left(c);
                    if (b >= k) {
                        auto canonicalMer = mer.get_canonical();
                        auto key = canonicalMer.get_bits(0, 2*k);

                        uint64_t val;
                        auto found = merIntMap.ary()->get_val_for_key(canonicalMer, &val);
                        if (!found) {
                            merIntMap.add(canonicalMer, numDistinctKmers);
                            kmerInfos.emplace_back(txpIndex, 0, 1);
                            ++numDistinctKmers;
                        } else {
                            kmerInfos[val].count++;
                        }

                        /*
                           start = std::chrono::high_resolution_clock::now();
                           auto it = cmap.find(key);
                           end = std::chrono::high_resolution_clock::now();
                           elapsedNS += end - start;
                        // If we found the k-mer, increment the count
                        if (it != cmap.end()) {
                        it->second.count++;
                        } else { // Otherwise, add it
                        cmap[key] = KmerInfo(txpIndex, 0, 1);
                        }
                        */
                        // No matter what, our k-mer count increased
                        numKmers++;
                    }
                }
                transcriptSeqs.push_back(read.seq);
                if (n % 10000 == 0) {
                    std::cerr << "\r\rcounted k-mers for " << n << " transcripts";
                }
            }
        }
        std::cerr << "\n";
    }
    parser->stop();
    std::cerr << "Clipped poly-A tails from " << numPolyAsClipped << " transcripts\n";

    std::ofstream txpLenStream(outputDir + "txplens.bin", std::ios::binary);
    {
        cereal::BinaryOutputArchive txpLenArchive(txpLenStream);
        txpLenArchive(transcriptLengths);
    }
    txpLenStream.close();
    transcriptLengths.clear();
    transcriptLengths.shrink_to_fit();

    constexpr uint32_t uint32Invalid = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> transcriptIDs(numKmers, uint32Invalid);

    std::cerr << "\n[Step 2 of 4] : marking k-mers\n";
    // Mark the transcript in which each occurence oc a k-mer appears
    // in the transcriptIDs vector.


    bool isRC{false};
    int32_t pos{0};
    uint32_t offset{0};
    uint32_t transcriptID{0};
    {
    ScopedTimer timer;


    for (auto& transcriptSeq : transcriptSeqs) {
        auto readLen = transcriptSeq.length();
        rapmap::utils::my_mer mer;
        mer.polyT();
        std::vector<KmerBinT> kmers;
        for (size_t b = 0; b < readLen; ++b) {
            int c = jellyfish::mer_dna::code(transcriptSeq[b]);
            mer.shift_left(c);
            if (b >= k) {
                auto canonicalMer = mer.get_canonical();
                uint64_t kmerIndex;
                auto found = merIntMap.ary()->get_val_for_key(canonicalMer, &kmerIndex);
                // Should ALWAYS find the key
                assert(found);


                auto& v = kmerInfos[kmerIndex];
                // use the highest bit to mark if we've seen this k-mer yet or not
                // If we haven't seen this k-mer yet
                if (!wasSeen(v.count)) {
                   // Where we start looking for transcripts for this k-mer
                   v.offset = offset;
                   offset += v.count;
                   // The number of transcripts we've currently added
                   v.count = 0;
                   markSeen(v.count);
                }

                // Note: We allow duplicate transcripts here --- they will always be adjacent
                auto lastOffset = v.offset + unmarked(v.count);
                transcriptIDs[lastOffset] = transcriptID;
                v.count++;
            }


        }
        if (transcriptID % 10000 == 0) {
            std::cerr << "\r\rmarked kmers for " << transcriptID << " transcripts";
        }
        ++transcriptID;
    }
    	std::cerr << "\n";
    }

    //printVec(transcriptIDs);
    // A hash to quickly and easily determine the equivalence classes
    std::unordered_map<std::vector<uint32_t>, uint32_t, VectorHasher> eqClassMap;
    // Holds the members of each equivalence class in the order in which classes
    // are assigned.  The final size should be \sum_{c \in eqclasses} |c|.
    std::vector<uint32_t> eqClassLabelVec;
    // Holds pointer information about a k-mer's equivalence class.
    // Specifically, where the label for the eq can be found
    // as an offset and length into eqClassTxpVec, and where the
    std::vector<rapmap::utils::EqClass> eqClasses;

    uint32_t eqClassVecSize{0};

    std::cerr << "\n[Step 3 of 4] : building k-mers equivalence classes\n";
    // Compute the equivalence classes for the k-mers
    {
        ScopedTimer timer;
        const auto ary = merIntMap.ary();
        auto hashIt = ary->iterator_all<eager_iterator>();
        std::vector<uint32_t> tlist;
        while (hashIt.next()) {
            auto& key = hashIt.key();
            auto& val = kmerInfos[hashIt.val()];
            //auto& val = kv.second;
            auto offset = val.offset;
            auto num = unmarked(val.count);

            tlist.clear();
            tlist.reserve(num);

            for (size_t idx = offset; idx < offset + num; ++idx) {
                auto tid = transcriptIDs[idx];
                // We won't consider duplicate transcript IDs when building the
                // equivalence classes
		//
		if (tlist.size() > 0 and tlist.back() > tid) {
		    std::cerr << "Non monotnoically increasing transcript id!\n";
		}
                if (tlist.size() == 0 or tlist.back() != tid) {
                    tlist.push_back(tid);
                }
            }

            auto eqIt = eqClassMap.find(tlist);
            uint32_t eqId{0};
            // If there is no such equivalence class yet, then add it
            if (eqIt == eqClassMap.end()) {
                eqId = eqClassMap.size();
                eqClassMap[tlist] = eqId;
                // The label of this eq-class starts at eqClassVecSize and is
                // tlist.size() elements long.
                eqClasses.emplace_back(eqClassVecSize, tlist.size());
                // Insert the label information into eqClassTxpVec
                eqClassLabelVec.insert(eqClassLabelVec.end(), tlist.begin(), tlist.end());
                eqClassVecSize += tlist.size();
            } else {
                eqId = eqIt->second;
            }
            // Set the equivalence class ID here for this transcript
            val.eqId = eqId;
            val.count = 0;
        }
    	std::cerr << "done! There were " << eqClassMap.size() << " classes\n";
    }
    // reuse the transcript IDs vector
    auto& posVec = transcriptIDs;


    constexpr uint8_t maxJump = std::numeric_limits<uint8_t>::max();
    // Also, attempt to build *jump* tables here!
    // How far we can move "forward" before hitting a new eq. class
    std::vector<uint8_t> fwdJump(numDistinctKmers, maxJump);
    // How far we can move "forward" backward hitting a new eq. class
    std::vector<uint8_t> revJump(numDistinctKmers, maxJump);
    int32_t lastBreak{0};

    std::cerr << "\n[Step 4 of 4] : finalizing index\n";
    transcriptID = 0;
    {
        ScopedTimer finalizeTimer;
        // Local vector to hold k-mers per transcript
        btree::btree_map<rapmap::utils::my_mer,
                         std::vector<PosInfo>> posHash;//std::vector<PosInfo> posInfos;

        // k-mers in the forward orientation w.r.t the reference txp
        std::vector<JumpCell> fwdJumpQueue;
        // k-mers in the reverse orientation w.r.t the reference txp
        std::vector<JumpCell> revJumpQueue;

        for (auto& transcriptSeq : transcriptSeqs) {
	    // We can always jump to the beginning of a
	    // new transcript
	    lastBreak = 0;

            fwdJumpQueue.clear();
            revJumpQueue.clear();
            posHash.clear();
            auto readLen = transcriptSeq.length();
            uint32_t currEqClass;
            uint32_t prevEqClass = std::numeric_limits<uint32_t>::max();
            rapmap::utils::my_mer mer;
            mer.polyT();
            for (size_t b = 0; b < readLen; ++b) {
                int c = jellyfish::mer_dna::code(transcriptSeq[b]);
                mer.shift_left(c);
                if (b >= k) {
                    auto canonicalMer = mer.get_canonical();
                    bool isRC = (mer != canonicalMer);

                    uint64_t kmerIndex;
                    auto found = merIntMap.ary()->get_val_for_key(
                                            canonicalMer, &kmerIndex);

                    auto& val = kmerInfos[kmerIndex];
                    currEqClass = val.eqId;

                    // Record the position of this k-mer in the transcript
                    uint32_t pos = b - k;
                    if (pos > readLen) {
                        std::cerr << "Pos is " << pos << ", but transcript length is " << readLen << "\n";
                    }

                    // === Jumping
                    // if we hit a node with in-degree > 1 or out-degree > 1
                    // then this defines the new breakpoint.  At this time, we
                    // clear out the queues and mark the appropriate skips for
                    // each k-mer we encountered.
                    if ( currEqClass != prevEqClass ) {
                        // For each k-mer in the forward direction, it can
                        // skip forward to this breakpoint, which is at
                        // position pos.
                        emptyJumpQueue(fwdJumpQueue, lastBreak, pos, merIntMap,
                                fwdJump, revJump);
                        // The only difference here is that we reverse the
                        // revJump and fwdJump arguments, since these are RC
                        // mers.
                        emptyJumpQueue(revJumpQueue, lastBreak, pos, merIntMap,
                                revJump, fwdJump);
                        lastBreak = pos;
                    }
                    // Does this k-mer exists in the table in the forward
                    // or reverse complement direction.
                    if (isRC) {
                        revJumpQueue.emplace_back(kmerIndex, pos, isRC);
                    } else {
                        fwdJumpQueue.emplace_back(kmerIndex, pos, isRC);
                    }
                    prevEqClass = currEqClass;
                    // === Jumping

                    //posInfos.emplace_back(canonicalMer, isRC, pos);
                    if (posHash[canonicalMer].size() > 0) {
                        if (pos < posHash[canonicalMer].back().pos) {
                            std::cerr << "NON-MONOTONIC POS\n";
                        }
                    }
                    posHash[canonicalMer].emplace_back(kmerIndex, isRC, pos);
                }
            }

            // === Jumping
            // Empty anything remaining out of the jump queues
            //
            // The last k-mer in the transcript is, by definition a breakpoint.
            // So, empty the queues.
            emptyJumpQueue(fwdJumpQueue, lastBreak, pos, merIntMap,
                    fwdJump, revJump);
            // The only difference here is that we reverse the
            // revJump and fwdJump arguments, since these are RC
            // mers.
            emptyJumpQueue(revJumpQueue, lastBreak, pos, merIntMap,
                    revJump, fwdJump);
            // === Jumping

            for (auto kv = posHash.begin(); kv != posHash.end(); ++kv) {
                    auto mer = kv->first;
                    auto& list = kv->second;
                    // Should ALWAYS find the key
                    assert(found);
                    auto& val = kmerInfos[list.front().merID];
                    uint32_t offset;
                    markNewTxpBit(list.front().pos);
                    for (auto& pi : list) {
                        offset = val.offset + val.count;
                        if (pi.isRC) {
                            markRCBit(pi.pos);
                        }
                        transcriptIDs[offset] = pi.pos;
                        ++val.count;
                    }
                }
            /*
            std::sort(posInfos.begin(), posInfos.end(),
                      [] (const PosInfo& a, const PosInfo& b) -> bool {
                        if (a.mer < b.mer) {
                            return true;
                        } else if (a.mer == b.mer) {
                            return a.pos < b.pos;
                        } else {
                            return false;
                        }
                      });
            PosInfo* prev{nullptr};
            for (auto& posInfo : posInfos) {
                    uint64_t kmerIndex;
                    auto found = merIntMap.ary()->get_val_for_key(
                                    posInfo.mer, &kmerIndex);
                    // Should ALWAYS find the key
                    assert(found);
                    auto& val = kmerInfos[kmerIndex];
                    auto offset = val.offset + val.count;
                    // Check if this offset is for a new transcript in the
                    // position list and, if so, set the proper flag bit.
                    if ( prev == nullptr or prev->mer != posInfo.mer) {
                        markNewTxpBit(posInfo.pos);
                    } else {
                        // These are the same k-mer so the pos better
                        // be increasing!
                        if (prev != nullptr) {
                            if ((prev->pos & uint30LowBitMask) >= posInfo.pos) {
                                std::cerr << "prev pos = " << (prev->pos & uint30LowBitMask)
                                          << ", but curr pos = " << posInfo.pos
                                          << "\n";
                            }
                        }
                    }
                    if ( posInfo.isRC ) {
                        markRCBit(posInfo.pos);
                    }
        		    prev = &posInfo;
                    transcriptIDs[offset] = posInfo.pos;
                    ++val.count;
                }
            */
            if (transcriptID % 10000 == 0) {
                std::cerr << "\r\rfinalized kmers for " << transcriptID << " transcripts";
            }
            ++transcriptID;
        }
	std::cerr << "\n";
    }

    /*
    std::cerr << "[DEBUG]: Verifying Index\n";
    {
        ScopedTimer st;
    auto hashIt = merIntMap.ary()->iterator_all<eager_iterator>();
    std::vector<uint32_t> tlist;
    while (hashIt.next()) {
        auto& key = hashIt.key();
        auto& val = kmerInfos[hashIt.val()];
        if (!(*(transcriptIDs.begin() + val.offset) & uint32HighestBitMask)) {
            std::cerr << "found un-marked k-mer position at beginning of list!\n";
        }
    }
    }
    */
    // merIntMap --- jf hash
    // kmerInfos --- info for each k-mer
    // std::vector<uint32_t> eqClassLabelVec --- transcripts for each eq class
    // std::vector<EqClass> eqClasses --- where each eq class starts and ends
    // transcriptIDs --- position vec

    std::cerr << "Writing the index to " << outputDir << "\n";

    using JFFileHeader = jellyfish::file_header;
    using JFDumper = jellyfish::binary_dumper<MerMapT::array>;

    SpecialHeader fh;
    fh.update_from_ary(*merIntMap.ary());
    fh.canonical(true);
    fh.format("gus/special"); // Thanks, Guillaume
    fh.counter_len(8*sizeof(uint32_t)); // size of counter in bits
    fh.fill_standard();
    //fh.set_cmdline(argc, argv);

    std::ofstream jfos(outputDir + "rapidx.jfhash");
    fh.write(jfos);
    merIntMap.ary()->write(jfos);
    jfos.close();

    // === Dump the jump tables to disk and reclaim the space
    std::ofstream fwdJumpStream(outputDir + "fwdjump.bin", std::ios::binary);
    {
        cereal::BinaryOutputArchive fwdJumpArchive(fwdJumpStream);
        fwdJumpArchive(fwdJump);
    }
    fwdJumpStream.close();
    fwdJump.clear();
    fwdJump.shrink_to_fit();

    std::ofstream revJumpStream(outputDir + "revjump.bin", std::ios::binary);
    {
        cereal::BinaryOutputArchive revJumpArchive(revJumpStream);
        revJumpArchive(revJump);
    }
    revJumpStream.close();
    revJump.clear();
    revJump.shrink_to_fit();
    // === Done dumping the jump tables





    std::ofstream kinfoStream(outputDir + "kinfo.bin", std::ios::binary);
    {
        cereal::BinaryOutputArchive kinfoArchive(kinfoStream);
        kinfoArchive(kmerInfos);
    }
    kinfoStream.close();

  std::ofstream eqLabelStream(outputDir + "eqlab.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive eqLabelArchive(eqLabelStream);
      eqLabelArchive(eqClassLabelVec);
  }
  eqLabelStream.close();

  std::ofstream eqClassStream(outputDir + "eqclass.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive eqClassArchive(eqClassStream);
      eqClassArchive(eqClasses);
  }
  eqClassStream.close();

  std::ofstream posStream(outputDir + "pos.bin", std::ios::binary);
  {
      cereal::BinaryOutputArchive posArchive(posStream);
      posArchive(transcriptIDs);
  }
  posStream.close();
  std::ofstream txpStream(outputDir + "txpnames.bin", std::ios::binary);
  {
    cereal::BinaryOutputArchive txpArchive(txpStream);
    txpArchive(transcriptNames);
  }
  txpStream.close();

  std::string indexVersion = "p0";
  IndexHeader header(IndexType::PSEUDO, indexVersion, true, k);
  // Finally (since everything presumably succeeded) write the header
  std::ofstream headerStream(outputDir + "header.json");
  {
    cereal::JSONOutputArchive archive(headerStream);
    archive(header);
  }
  headerStream.close();


  std::cerr << "transcriptIDs.size() = " << transcriptIDs.size() << "\n";
  std::cerr << "parsed " << transcriptNames.size() << " transcripts\n";
  std::cerr << "There were " << numDistinctKmers << " distinct k-mers (canonicalized)\n";

  // Data structure idea:
  // k-mer => equivalence class, position array offset
  // equivalence class = ordered (unique) list of transcripts
  // position array = *self-delimited* list of positions with the same order as txps in the equivalence class
}



int rapMapIndex(int argc, char* argv[]) {
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
    if (k % 2 == 0) {
        std::cerr << "Error: k must be an odd value, you chose " << k << '\n';
        std::exit(1);
    } else if (k > 31) {
        std::cerr << "Error: k must not be larger than 31, you chose " << k << '\n';
        std::exit(1);
    }
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

    size_t numThreads{1};

    std::unique_ptr<single_parser> transcriptParserPtr{nullptr};
    //transcriptParserPtr.reset(
    //    new single_parser(4 * numThreads, maxReadGroup, concurrentFile, streams));

    size_t numProd = 1;
    transcriptParserPtr.reset(
			      new single_parser(transcriptFiles, numThreads, numProd));

    transcriptParserPtr->start();
    std::mutex iomutex;
    processTranscripts(transcriptParserPtr.get(), indexDir, iomutex);
    return 0;
}


