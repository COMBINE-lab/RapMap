#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/archives/binary.hpp>
#include "RapMapUtils.hpp"
#include "RapMapSAIndex.hpp"
#include "RapMapIndex.hpp"
#include "PairAlignmentFormatter.hpp"
#include "SingleAlignmentFormatter.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

namespace rapmap {
    namespace utils {
        std::vector<std::string> tokenize(const std::string &s, char delim) {
            std::stringstream ss(s);
            std::string item;
            std::vector<std::string> elems;
            while (std::getline(ss, item, delim)) {
                elems.push_back(item);
            }
            return elems;
        }


		// positions are stored in a packed format, where the highest
		// 2-bits encode if this position refers to a new transcript
		// and whether or not the k-mer from the hash matches this txp
		// in the forward or rc direction.
		void decodePosition(uint32_t p, int32_t& pOut, bool& newTxp, bool& isRC) {
			uint32_t highBits = (p >> 30);
			pOut = (p & 0x3fffffff);
			newTxp = (highBits & 0x1);
			isRC = (highBits & 0x2);
		}



        static constexpr int8_t rc_table[128] = {
            78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 15
            78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 31
            78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 787
            78, 78,  78, 78,  78,  78,  78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 63
            78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // 79
            78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78, // 95
            78, 84, 78, 71, 78,  78,  78, 67, 78, 78, 78, 78,  78, 78, 78, 78, // 101
            78, 78,  78, 78,  65, 65, 78, 78,  78, 78, 78, 78,  78, 78, 78, 78  // 127
        };

        // Adapted from
        // https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
        void reverseRead(std::string& seq,
                std::string& qual,
                std::string& readWork,
                std::string& qualWork) {

            readWork.resize(seq.length(), 'A');
            qualWork.resize(qual.length(), 'I');
            int32_t end = seq.length()-1, start = 0;
            //readWork[end] = '\0';
            //qualWork[end] = '\0';
            -- end;
            while (LIKELY(start < end)) {
                readWork[start] = (char)rc_table[(int8_t)seq[end]];
                readWork[end] = (char)rc_table[(int8_t)seq[start]];
                qualWork[start] = qual[end];
                qualWork[end] = qual[start];
                ++ start;
                -- end;
            }
            // If odd # of bases, we still have to complement the middle
            if (start == end) {
                readWork[start] = (char)rc_table[(int8_t)seq[start]];
                // but don't need to mess with quality
                // qualWork[start] = qual[start];
            }
            //std::swap(seq, readWork);
            //std::swap(qual, qualWork);
        }

        template <typename ReadT, typename IndexT>
        uint32_t writeAlignmentsToStream(
                ReadT& r,
                SingleAlignmentFormatter<IndexT>& formatter,
                HitCounters& hctr,
                std::vector<rapmap::utils::QuasiAlignment>& hits,
                fmt::MemoryWriter& sstream
                ) {
                // Convenient variable name bindings
                auto& txpNames = formatter.index->txpNames;
                auto& txpLens = formatter.index->txpLens;

                auto& readTemp = formatter.readTemp;
                auto& qualTemp = formatter.qualTemp;
                auto& cigarStr = formatter.cigarStr;

                uint16_t flags;

                auto& readName = r.header;
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                auto before = readName.find_first_of(':');
                before = readName.find_first_of(':', before+1);
                auto after = readName.find_first_of(':', before+1);
                const auto& txpName = readName.substr(before+1, after-before-1);
#endif //__DEBUG__
                // If the read name contains multiple space-separated parts, print
                // only the first
                size_t splitPos = readName.find(' ');
                if (splitPos < readName.length()) {
                    readName[splitPos] = '\0';
                }


                std::string numHitFlag = fmt::format("NH:i:{}", hits.size());
                uint32_t alnCtr{0};
                bool haveRev{false};
                for (auto& qa : hits) {
                    auto& transcriptName = txpNames[qa.tid];
                    // === SAM
                    rapmap::utils::getSamFlags(qa, flags);
                    if (alnCtr != 0) {
                        flags |= 0x900;
                    }

                    std::string* readSeq = &(r.seq);
                    std::string* qstr = &(r.qual);

                    if (!qa.fwd) {
                        if (!haveRev) {
                            rapmap::utils::reverseRead(*readSeq, *qstr,
                                                       readTemp, qualTemp);
                            haveRev = true;
                        }
                        readSeq = &(readTemp);
                        qstr = &(qualTemp);
                    }

                   rapmap::utils::adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], cigarStr);

                    sstream << readName.c_str() << '\t' // QNAME
                        << flags << '\t' // FLAGS
                        << transcriptName << '\t' // RNAME
                        << qa.pos + 1 << '\t' // POS (1-based)
                        << 255 << '\t' // MAPQ
                        << cigarStr.c_str() << '\t' // CIGAR
                        << '*' << '\t' // MATE NAME
                        << 0 << '\t' // MATE POS
                        << qa.fragLen << '\t' // TLEN
                        << *readSeq << '\t' // SEQ
                        << *qstr << '\t' // QSTR
                        << numHitFlag << '\n';
                    ++alnCtr;
                    // === SAM
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    if (txpNames[qa.tid] == txpName) { ++hctr.trueHits; }
#endif //__DEBUG__
                }
                return alnCtr;
            }

        // For reads paired *in sequencing*
        template <typename ReadPairT, typename IndexT>
        uint32_t writeAlignmentsToStream(
                ReadPairT& r,
                PairAlignmentFormatter<IndexT>& formatter,
                HitCounters& hctr,
                std::vector<rapmap::utils::QuasiAlignment>& jointHits,
                fmt::MemoryWriter& sstream
                ) {
                // Convenient variable name bindings
                auto& txpNames = formatter.index->txpNames;
                auto& txpLens = formatter.index->txpLens;

                auto& read1Temp = formatter.read1Temp;
                auto& read2Temp = formatter.read2Temp;
                auto& qual1Temp = formatter.qual1Temp;
                auto& qual2Temp = formatter.qual2Temp;
                auto& cigarStr1 = formatter.cigarStr1;
                auto& cigarStr2 = formatter.cigarStr2;

                uint16_t flags1, flags2;

                auto& readName = r.first.header;
                // If the read name contains multiple space-separated parts,
                // print only the first
                size_t splitPos = readName.find(' ');
                if (splitPos < readName.length()) {
                    readName[splitPos] = '\0';
                } else {
                    splitPos = readName.length();
                }

                // trim /1 from the pe read
                if (splitPos > 2 and readName[splitPos - 2] == '/') {
                    readName[splitPos - 2] = '\0';
                }

                auto& mateName = r.second.header;
                // If the read name contains multiple space-separated parts,
                // print only the first
                splitPos = mateName.find(' ');
                if (splitPos < mateName.length()) {
                    mateName[splitPos] = '\0';
                } else {
                    splitPos = mateName.length();
                }

                // trim /2 from the pe read
                if (splitPos > 2 and mateName[splitPos - 2] == '/') {
                    mateName[splitPos - 2] = '\0';
                }

                /*
                // trim /1 and /2 from pe read names
                if (readName.length() > 2 and
                        readName[readName.length() - 2] == '/') {
                    readName[readName.length() - 2] = '\0';
                }
                if (mateName.length() > 2 and
                        mateName[mateName.length() - 2] == '/') {
                    mateName[mateName.length() - 2] = '\0';
                }
                */

                std::string numHitFlag = fmt::format("NH:i:{}", jointHits.size());
                uint32_t alnCtr{0};
				uint32_t trueHitCtr{0};
				QuasiAlignment* firstTrueHit{nullptr};
                bool haveRev1{false};
                bool haveRev2{false};
                bool* haveRev = nullptr;
                for (auto& qa : jointHits) {
                    auto& transcriptName = txpNames[qa.tid];
                    // === SAM
                    if (qa.isPaired) {
                        rapmap::utils::getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            flags2 |= 0x900;
                        }
                        rapmap::utils::adjustOverhang(qa, txpLens[qa.tid], cigarStr1, cigarStr2);

                        // Reverse complement the read and reverse
                        // the quality string if we need to
                        std::string* readSeq1 = &(r.first.seq);
                        std::string* qstr1 = &(r.first.qual);
                        if (!qa.fwd) {
                            if (!haveRev1) {
                                rapmap::utils::reverseRead(*readSeq1, *qstr1,
                                        read1Temp, qual1Temp);
                                haveRev1 = true;
                            }
                            readSeq1 = &(read1Temp);
                            qstr1 = &(qual1Temp);
                        }

                        std::string* readSeq2 = &(r.second.seq);
                        std::string* qstr2 = &(r.second.qual);
                        if (!qa.mateIsFwd) {
                            if (!haveRev2) {
                                rapmap::utils::reverseRead(*readSeq2, *qstr2,
                                        read2Temp, qual2Temp);
                                haveRev2 = true;
                            }
                            readSeq2 = &(read2Temp);
                            qstr2 = &(qual2Temp);
                        }

                        sstream << readName.c_str() << '\t' // QNAME
                                << flags1 << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.pos + 1 << '\t' // POS (1-based)
                                << 1 << '\t' // MAPQ
                                << cigarStr1.c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << qa.matePos + 1 << '\t' // PNEXT
                                << qa.fragLen << '\t' // TLEN
                                << *readSeq1 << '\t' // SEQ
                                << *qstr1 << '\t' // QUAL
                                << numHitFlag << '\n';

                        sstream << mateName.c_str() << '\t' // QNAME
                                << flags2 << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.matePos + 1 << '\t' // POS (1-based)
                                << 1 << '\t' // MAPQ
                                << cigarStr2.c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << qa.pos + 1 << '\t' // PNEXT
                                << qa.fragLen << '\t' // TLEN
                                << *readSeq2 << '\t' // SEQ
                                << *qstr2 << '\t' // QUAL
                                << numHitFlag << '\n';
                    } else {
                        rapmap::utils::getSamFlags(qa, true, flags1, flags2);
                        if (alnCtr != 0) {
                            flags1 |= 0x900; flags2 |= 0x900;
                        } else {
                            // If this is the first alignment for this read
                            // If the left end is mapped, set 0x900 on the right end
                            if (qa.mateStatus == MateStatus::PAIRED_END_LEFT) {
                                flags2 |= 0x900;
                            } else {
                            // Otherwise, set 0x900 on the left end
                                flags1 |= 0x900;
                            }
                        }

                        std::string* readSeq{nullptr};
                        std::string* unalignedSeq{nullptr};

                        uint32_t flags, unalignedFlags;
                        std::string* qstr{nullptr};
                        std::string* unalignedQstr{nullptr};
                        std::string* alignedName{nullptr};
                        std::string* unalignedName{nullptr};
                        std::string* readTemp{nullptr};
                        std::string* qualTemp{nullptr};

                        rapmap::utils::FixedWriter* cigarStr;
                        if (qa.mateStatus == MateStatus::PAIRED_END_LEFT) { // left read
                            alignedName = &readName;
                            unalignedName = &mateName;

                            readSeq = &(r.first.seq);
                            unalignedSeq = &(r.second.seq);

                            qstr = &(r.first.qual);
                            unalignedQstr = &(r.second.qual);

                            flags = flags1;
                            unalignedFlags = flags2;

                            cigarStr = &cigarStr1;

                            haveRev = &haveRev1;
                            readTemp = &read1Temp;
                            qualTemp = &qual1Temp;
                        } else { // right read
                            alignedName = &mateName;
                            unalignedName = &readName;

                            readSeq = &(r.second.seq);
                            unalignedSeq = &(r.first.seq);

                            qstr = &(r.second.qual);
                            unalignedQstr = &(r.first.qual);

                            flags = flags2;
                            unalignedFlags = flags1;

                            cigarStr = &cigarStr2;
                            haveRev = &haveRev2;
                            readTemp = &read2Temp;
                            qualTemp = &qual2Temp;
                        }

                        // Reverse complement the read and reverse
                        // the quality string if we need to
                        if (!qa.fwd) {
                            if (!(*haveRev)) {
                                rapmap::utils::reverseRead(*readSeq, *qstr,
                                        *readTemp, *qualTemp);
                                *haveRev = true;
                            }
                            readSeq = readTemp;
                            qstr = qualTemp;
                        }

                        /*
                        if (!qa.fwd) {
                            rapmap::utils::reverseRead(*readSeq, *qstr,
                                        read1Temp, qual1Temp);
                        }
                        */

                        rapmap::utils::adjustOverhang(qa.pos, qa.readLen, txpLens[qa.tid], *cigarStr);
                        sstream << alignedName->c_str() << '\t' // QNAME
                                << flags << '\t' // FLAGS
                                << transcriptName << '\t' // RNAME
                                << qa.pos + 1 << '\t' // POS (1-based)
                                << 1 << '\t' // MAPQ
                                << cigarStr->c_str() << '\t' // CIGAR
                                << '=' << '\t' // RNEXT
                                << qa.pos+1 << '\t' // PNEXT (only 1 read in templte)
                                << 0 << '\t' // TLEN (spec says 0, not read len)
                                << *readSeq << '\t' // SEQ
                                << *qstr << '\t' // QUAL
                                << numHitFlag << '\n';

                        // Output the info for the unaligned mate.
                        sstream << unalignedName->c_str() << '\t' // QNAME
                            << unalignedFlags << '\t' // FLAGS
                            << transcriptName << '\t' // RNAME (same as mate)
                            << qa.pos + 1 << '\t' // POS (same as mate)
                            << 0 << '\t' // MAPQ
                            << readSeq->length() << 'M' << '\t' // CIGAR
                            << '=' << '\t' // RNEXT
                            << qa.pos + 1 << '\t' // PNEXT (only 1 read in template)
                            << 0 << '\t' // TLEN (spec says 0, not read len)
                            << *unalignedSeq << '\t' // SEQ
                            << *unalignedQstr << '\t' // QUAL
                            << numHitFlag << '\n';
                    }
                    ++alnCtr;
                    // == SAM
#if defined(__DEBUG__) || defined(__TRACK_CORRECT__)
                    if (transcriptName == trueTxpName) {
							if (trueHitCtr == 0) {
									++hctr.trueHits;
									++trueHitCtr;
									firstTrueHit = &qa;
							} else {
									++trueHitCtr;
									std::cerr << "Found true hit " << trueHitCtr << " times!\n";
									std::cerr << transcriptName << '\t' << firstTrueHit->pos
											<< '\t' << firstTrueHit->fwd << '\t' << firstTrueHit->fragLen
											<< '\t' << (firstTrueHit->isPaired ? "Paired" : "Orphan") << '\t';
								    printMateStatus(firstTrueHit->mateStatus);
								    std::cerr << '\n';
									std::cerr << transcriptName << '\t' << qa.pos
											  << '\t' << qa.fwd << '\t' << qa.fragLen
										      << '\t' << (qa.isPaired ? "Paired" : "Orphan") << '\t';
								    printMateStatus(qa.mateStatus);
								    std::cerr << '\n';
							}
					}
#endif //__DEBUG__
                }
                return alnCtr;
        }



        // Is there a smarter way to do save / load here?
        /*
        template <typename Archive, typename MerT>
            void save(Archive& archive, const MerT& mer) const {
                auto key = mer.get_bits(0, 2*mer.k());
                archive(key);
            }

        template <typename Archive>
            void load(Archive& archive, const MerT& mer) {
                mer.polyT();
                uint64_t bits;
                archive(bits);
                auto k = mer.k();
                mer.set_bits(0, 2*k, bits);
            }
        */
    }
}

// Explicit instantiations
template uint32_t rapmap::utils::writeAlignmentsToStream<std::pair<header_sequence_qual, header_sequence_qual>, RapMapSAIndex*>(
                std::pair<header_sequence_qual, header_sequence_qual>& r,
                PairAlignmentFormatter<RapMapSAIndex*>& formatter,
                rapmap::utils::HitCounters& hctr,
                std::vector<rapmap::utils::QuasiAlignment>& jointHits,
                fmt::MemoryWriter& sstream
                );

template uint32_t rapmap::utils::writeAlignmentsToStream<std::pair<header_sequence_qual, header_sequence_qual>, RapMapIndex*>(
                std::pair<header_sequence_qual, header_sequence_qual>& r,
                PairAlignmentFormatter<RapMapIndex*>& formatter,
                rapmap::utils::HitCounters& hctr,
                std::vector<rapmap::utils::QuasiAlignment>& jointHits,
                fmt::MemoryWriter& sstream
                );
template uint32_t rapmap::utils::writeAlignmentsToStream<jellyfish::header_sequence_qual, RapMapSAIndex*>(
                jellyfish::header_sequence_qual& r,
                SingleAlignmentFormatter<RapMapSAIndex*>& formatter,
                rapmap::utils::HitCounters& hctr,
                std::vector<rapmap::utils::QuasiAlignment>& jointHits,
                fmt::MemoryWriter& sstream
                );
template uint32_t rapmap::utils::writeAlignmentsToStream<jellyfish::header_sequence_qual, RapMapIndex*>(
                jellyfish::header_sequence_qual& r,
                SingleAlignmentFormatter<RapMapIndex*>& formatter,
                rapmap::utils::HitCounters& hctr,
                std::vector<rapmap::utils::QuasiAlignment>& jointHits,
                fmt::MemoryWriter& sstream
                );


