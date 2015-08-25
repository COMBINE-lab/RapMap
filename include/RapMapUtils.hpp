#ifndef __RAP_MAP_UTILS_HPP__
#define __RAP_MAP_UTILS_HPP__

#include <atomic>
#include <memory>
#include "xxhash.h"
#include <cereal/archives/binary.hpp>
#include "jellyfish/mer_dna.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"
#include "PairSequenceParser.hpp"

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

// Must be forward-declared
template <typename IndexT>
class PairAlignmentFormatter;
template <typename IndexT>
class SingleAlignmentFormatter;

// Forward-declare because the C++ compiler is dumb
class RapMapSAIndex;
class RapMapIndex;

namespace rapmap {
    namespace utils {

    using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

    constexpr uint32_t newTxpSetMask = 0x80000000;
    constexpr uint32_t rcSetMask = 0x40000000;

    // Positions are stored in a packed format, where the highest
    // 2-bits encode if this position refers to a new transcript
    // and whether or not the k-mer from the hash matches this txp
    // in the forward or RC direction.
    void decodePosition(uint32_t p, uint32_t& pout, bool& newTxp, bool& isRC);

    template <typename IndexT>
        void writeSAMHeader(IndexT& rmi, std::shared_ptr<spdlog::logger> out) {
            fmt::MemoryWriter hd;
            hd.write("@HD\tVN:0.1\tSO:unknown\n");

            auto& txpNames = rmi.txpNames;
            auto& txpLens = rmi.txpLens;

            auto numRef = txpNames.size();
            for (size_t i = 0; i < numRef; ++i) {
                hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
            }
            // Eventuall output a @PG line
            //hd.format("@PG\t");
            std::string headerStr(hd.str());
            // Don't include the last '\n', since the logger will do it for us.
            headerStr.pop_back();
            out->info() << headerStr;
        }

    template <typename IndexT>
        void writeSAMHeader(IndexT& rmi, std::ostream& outStream) {
            fmt::MemoryWriter hd;
            hd.write("@HD\tVN:0.1\tSO:unknown\n");

            auto& txpNames = rmi.txpNames;
            auto& txpLens = rmi.txpLens;

            auto numRef = txpNames.size();
            for (size_t i = 0; i < numRef; ++i) {
                hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
            }
            // Eventuall output a @PG line
            //hd.format("@PG\t");
            outStream << hd.str();
        }

    // from http://stackoverflow.com/questions/9435385/split-a-string-using-c11
    std::vector<std::string> tokenize(const std::string &s, char delim);

    // from https://github.com/cppformat/cppformat/issues/105
    class FixedBuffer : public fmt::Buffer<char> {
        public:
            FixedBuffer(char *array, std::size_t size)
                : fmt::Buffer<char>(array, size) {}

        protected:
            void grow(std::size_t size) {
                throw std::runtime_error("buffer overflow");
            }
    };

    class FixedWriter : public fmt::Writer {
        private:
            FixedBuffer buffer_;
        public:
            FixedWriter(char *array, std::size_t size)
                : fmt::Writer(buffer_), buffer_(array, size) {}
    };



    struct SAInterval {
        uint32_t begin;
        uint32_t end;
        template <typename Archive>
            void load(Archive& ar) { ar(begin, end); }

        template <typename Archive>
            void save(Archive& ar) const { ar(begin, end); }
    };


    struct HitCounters {
        std::atomic<uint64_t> peHits{0};
        std::atomic<uint64_t> seHits{0};
        std::atomic<uint64_t> trueHits{0};
        std::atomic<uint64_t> totHits{0};
        std::atomic<uint64_t> numReads{0};
        std::atomic<uint64_t> tooManyHits{0};
        std::atomic<uint64_t> lastPrint{0};
    };

    class JFMerKeyHasher{
        public:
            size_t operator()(const my_mer& m) const {
                auto k = rapmap::utils::my_mer::k();
                auto v = m.get_bits(0, 2*k);
                return XXH64(static_cast<void*>(&v), 8, 0);
            }
    };

    class KmerKeyHasher {
        public:
            size_t operator()(const uint64_t& m) const {
                //auto k = rapmap::utils::my_mer::k();
                //auto v = m.get_bits(0, 2*k);
                auto v = m;
                return XXH64(static_cast<void*>(&v), 8, 0);
            }
    };

    struct KmerInterval {
        uint64_t offset;
        uint32_t length;

        template <typename Archive>
            void save(Archive& arch) const {
                arch(offset, length);
            }

        template <typename Archive>
            void load(Archive& arch) {
                arch(offset, length);
            }
    };

    struct KmerInfo {
        KmerInfo () : eqId(0), offset(0), count(0) {}


        KmerInfo(uint32_t eqIdIn, uint32_t offsetIn, uint32_t countIn) :
            eqId(eqIdIn), offset(offsetIn), count(countIn) {}

        template <typename Archive>
        void load(Archive& ar) {
            ar(eqId, offset, count);
        }

        template <typename Archive>
        void save(Archive& ar) const {
            ar(eqId, offset, count);
        }
        uint32_t eqId = 0;
        uint32_t offset = 0;
        uint32_t count = 0;
    };


    template <class T>
    inline void hashCombine(std::size_t& seed, const T& v)
    {
            std::hash<T> hasher;
            seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }

    constexpr uint32_t uint32Invalid = std::numeric_limits<uint32_t>::max();
    using PositionList = std::vector<uint32_t>;
    using KmerInfoList = std::vector<KmerInfo>;

    enum class MateStatus : uint8_t {
        SINGLE_END = 0,
        PAIRED_END_LEFT = 1,
        PAIRED_END_RIGHT = 2,
        PAIRED_END_PAIRED = 3 };

    // Wraps the standard iterator of the Position list to provide
    // some convenient functionality.  In the future, maybe this
    // should be a proper iterator adaptor.
    struct PositionListHelper{
        using PLIt = PositionList::iterator;

        PositionListHelper(PLIt itIn, PLIt endIn) :
            it_(itIn), end_(endIn) {}
        // The underlying iterator shouldn't be advanced further
        inline bool done() { return it_ == end_; }

        // The actual postion on the transcript
        int32_t pos() const { return static_cast<int32_t>((*it_) & 0x3FFFFFFF); }

        // True if the position encoded was on the reverse complement strand
        // of the reference transcript, false otherwise.
        bool isRC() const { return (*it_) & 0x40000000; }

        // True if we hit the position list for a new transcript, false otherwise
        bool isNewTxp() const { return (*it_) & 0x80000000; }

        void advanceToNextTranscript() {
            if (it_ < end_) {
                do {
                    ++it_;
                } while (!isNewTxp() and it_ != end_);

            }
        }

        PLIt it_; // The underlying iterator
        PLIt end_; // The end of the container
    };


    struct QuasiAlignment {
  	QuasiAlignment() :
		tid(std::numeric_limits<uint32_t>::max()),
		pos(std::numeric_limits<uint32_t>::max()),
		fwd(true),
		readLen(std::numeric_limits<uint32_t>::max()),
		fragLen(std::numeric_limits<uint32_t>::max()),
		isPaired(false) {}

        QuasiAlignment(uint32_t tidIn, uint32_t posIn,
                bool fwdIn, uint32_t readLenIn,
                uint32_t fragLenIn = 0,
                bool isPairedIn = false) :
            tid(tidIn), pos(posIn), fwd(fwdIn),
            readLen(readLenIn), fragLen(fragLenIn),
            isPaired(isPairedIn) {}
        QuasiAlignment(QuasiAlignment&& other) = default;
        QuasiAlignment& operator=(const QuasiAlignment&) = default;
        // Only 1 since the mate must have the same tid
        // we won't call *chimeric* alignments here.
        uint32_t tid;
        // Left-most position of the hit
        int32_t pos;
        // left-most position of the mate
        int32_t matePos;
        // Is the read from the forward strand
        bool fwd;
        // Is the mate from the forward strand
        bool mateIsFwd;
        // The fragment length (template length)
        // This is 0 for single-end or orphaned reads.
        uint32_t fragLen;
        // The read's length
        uint32_t readLen;
        // The mate's length
        uint32_t mateLen;
        // Is this a paired *alignment* or not
        bool isPaired;
        MateStatus mateStatus;
    };

    struct HitInfo {
        HitInfo(KmerInfoList::iterator kit, uint32_t merIDIn,
                int32_t queryPosIn, bool queryRCIn) :
                kinfo(kit), merID(merIDIn), queryPos(queryPosIn),
                queryRC(queryRCIn) {}

        KmerInfoList::iterator kinfo;
        uint32_t merID;
        int32_t queryPos;
        bool queryRC;
    };

    struct SAIntervalHit {
        SAIntervalHit(int beginIn, int endIn, uint32_t lenIn, uint32_t queryPosIn, bool queryRCIn) :
            begin(beginIn), end(endIn), len(lenIn), queryPos(queryPosIn), queryRC(queryRCIn) {}

	int span() { return end - begin; }
        int begin, end;
        uint32_t len, queryPos;
        bool queryRC;
    };

    struct SATxpQueryPos {
	SATxpQueryPos(uint32_t posIn, uint32_t qposIn, bool queryRCIn, bool activeIn = false) :
		pos(posIn), queryPos(qposIn), queryRC(queryRCIn), active(activeIn) {}
	uint32_t pos, queryPos;
	bool queryRC, active;
    };

    struct ProcessedSAHit {
	    ProcessedSAHit() : tid(std::numeric_limits<uint32_t>::max()), active(false), numActive(1) {}

	    ProcessedSAHit(uint32_t txpIDIn, uint32_t txpPosIn, uint32_t queryPosIn, bool queryRCIn) :
		    tid(txpIDIn), active(false), numActive(1)
	    {
		tqvec.emplace_back(txpPosIn, queryPosIn, queryRCIn);
	    }

	    uint32_t tid;
	    std::vector<SATxpQueryPos> tqvec;
            bool active;
	    uint32_t numActive;
    };

    struct SAHitInfo {
	    SAHitInfo(uint32_t txpIDIn, uint32_t txpPosIn, uint32_t queryPosIn, bool queryRCIn) :
		    tid(txpIDIn), pos(txpPosIn), queryPos(queryPosIn), queryRC(queryRCIn) {}
	    uint32_t tid;
	    uint32_t pos;
	    uint32_t queryPos;
	    bool queryRC;
    };

    struct TxpQueryPos {
        TxpQueryPos(PositionListHelper& ph, int32_t queryPosIn, bool queryRCIn) :
                txpPosInfo(ph), queryPos(queryPosIn), queryRC(queryRCIn) {}
        // Iterator for the beginning of the position list
        // of a given k-mer into a given transcript.
        PositionListHelper txpPosInfo;
        // The position of the k-mer on the query.
        int32_t queryPos;
        bool queryRC;
    };

    struct ProcessedHit {
        ProcessedHit() : tid(std::numeric_limits<uint32_t>::max()) {}
        ProcessedHit(uint32_t tidIn,
                     PositionListHelper ph, int32_t queryPos, bool queryRC) :
                     tid(tidIn) {
                         tqvec.emplace_back(ph, queryPos, queryRC);
                     }


        uint32_t tid; // transcript id
        // A vector of iterators into the position list
        // for the k-mers hitting this transcript
        std::vector<TxpQueryPos> tqvec;
    };


    struct EqClass {
        EqClass() :
            txpListStart(uint32Invalid), txpListLen(uint32Invalid) {}
        EqClass(uint32_t txpListStartIn, uint32_t txpListLenIn) :
            txpListStart(txpListStartIn), txpListLen(txpListLenIn) {}

        template <typename Archive>
        void load (Archive& ar) {
            ar(txpListStart, txpListLen);
        }

        template <typename Archive>
        void save (Archive& ar) const {
            ar(txpListStart, txpListLen);
        }

        uint32_t txpListStart;
        uint32_t txpListLen;
    };

    inline void printMateStatus(rapmap::utils::MateStatus ms) {
        switch(ms) {
            case rapmap::utils::MateStatus::SINGLE_END:
                std::cerr << "SINGLE END";
                break;
            case rapmap::utils::MateStatus::PAIRED_END_LEFT:
                std::cerr << "PAIRED END (LEFT)";
                break;
            case rapmap::utils::MateStatus::PAIRED_END_RIGHT:
                std::cerr << "PAIRED END (RIGHT)";
                break;
            case rapmap::utils::MateStatus::PAIRED_END_PAIRED:
                std::cerr << "PAIRED END (PAIRED)";
                break;
        }
    }


    // Declarations for functions dealing with SAM formatting and output
    //
    inline void adjustOverhang(int32_t& pos, uint32_t readLen,
		    uint32_t txpLen, FixedWriter& cigarStr) {
	    cigarStr.clear();
	    if (pos + readLen < 0) {
            cigarStr.write("{}S", readLen);
            pos = 0;
        } else if (pos < 0) {
		    int32_t matchLen = readLen + pos;
            int32_t clipLen = readLen - matchLen;
		    cigarStr.write("{}S{}M", clipLen, matchLen);
		    // Now adjust the mapping position
		    pos = 0;
	    } else if (pos > txpLen) {
            cigarStr.write("{}S", readLen);
        } else if (pos + readLen > txpLen) {
		    int32_t matchLen = txpLen - pos;
		    int32_t clipLen = readLen - matchLen;
		    cigarStr.write("{}M{}S", matchLen, clipLen);
	    } else {
		    cigarStr.write("{}M", readLen);
	    }
    }

    inline void adjustOverhang(QuasiAlignment& qa, uint32_t txpLen,
		    FixedWriter& cigarStr1, FixedWriter& cigarStr2) {
	    if (qa.isPaired) { // both mapped
		    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);
		    adjustOverhang(qa.matePos, qa.mateLen, txpLen, cigarStr2);
	    } else if (qa.mateStatus == MateStatus::PAIRED_END_LEFT ) {
		    // left read mapped
		    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr1);
	    } else if (qa.mateStatus == MateStatus::PAIRED_END_RIGHT) {
		    // right read mapped
		    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr2);
	    }
    }



        // get the sam flags for the quasialignment qaln.
        // peinput is true if the read is paired in *sequencing*; false otherwise
        // the sam flags for mate 1 are written into flags1 and for mate2 into flags2
        inline void getSamFlags(const QuasiAlignment& qaln,
                uint16_t& flags) {
            constexpr uint16_t pairedInSeq = 0x1;
            constexpr uint16_t mappedInProperPair = 0x2;
            constexpr uint16_t unmapped = 0x4;
            constexpr uint16_t mateUnmapped = 0x8;
            constexpr uint16_t isRC = 0x10;
            constexpr uint16_t mateIsRC = 0x20;
            constexpr uint16_t isRead1 = 0x40;
            constexpr uint16_t isRead2 = 0x80;
            constexpr uint16_t isSecondaryAlignment = 0x100;
            constexpr uint16_t failedQC = 0x200;
            constexpr uint16_t isPCRDup = 0x400;
            constexpr uint16_t supplementaryAln = 0x800;

            flags = 0;
            // Not paired in sequencing
            // flags1 = (peInput) ? pairedInSeq : 0;
            // flags |= properlyAligned;
            // we don't output unmapped yet
            // flags |= unmapped
            // flags |= mateUnmapped
            flags |= (qaln.fwd) ? 0 : isRC;
            // Mate flag meaningless
            // flags1 |= (qaln.mateIsFwd) ? 0 : mateIsRC;
            // flags |= isRead1;
            //flags2 |= isRead2;
        }

        // get the sam flags for the quasialignment qaln.
        // peinput is true if the read is paired in *sequencing*; false otherwise
        // the sam flags for mate 1 are written into flags1 and for mate2 into flags2
        inline void getSamFlags(const QuasiAlignment& qaln,
                bool peInput,
                uint16_t& flags1,
                uint16_t& flags2) {
            constexpr uint16_t pairedInSeq = 0x1;
            constexpr uint16_t properlyAligned = 0x2;
            constexpr uint16_t unmapped = 0x4;
            constexpr uint16_t mateUnmapped = 0x8;
            constexpr uint16_t isRC = 0x10;
            constexpr uint16_t mateIsRC = 0x20;
            constexpr uint16_t isRead1 = 0x40;
            constexpr uint16_t isRead2 = 0x80;
            constexpr uint16_t isSecondaryAlignment = 0x100;
            constexpr uint16_t failedQC = 0x200;
            constexpr uint16_t isPCRDup = 0x400;
            constexpr uint16_t supplementaryAln = 0x800;

            flags1 = flags2 = 0;
            flags1 = (peInput) ? pairedInSeq : 0;
            flags1 |= (qaln.isPaired) ? properlyAligned : 0;
            flags2 = flags1;
            // we don't output unmapped yet
            bool read1Unaligned = qaln.mateStatus == MateStatus::PAIRED_END_RIGHT;
            bool read2Unaligned = qaln.mateStatus == MateStatus::PAIRED_END_LEFT;
            // If read 1 is unaligned, flags1 gets "unmapped" and flags2 gets "mate unmapped"
            flags1 |= (read1Unaligned) ? unmapped : 0;
            flags2 |= (read1Unaligned) ? mateUnmapped : 0;
            // If read 2 is unaligned, flags2 gets "unmapped" and flags1 gets "mate unmapped"
            flags2 |= (read2Unaligned) ? unmapped : 0;
            flags1 |= (read2Unaligned) ? mateUnmapped : 0;

            flags1 |= (qaln.fwd) ? 0 : isRC;
            flags1 |= (qaln.mateIsFwd) ? 0 : mateIsRC;
            flags2 |= (qaln.mateIsFwd) ? 0 : isRC;
            flags2 |= (qaln.fwd) ? 0 : mateIsRC;
            flags1 |= isRead1;
            flags2 |= isRead2;
        }

	// Adapted from
        // https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/blob/8c9933a1685e0ab50c7d8b7926c9068bc0c9d7d2/src/main.c#L36
        void reverseRead(std::string& seq,
                std::string& qual,
                std::string& readWork,
                std::string& qualWork);

        template <typename ReadPairT, typename IndexT>
        uint32_t writeAlignmentsToStream(
                ReadPairT& r,
                PairAlignmentFormatter<IndexT>& formatter,
                HitCounters& hctr,
                std::vector<QuasiAlignment>& jointHits,
                fmt::MemoryWriter& sstream);

        template <typename ReadT, typename IndexT>
        uint32_t writeAlignmentsToStream(
                ReadT& r,
                SingleAlignmentFormatter<IndexT>& formatter,
                HitCounters& hctr,
                std::vector<QuasiAlignment>& jointHits,
                fmt::MemoryWriter& sstream);


        inline void mergeLeftRightHits(
                std::vector<QuasiAlignment>& leftHits,
                std::vector<QuasiAlignment>& rightHits,
                std::vector<QuasiAlignment>& jointHits,
                uint32_t readLen,
                uint32_t maxNumHits,
                bool& tooManyHits,
                HitCounters& hctr) {
            if (leftHits.size() > 0) {
                auto leftIt = leftHits.begin();
                auto leftEnd = leftHits.end();
                auto leftLen = std::distance(leftIt, leftEnd);
                if (rightHits.size() > 0) {
                    auto rightIt = rightHits.begin();
                    auto rightEnd = rightHits.end();
                    auto rightLen = std::distance(rightIt, rightEnd);
                    size_t numHits{0};
                    jointHits.reserve(std::min(leftLen, rightLen));
                    uint32_t leftTxp, rightTxp;
                    while (leftIt != leftEnd && rightIt != rightEnd) {
                        leftTxp = leftIt->tid;
                        rightTxp = rightIt->tid;
                        if (leftTxp < rightTxp) {
                            ++leftIt;
                        } else {
                            if (!(rightTxp < leftTxp)) {
                                int32_t startRead1 = leftIt->pos;
                                int32_t startRead2 = rightIt->pos;
                                int32_t fragStartPos = std::min(leftIt->pos, rightIt->pos);
                                int32_t fragEndPos = std::max(leftIt->pos, rightIt->pos) + readLen;
                                uint32_t fragLen = fragEndPos - fragStartPos;
                                jointHits.emplace_back(leftTxp,
                                        startRead1,
                                        leftIt->fwd,
                                        leftIt->readLen,
                                        fragLen, true);
                                // Fill in the mate info
                                auto& qaln = jointHits.back();
                                qaln.mateLen = rightIt->readLen;
                                qaln.matePos = startRead2;
                                qaln.mateIsFwd = rightIt->fwd;
                                jointHits.back().mateStatus = MateStatus::PAIRED_END_PAIRED;

                                ++numHits;
                                if (numHits > maxNumHits) { tooManyHits = true; break; }
                                ++leftIt;
                            }
                            ++rightIt;
                        }
                    }
                }
                if (tooManyHits) { jointHits.clear(); ++hctr.tooManyHits; }
            }

            // If we had proper paired hits
            if (jointHits.size() > 0) {
                hctr.peHits += jointHits.size();
                //orphanStatus = 0;
            } else if (leftHits.size() + rightHits.size() > 0 and !tooManyHits) {
                // If there weren't proper paired hits, then either
                // there were too many hits, and we forcibly discarded the read
                // or we take the single end hits.
                auto numHits = leftHits.size() + rightHits.size();
                hctr.seHits += numHits;
                //orphanStatus = 0;
                //orphanStatus |= (leftHits.size() > 0) ? 0x1 : 0;
                //orphanStatus |= (rightHits.size() > 0) ? 0x2 : 0;
                jointHits.insert(jointHits.end(),
                        std::make_move_iterator(leftHits.begin()),
                        std::make_move_iterator(leftHits.end()));
                jointHits.insert(jointHits.end(),
                        std::make_move_iterator(rightHits.begin()),
                        std::make_move_iterator(rightHits.end()));
            }
        }

    /*
    template <typename Archive>
    void save(Archive& archive, const my_mer& mer);

    template <typename Archive>
    void load(Archive& archive, my_mer& mer);
    */
    }
}


#endif // __RAP_MAP_UTILS_HPP__
