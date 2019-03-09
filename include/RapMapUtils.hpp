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

#ifndef __RAP_MAP_UTILS_HPP__
#define __RAP_MAP_UTILS_HPP__

#include <atomic>
#include <cmath>
#include <memory>
#include "xxhash.h"
#include "sparsepp/spp.h"
#include "SparseHashSerializer.hpp"
#include <cereal/archives/binary.hpp>
#include "Kmer.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"
#include "chobo/small_vector.hpp"
#include "RapMapConfig.hpp"
#include "nonstd/optional.hpp"
#include "FastxParser.hpp"

#ifdef RAPMAP_SALMON_SUPPORT
#include "LibraryFormat.hpp"
#endif

#ifndef __DEFINE_LIKELY_MACRO__
#define __DEFINE_LIKELY_MACRO__
#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x), 1)
#define UNLIKELY(x) __builtin_expect((x), 0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif
#endif

// Must be forward-declared
template <typename IndexT>
struct PairAlignmentFormatter;
template <typename IndexT>
struct SingleAlignmentFormatter;

// Forward-declare because the C++ compiler is dumb
// class RapMapIndex;

template<typename KeyT, typename ValT, typename HasherT>
//using RegHashT = google::dense_hash_map<KeyT, ValT, HasherT>;
using RegHashT = spp::sparse_hash_map<KeyT, ValT, HasherT>;

template<typename KeyT, typename ValT>
class FrugalBooMap;

template<typename KeyT, typename ValT>
using PerfectHashT = FrugalBooMap<KeyT, ValT>;

namespace rapmap {
    namespace utils {

    using my_mer = combinelib::kmers::Kmer<32,1>;

    constexpr uint32_t newTxpSetMask = 0x80000000;
    constexpr uint32_t rcSetMask = 0x40000000;

    class MappingConfig {
    public:
      bool consistentHits{false};
      bool doChaining{false};
      float consensusFraction{1.0};
      bool considerMultiPos{false};
      bool allowDovetail{false};
    };

    // Positions are stored in a packed format, where the highest
    // 2-bits encode if this position refers to a new transcript
    // and whether or not the k-mer from the hash matches this txp
    // in the forward or RC direction.
    void decodePosition(uint32_t p, uint32_t& pout, bool& newTxp, bool& isRC);

    template <typename IndexT>
        void writeSAMHeader(IndexT& rmi, std::shared_ptr<spdlog::logger> out) {
            fmt::MemoryWriter hd;
	    hd.write("@HD\tVN:1.0\tSO:unknown\n");

            auto& txpNames = rmi.txpNames;
            auto& txpLens = rmi.txpLens;

            auto numRef = txpNames.size();
            for (size_t i = 0; i < numRef; ++i) {
                hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
            }
            // Eventually output a @PG line
            hd.write("@PG\tID:rapmap\tPN:rapmap\tVN:{}\n", rapmap::version);
            std::string headerStr(hd.str());
            // Don't include the last '\n', since the logger will do it for us.
            headerStr.pop_back();
            out->info(headerStr);
        }

    template <typename IndexT>
        void writeSAMHeader(IndexT& rmi, std::ostream& outStream) {
            fmt::MemoryWriter hd;
	    hd.write("@HD\tVN:1.0\tSO:unknown\n");

            auto& txpNames = rmi.txpNames;
            auto& txpLens = rmi.txpLens;

            auto numRef = txpNames.size();
            for (size_t i = 0; i < numRef; ++i) {
                hd.write("@SQ\tSN:{}\tLN:{:d}\n", txpNames[i], txpLens[i]);
            }
            // Eventually output a @PG line
            hd.write("@PG\tID:rapmap\tPN:rapmap\tVN:0.3.1\n");
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

    /**
     * Stores both the key (k-mer)
     * and the interval to which it corresponds.
     * This is useful if the hash itself doesn't validate
     * the key (e.g. a minimum perfect hash).
     **/
    template <typename IndexT>
    struct SAIntervalWithKey {
        uint64_t kmer;
      //  SAInterval<IndexT> second;
        IndexT begin_;
        IndexT end_;

        inline IndexT begin() const { return begin_; }
        inline IndexT end() const { return end_; }

        template <typename Archive>
            void load(Archive& ar) { ar(kmer, begin_, end_); }

        template <typename Archive>
            void save(Archive& ar) const { ar(kmer, begin_, end_); }
    };

    template <typename IndexT>
    struct SAInterval {
      /*
        SAInterval(IndexT beginIn, IndexT endIn) : begin(beginIn), end(endIn) {}
	SAInterval(std::initializer_list<IndexT> il) {
	  auto it = il.begin();
	  begin = *(it);
	  ++it;
	  end = *(il.begin());
	}
	*/
        using index_type = IndexT;
        IndexT begin_;
        IndexT end_;
        
        inline IndexT begin() const { return begin_; }
        inline IndexT end() const { return end_; }

        template <typename Archive>
        void load(Archive& ar) { ar(begin_, end_); }
        //void load(Archive& ar) { ar(begin_, len_); }

        template <typename Archive>
        void save(Archive& ar) const { ar(begin_, end_); }
        //void save(Archive& ar) const { ar(begin_, len_); }
    };


    struct HitCounters {
        std::atomic<uint64_t> peHits{0};
        std::atomic<uint64_t> seHits{0};
        std::atomic<uint64_t> trueHits{0};
        std::atomic<uint64_t> totHits{0};
        std::atomic<uint64_t> numReads{0};
        std::atomic<uint64_t> tooManyHits{0};
        std::atomic<uint64_t> lastPrint{0};
        std::atomic<uint64_t> numDovetails{0};
    };

    class JFMerKeyHasher{
        public:
            size_t operator()(const my_mer& m) const {
                auto v = m.word(0);//get_bits(0, 2*k);
                return XXH64(static_cast<void*>(&v), 8, 0);
            }
    };

    class KmerKeyHasher {
        //spp::spp_hash<uint64_t> hasher;
        public:
        //inline size_t operator()(const uint64_t& m) const { //{ return hasher(m); }
        inline size_t operator()(const rapmap::utils::my_mer& m) const { //{ return hasher(m); }
                //auto k = rapmap::utils::my_mer::k();
                //auto v = m.get_bits(0, 2*k);
                //auto v = m;
            return XXH64(static_cast<void*>(const_cast<rapmap::utils::my_mer::base_type*>(m.data())), sizeof(m.word(0)) * m.nb_words(), 0);
            }
        inline size_t operator()(const uint64_t& m) const { //{ return hasher(m); }
            return XXH64(static_cast<void*>(const_cast<uint64_t*>(&m)), sizeof(m), 0);
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

    enum class ChainStatus : uint8_t {
      PERFECT = 0,
      UNGAPPED = 1,
      ALIGNED_ON_LEFT = 2,
      ALIGNED_ON_RIGHT = 3,
      REGULAR = 4
    };

      inline std::ostream& operator<<(std::ostream& os, ChainStatus s) {
        switch (s) {
        case ChainStatus::PERFECT:
          os << "PERFECT";
          break;
        case ChainStatus::UNGAPPED:
          os << "UNGAPPED";
          break;
        case ChainStatus::ALIGNED_ON_LEFT:
          os << "ALIGNED_ON_LEFT";
          break;
        case ChainStatus::ALIGNED_ON_RIGHT:
          os << "ALIGNED_ON_RIGHT";
          break;
        case ChainStatus::REGULAR:
          os << "REGULAR";
          break;
        default:
          os << "UNKNOWN";
          break;
        }
        return os;
      }

    class FragmentChainStatus {
      public:
        FragmentChainStatus() :
          left(static_cast<typename std::underlying_type<ChainStatus>::type>(ChainStatus::REGULAR)),
          right(static_cast<typename std::underlying_type<ChainStatus>::type>(ChainStatus::REGULAR))
        {}
        FragmentChainStatus(ChainStatus ls, ChainStatus rs) :
          left(static_cast<typename std::underlying_type<ChainStatus>::type>(ls)),
          right(static_cast<typename std::underlying_type<ChainStatus>::type>(rs))
        {}

        void setLeft(ChainStatus s) {
          left = static_cast<typename std::underlying_type<ChainStatus>::type>(s);
        }
        void setRight(ChainStatus s) {
          right = static_cast<typename std::underlying_type<ChainStatus>::type>(s);
        }

      ChainStatus getLeft() const {
        return static_cast<ChainStatus>(left);
      }

      ChainStatus getRight() const {
        return static_cast<ChainStatus>(right);
      }

      friend inline std::ostream& operator<<(std::ostream& os, const FragmentChainStatus& s) {
        os << "fragment chain status {" << s.getLeft() << ", " << s.getRight() << "}";
        return os;
      }

      private:
        uint8_t left : 4, right : 4;
    };

      enum class MateStatus : uint8_t {
        SINGLE_END = 0,
        PAIRED_END_LEFT = 1,
        PAIRED_END_RIGHT = 2,
        PAIRED_END_PAIRED = 3,
        NOTHING = std::numeric_limits<uint8_t>::max()
    };

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
		pos(std::numeric_limits<int32_t>::max()),
		fwd(true),
		fragLen(std::numeric_limits<uint32_t>::max()),
		readLen(std::numeric_limits<uint32_t>::max()),
		isPaired(false)
#ifdef RAPMAP_SALMON_SUPPORT
        ,format(LibraryFormat::formatFromID(0))
#endif // RAPMAP_SALMON_SUPPORT
        {}

        QuasiAlignment(uint32_t tidIn, int32_t posIn,
                bool fwdIn, uint32_t readLenIn,
                uint32_t fragLenIn = 0,
                bool isPairedIn = false) :
            tid(tidIn), pos(posIn), fwd(fwdIn),
            fragLen(fragLenIn), readLen(readLenIn), 
            isPaired(isPairedIn)
#ifdef RAPMAP_SALMON_SUPPORT
        ,format(LibraryFormat::formatFromID(0))
#endif // RAPMAP_SALMON_SUPPORT
        {}
        QuasiAlignment(QuasiAlignment&& other) = default;
        QuasiAlignment& operator=(QuasiAlignment&) = default;
        QuasiAlignment& operator=(QuasiAlignment&& o) = default;
        QuasiAlignment(const QuasiAlignment& o) = default;
        QuasiAlignment(QuasiAlignment& o) = default;

      inline void setChainScore(double chainScoreIn) {
        chainScore_ = chainScoreIn;
      }

      inline double chainScore() const {
        return chainScore_;
      }

      inline uint32_t transcriptID() const { return tid; }
      inline double score() const { return score_; }
      inline void score(double scoreIn) { score_ = scoreIn; }
      inline int32_t alnScore() const { return alnScore_; }
      inline void alnScore(int32_t alnScoreIn) { alnScore_ = alnScoreIn; }
      inline uint32_t fragLength() const { return fragLen; }
      inline int32_t hitPos() { return std::min(pos, matePos); }

// Some convenience functions to allow salmon interop
#ifdef RAPMAP_SALMON_SUPPORT
      inline uint32_t fragLengthPedantic(uint32_t txpLen) const {
        if (mateStatus != rapmap::utils::MateStatus::PAIRED_END_PAIRED
            or fwd == mateIsFwd) {
          return 0;
        }
        int32_t p1 = fwd ? pos : matePos;
        int32_t sTxpLen = static_cast<int32_t>(txpLen);
        p1 = (p1 < 0) ? 0 : p1;
        p1 = (p1 > sTxpLen) ? sTxpLen : p1;
        int32_t p2 = fwd ? matePos + mateLen : pos + readLen;
        p2 = (p2 < 0) ? 0 : p2;
        p2 = (p2 > sTxpLen) ? sTxpLen : p2;

        return (p1 > p2) ? p1 - p2 : p2 - p1;
      }

      double logProb{HUGE_VAL};
      double logBias{HUGE_VAL};
      inline LibraryFormat libFormat() { return format; }
      LibraryFormat format;
#endif // RAPMAP_SALMON_SUPPORT
       bool hasMultiPos{false};
       chobo::small_vector<int32_t> allPositions;
       chobo::small_vector<int32_t> oppositeStrandPositions;

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
        // numeric score associated with this mapping
        double score_{1.0};
        // actual ``alignment'' score associated with this mapping.
        int32_t alnScore_{0};
        // If one or both of the reads is a complete match (no mismatch, indels), say what kind.
        FragmentChainStatus chainStatus;
        double chainScore_{std::numeric_limits<double>::lowest()};
      //int32_t queryOffset{-1};
      //MateStatus completeMatchType{MateStatus::NOTHING};
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

    template <typename OffsetT>
    struct SAIntervalHit {
        SAIntervalHit(OffsetT beginIn, OffsetT endIn, uint32_t lenIn, uint32_t queryPosIn, bool queryRCIn) :
            begin(beginIn), end(endIn), len(lenIn), queryPos(queryPosIn), queryRC(queryRCIn) {}

	      OffsetT span() { return end - begin; }
        OffsetT begin, end;
        uint32_t len, queryPos;
        bool queryRC;
    };

      struct SATxpQueryPos {
        SATxpQueryPos(uint32_t posIn, uint32_t qposIn, bool queryRCIn, /*bool activeIn = false,*/ int32_t lenIn = -1) :
          pos(posIn), queryPos(qposIn), queryRC(queryRCIn), /*active(activeIn),*/ len(lenIn){}
        uint32_t pos, queryPos;
        bool queryRC;//, active;
        int32_t len;
      };

    struct ProcessedSAHit {
	    ProcessedSAHit() : tid(std::numeric_limits<uint32_t>::max()), active(false), numActive(1) {}

	    ProcessedSAHit(uint32_t txpIDIn, uint32_t txpPosIn, uint32_t queryPosIn, bool queryRCIn, uint32_t lenIn) :
		    tid(txpIDIn), active(false), numActive(1)
	    {
        tqvec.emplace_back(txpPosIn, queryPosIn, queryRCIn, lenIn);
	    }

        /**
         * This enforces a more stringent consistency check on
         * the hits for this transcript.  The hits must be co-linear
         * with respect to the query and target.
         *
         * input: numToCheck --- the number of hits to check in sorted order
         *                       hits after the last of these need not be consistent.
         * return: numToCheck if the first numToCheck hits are consistent;
         *         -1 otherwise
         **/
        int32_t checkConsistent(size_t readLen, int32_t numToCheck) {
            auto numHits = tqvec.size();

            // special case for only 1 or two hits (common)
            if (numHits == 1) {
                return numToCheck;
            } else if (numHits == 2) {
                auto& h1 = (tqvec[0].queryPos < tqvec[1].queryPos) ? tqvec[0] : tqvec[1];
                auto& h2 = (tqvec[0].queryPos < tqvec[1].queryPos) ? tqvec[1] : tqvec[0];
                if (h2.pos > h1.pos) {
                    int32_t distortion = (h2.pos - h1.pos) - (h2.queryPos - h1.queryPos);
                    return (distortion > -10 and distortion < 10) ? numToCheck : -1;
                } else {
                    return -1;
                }
                //return (h2.pos > h1.pos) ? (numToCheck) : -1;
            } else {
                // first, sort by query position
                std::sort(tqvec.begin(), tqvec.end(),
                          [](const SATxpQueryPos& q1, const SATxpQueryPos& q2) -> bool {
                              return q1.queryPos < q2.queryPos;
                          });

                int32_t lastRefPos{std::numeric_limits<int32_t>::min()};
                int32_t lastQueryPos{std::numeric_limits<int32_t>::min()};
                bool firstHit{true};
                //int32_t maxDistortion{0};
                for (int32_t i = 0; i < numToCheck; ++i) {
                    int32_t refPos = static_cast<int32_t>(tqvec[i].pos);
                    int32_t queryPos = static_cast<int32_t>(tqvec[i].queryPos);
                    if (refPos > lastRefPos) {
                        int32_t distortion = 
                            firstHit ? 0 : ((refPos - lastRefPos) - (queryPos - lastQueryPos));
                        firstHit = false;
                        if (distortion < -10 or distortion > 10) {
                            return i;
                        }
                        lastRefPos = refPos;
                        lastQueryPos = queryPos;
                    } else {
                        return i;
                    }
                }
                return numToCheck;
            }
        }

	    uint32_t tid;
	    std::vector<SATxpQueryPos> tqvec;
      bool active;
      uint32_t numActive;
      uint32_t lastActiveInterval{1};
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
          case rapmap::utils::MateStatus::NOTHING:
            std::cerr << "NOTHING";
            break;
        }
    }


    // Declarations for functions dealing with SAM formatting and output
    //
    inline void adjustOverhang(int32_t& pos, uint32_t readLen,
		    uint32_t txpLen, FixedWriter& cigarStr) {
      int32_t sTxpLen = static_cast<int32_t>(txpLen);
      int32_t sReadLen = static_cast<int32_t>(readLen);
	    cigarStr.clear();
	    if (pos + static_cast<int32_t>(readLen) < 0) {
            cigarStr.write("{}S", readLen);
            pos = 0;
        } else if (pos < 0) {
		    int32_t matchLen = readLen + pos;
            int32_t clipLen = readLen - matchLen;
		    cigarStr.write("{}S{}M", clipLen, matchLen);
		    // Now adjust the mapping position
		    pos = 0;
	    } else if (pos > sTxpLen) {
            cigarStr.write("{}S", readLen);
        } else if (pos + sReadLen > sTxpLen) {
		    int32_t matchLen = sTxpLen - pos;
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
		    // right read un-mapped will just be read length * S
		    cigarStr2.clear();
		    cigarStr2.write("{}S", qa.mateLen);
	    } else if (qa.mateStatus == MateStatus::PAIRED_END_RIGHT) {
		    // right read mapped
		    adjustOverhang(qa.pos, qa.readLen, txpLen, cigarStr2);
		    // left read un-mapped will just be read length * S
		    cigarStr1.clear();
		    cigarStr1.write("{}S", qa.readLen);
	    }
    }



        // get the sam flags for the quasialignment qaln.
        // peinput is true if the read is paired in *sequencing*; false otherwise
        // the sam flags for mate 1 are written into flags1 and for mate2 into flags2
        inline void getSamFlags(const QuasiAlignment& qaln,
                uint16_t& flags) {
          /*
            constexpr uint16_t pairedInSeq = 0x1;
            constexpr uint16_t mappedInProperPair = 0x2;
            constexpr uint16_t unmapped = 0x4;
            constexpr uint16_t mateUnmapped = 0x8;
          */
            constexpr uint16_t isRC = 0x10;
            /*
            constexpr uint16_t mateIsRC = 0x20;
            constexpr uint16_t isRead1 = 0x40;
            constexpr uint16_t isRead2 = 0x80;
            constexpr uint16_t isSecondaryAlignment = 0x100;
            constexpr uint16_t failedQC = 0x200;
            constexpr uint16_t isPCRDup = 0x400;
            constexpr uint16_t supplementaryAln = 0x800;
            */

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
            //constexpr uint16_t isSecondaryAlignment = 0x100;
            //constexpr uint16_t failedQC = 0x200;
            //constexpr uint16_t isPCRDup = 0x400;
            //constexpr uint16_t supplementaryAln = 0x800;

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

        void reverseRead(std::string& seq,
                         std::string& readWork);


        std::string reverseComplement(std::string& seq);


        uint32_t writeUnalignedSingleToStream(fastx_parser::ReadSeq& r, fmt::MemoryWriter& sstream);
        uint32_t writeUnalignedPairToStream(fastx_parser::ReadPair& r, fmt::MemoryWriter& sstream);

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

        inline MateStatus mergeMatchType(MateStatus leftT, MateStatus rightT) {
          if (leftT == MateStatus::NOTHING) {
            return rightT;
          }
          if (rightT == MateStatus::NOTHING) {
            return leftT;
          }
          return MateStatus::PAIRED_END_PAIRED;
        }

      enum class MergeResult : uint8_t {
        HAD_NONE,
        HAD_EMPTY_INTERSECTION,
        HAD_CONCORDANT,
        HAD_DISCORDANT,
        HAD_ONLY_LEFT,
        HAD_ONLY_RIGHT,
      };

        inline MergeResult mergeLeftRightHitsFuzzy(
                bool leftMatches,
                bool rightMatches,
                std::vector<QuasiAlignment>& leftHits,
                std::vector<QuasiAlignment>& rightHits,
                std::vector<QuasiAlignment>& jointHits,
                rapmap::utils::MappingConfig& mc,
                uint32_t readLen,
                uint32_t maxNumHits,
                bool& tooManyHits,
                HitCounters& hctr) {

          using rapmap::utils::MergeResult;
          MergeResult mergeRes{MergeResult::HAD_NONE};

          const constexpr int32_t dovetailPenalty = std::numeric_limits<int32_t>::max() / 2;
          bool considerMultiPos = mc.considerMultiPos;
          bool allowDovetail = mc.allowDovetail;
            if (leftHits.empty()) {
                if (!leftMatches) {
                    if (!rightHits.empty()) {
                        jointHits.insert(jointHits.end(),
                                std::make_move_iterator(rightHits.begin()),
                                std::make_move_iterator(rightHits.end()));
                        hctr.seHits += rightHits.size();
                        mergeRes = MergeResult::HAD_ONLY_RIGHT;
                    }
                }
            } else if (rightHits.empty()) {
                if (!rightMatches) {
                    if (!leftHits.empty()) {
                        jointHits.insert(jointHits.end(),
                                std::make_move_iterator(leftHits.begin()),
                                std::make_move_iterator(leftHits.end()));
                        hctr.seHits += leftHits.size();
                        mergeRes = MergeResult::HAD_ONLY_LEFT;
                    }
                }
            } else {
                bool bestMappingIsDovetail = true;
                bool hadOppositeStrandMapping = false;
                constexpr const int32_t signedZero{0};
                uint32_t sameTxpCount{0};
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
                              ++sameTxpCount;

                              // returned tuple is fwPos, rcPos, gapLength
                              auto findBestHitFWRC = [signedZero, considerMultiPos, allowDovetail, &hadOppositeStrandMapping](
                                                                       chobo::small_vector<int32_t>& fwdHits,
                                                                       chobo::small_vector<int32_t>& rcHits,
                                                                       int32_t fwdReadLen,
                                                                       bool& bestHitIsDovetail) ->
                                nonstd::optional<std::tuple<int32_t, int32_t, int32_t>> {

                                // If either of the position vectors is empty, there can be no valid
                                // mapping.
                                if (fwdHits.empty() or rcHits.empty()) {
                                  bestHitIsDovetail = true;
                                  return nonstd::nullopt;
                                }
                                hadOppositeStrandMapping = true;

                                // Remember the pair of positions that gives us the best gap
                                // here, a gap of 0 is "optimal".
                                int32_t bestGap = std::numeric_limits<int32_t>::max();
                                auto bestFWPosIt = fwdHits.begin();
                                auto bestRCPosIt = rcHits.begin();

                                // NOTE: Do we need an explicit fast path here?
                                // if (considerMultiPos and (fwdHits.size() > 1 or rcHits.size() > 1))

                                // Given a left position and a right position, if they produce a better
                                // gap than the current best, then update the best gap and remember these
                                // positions that produced it.
                                auto updateBestGap = [&bestGap, &bestFWPosIt, &bestRCPosIt,
                                                      signedZero, fwdReadLen, allowDovetail](
                                                                                   chobo::small_vector<int32_t>::iterator fwdPosIt,
                                                                                   chobo::small_vector<int32_t>::iterator rcPosIt
                                                                                   ) {
                                  // The gap between the end of the first read and the start of the
                                  // second (we take the absolute value so it is always non-negative,
                                  // even if they overlap).

                                  // we expect the rc read to be "downstream" of the fwd read.
                                  constexpr int32_t maxGap = std::numeric_limits<int32_t>::max();

                                  int32_t gap = maxGap;
                                  int32_t fwdPos = *fwdPosIt;
                                  int32_t rcPos = *rcPosIt;

                                  // If the rc read is already downstream of the fwd read, then it's not dovetailed
                                  // so figure out the gap.
                                  if (rcPos >= fwdPos) {
                                    // NOTE: Can think harder about what the best measure of the gap penalty is.
                                    gap = std::abs(rcPos - (fwdPos + static_cast<int32_t>(fwdReadLen)));
                                  } else {
                                    // If the rc read is upstream of the fwd read, then compute the gap with dovetail penalty
                                    // if dovetails are allowed; otherwise leave it as maxGap.
                                    gap = (fwdPos - rcPos) + dovetailPenalty;
                                  }

                                  // if this is the best gap so far
                                  if (gap < bestGap) {
                                    bestGap = gap;
                                    bestFWPosIt = fwdPosIt;
                                    bestRCPosIt = rcPosIt;
                                  }
                                };


                                auto rcBeg = rcHits.begin(); auto rcEnd = rcHits.end();

                                // for every position the forward read could start
                                for (auto pIt = fwdHits.begin(); pIt != fwdHits.end(); ++pIt) {
                                  auto p1 = *pIt;

                                  // find the closest position for the rc read
                                  auto lbIt = std::lower_bound(rcBeg, rcEnd, p1);

                                  // p1 is greater than every position where the rc read can start
                                  if (lbIt == rcEnd) {
                                    auto closestIt = lbIt - 1;
                                    updateBestGap(pIt, closestIt);
                                  } else if (lbIt == rcBeg) {
                                    // every position where the rc read can start is greater than p1
                                    updateBestGap(pIt, lbIt);
                                  } else {
                                    // check the current element
                                    updateBestGap(pIt, lbIt);
                                    // and the previous
                                    auto prevIt = lbIt - 1;
                                    updateBestGap(pIt, prevIt);
                                  }
                                }

                                // if we had a valid gap, return the best gap
                                bestHitIsDovetail = (bestGap >= dovetailPenalty);
                                if (allowDovetail) {
                                  return (bestGap ==  std::numeric_limits<int32_t>::max()) ? nonstd::nullopt :
                                  nonstd::optional<std::tuple<int32_t, int32_t, int32_t>>(std::make_tuple(*bestFWPosIt, *bestRCPosIt, bestGap));
                                } else {
                                  return (bestGap > dovetailPenalty) ? nonstd::nullopt :
                                  nonstd::optional<std::tuple<int32_t, int32_t, int32_t>>(std::make_tuple(*bestFWPosIt, *bestRCPosIt, bestGap));
                                }
                              };


                              // valid pairings have hits on opposite strands
                              auto& leftFwdHits = (leftIt->fwd) ? leftIt->allPositions : leftIt->oppositeStrandPositions;
                              auto& leftRCHits  = (leftIt->fwd) ? leftIt->oppositeStrandPositions : leftIt->allPositions;

                              auto& rightFwdHits = (rightIt->fwd) ? rightIt->allPositions : rightIt->oppositeStrandPositions;
                              auto& rightRCHits  = (rightIt->fwd) ? rightIt->oppositeStrandPositions : rightIt->allPositions;

                              bool bestHitIsDovetailFWRC{false};
                              bool bestHitIsDovetailRCFW{false};
                              auto bestFWRC = findBestHitFWRC(leftFwdHits, rightRCHits, static_cast<int32_t>(leftIt->readLen), bestHitIsDovetailFWRC);
                              auto bestRCFW = findBestHitFWRC(rightFwdHits, leftRCHits, static_cast<int32_t>(rightIt->readLen), bestHitIsDovetailRCFW);
                              bestMappingIsDovetail = (bestMappingIsDovetail and bestHitIsDovetailFWRC and bestHitIsDovetailRCFW);

                              bool foundValidHit{false};
                              bool leftFwd{false};
                              bool rightFwd{false};
                              int32_t bestGap{std::numeric_limits<int32_t>::max()};
                              int32_t leftPos = -1, rightPos = -1;
                              if (bestFWRC){
                                std::tie(leftPos, rightPos, bestGap) = *bestFWRC;
                                leftFwd = true; rightFwd = false;
                                foundValidHit = true;
                              }
                              if (bestRCFW) {
                                int32_t fwPos, rcPos, bestGapRCFW;
                                std::tie(fwPos, rcPos, bestGapRCFW) = *bestRCFW;
                                if (bestGapRCFW < bestGap) {
                                  leftPos = rcPos;
                                  rightPos = fwPos;
                                  leftFwd = false; rightFwd = true;
                                }
                                foundValidHit = true;
                              }

                              if (foundValidHit) {
                                // If we consider only a single position per transcript
                                int32_t startRead1 = std::max(leftPos, signedZero);
                                int32_t startRead2 = std::max(rightPos, signedZero);
                                bool read1First{(startRead1 < startRead2)};
                                int32_t fragStartPos = read1First ? startRead1 : startRead2;
                                int32_t fragEndPos = read1First ?
                                  (startRead2 + rightIt->readLen) : (startRead1 + leftIt->readLen);
                                uint32_t fragLen = fragEndPos - fragStartPos;
                                jointHits.emplace_back(leftTxp,
                                                       leftPos,
                                                       leftFwd,
                                                       leftIt->readLen,
                                                       fragLen, true);
                                // Fill in the mate info
                                auto& qaln = jointHits.back();
                                qaln.mateLen = rightIt->readLen;
                                qaln.matePos = rightPos;
                                qaln.mateIsFwd = rightFwd;
                                jointHits.back().mateStatus = MateStatus::PAIRED_END_PAIRED;
                                jointHits.back().chainStatus = FragmentChainStatus(leftIt->chainStatus.getLeft(), rightIt->chainStatus.getRight());
                                ++numHits;
                                mergeRes = MergeResult::HAD_CONCORDANT;
                                if (numHits > maxNumHits) { tooManyHits = true; break; }
                              }
                              ++leftIt;

                            } // END if (!(rightTxp < leftTxp))

                            ++rightIt;
                        }
                    }
                    //if (triedHit and jointHits.empty()) { tooManyHits = true;}
                }
                // If we had a potentially valid mapping (hits from different strands), but the
                // (best) only mappings for this fragment were dovetailed, then increment the
                // numDovetails counter here.
                if (bestMappingIsDovetail and hadOppositeStrandMapping) { ++hctr.numDovetails; }

                if (tooManyHits) { jointHits.clear(); ++hctr.tooManyHits; }

                if (mergeRes == MergeResult::HAD_NONE) {
                  // If we had hits on the same transcript, but our status isn't concordant,
                  // then we had discordant hits, otherwise, we had a null intersection.
                  mergeRes = (sameTxpCount > 0) ? MergeResult::HAD_DISCORDANT : MergeResult::HAD_EMPTY_INTERSECTION;
                }
            }

            // If we had proper paired hits
            if (jointHits.size() > 0) {
                hctr.peHits += jointHits.size();
                //orphanStatus = 0;
            }

            return mergeRes;
        }

        inline void mergeLeftRightHits(
                std::vector<QuasiAlignment>& leftHits,
                std::vector<QuasiAlignment>& rightHits,
                std::vector<QuasiAlignment>& jointHits,
                uint32_t readLen,
                uint32_t maxNumHits,
                bool& tooManyHits,
                HitCounters& hctr) {
            if (leftHits.size() > 0) {
                constexpr const int32_t signedZero{0};
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
                                int32_t startRead1 = std::max(leftIt->pos, signedZero);
                                int32_t startRead2 = std::max(rightIt->pos, signedZero);
                                bool read1First{(startRead1 < startRead2)};
                                int32_t fragStartPos = read1First ? startRead1 : startRead2;
                                int32_t fragEndPos = read1First ?
                                    (startRead2 + rightIt->readLen) : (startRead1 + leftIt->readLen);
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
                                jointHits.back().chainStatus = FragmentChainStatus(leftIt->chainStatus.getLeft(), rightIt->chainStatus.getRight());
                                //jointHits.back().completeMatchType = mergeMatchType(leftIt->completeMatchType, rightIt->completeMatchType);

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

    teplate <typename Archive>
    void load(Archive& archive, my_mer& mer);
    */
    }
}


#endif // __RAP_MAP_UTILS_HPP__
