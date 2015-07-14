#ifndef __RAP_MAP_UTILS_HPP__
#define __RAP_MAP_UTILS_HPP__

#include "xxhash.h"
#include <cereal/archives/binary.hpp>
#include "jellyfish/mer_dna.hpp"
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



    /*
    template <typename Archive>
    void save(Archive& archive, const my_mer& mer);

    template <typename Archive>
    void load(Archive& archive, my_mer& mer);
    */
    }
}
#endif // __RAP_MAP_UTILS_HPP__
