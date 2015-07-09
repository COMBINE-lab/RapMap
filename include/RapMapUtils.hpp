#ifndef __RAP_MAP_UTILS_HPP__
#define __RAP_MAP_UTILS_HPP__

#include <cereal/archives/binary.hpp>
#include "jellyfish/mer_dna.hpp"
namespace rapmap {
    namespace utils {

    using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

    constexpr uint32_t newTxpSetMask = 0x80000000;
    constexpr uint32_t rcSetMask = 0x40000000;

    // positions are stored in a packed format, where the highest
    // 2-bits encode if this position refers to a new transcript
    // and whether or not the k-mer from the hash matches this txp
    // in the forward or rc direction.
    void decodePosition(uint32_t p, int32_t& pOut, bool& newTxp, bool& isRC);

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

    template <class T>
    inline void hashCombine(std::size_t& seed, const T& v)
    {
            std::hash<T> hasher;
            seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }

    constexpr uint32_t uint32Invalid = std::numeric_limits<uint32_t>::max();

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


    /*
    template <typename Archive>
    void save(Archive& archive, const my_mer& mer);

    template <typename Archive>
    void load(Archive& archive, my_mer& mer);
    */
    }
}
#endif // __RAP_MAP_UTILS_HPP__
