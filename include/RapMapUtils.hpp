#ifndef __RAP_MAP_UTILS_HPP__
#define __RAP_MAP_UTILS_HPP__

#include <cereal/archives/binary.hpp>
#include "jellyfish/mer_dna.hpp"
namespace rapmap {
    namespace utils {

    using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

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
    /*
    template <typename Archive>
    void save(Archive& archive, const my_mer& mer);

    template <typename Archive>
    void load(Archive& archive, my_mer& mer);
    */
    }
}
#endif // __RAP_MAP_UTILS_HPP__
