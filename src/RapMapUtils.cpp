#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/archives/binary.hpp>

namespace rapmap {
    namespace utils {
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
