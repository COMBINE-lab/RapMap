#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/archives/binary.hpp>

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
