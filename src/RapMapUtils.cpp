#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/archives/binary.hpp>
#include "RapMapUtils.hpp"

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

            readWork.resize(seq.length());
            qualWork.resize(qual.length());
            int32_t end = seq.length(), start = 0;
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
            std::swap(seq, readWork);
            std::swap(qual, qualWork);
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
