#ifndef __RAPMAP_SA_INDEX_HPP__
#define __RAPMAP_SA_INDEX_HPP__

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"
#include "RSDic.hpp"
#include "google/dense_hash_map"
#include "bit_array.h"
//#include "bitmap.h"
//#include "shared.h"
#include "rank9b.h"


#include <cstdio>
#include <vector>
#include <memory>

class RapMapSAIndex {
    public:
	struct BitArrayDeleter {
	   void operator()(BIT_ARRAY* b) {
	       if(b != nullptr) {
	           bit_array_free(b);
	       }
	   }
	};

	using BitArrayPointer = std::unique_ptr<BIT_ARRAY, BitArrayDeleter>;

        RapMapSAIndex() {}

	// Given a position, p, in the concatenated text,
	// return the corresponding transcript
	uint32_t transcriptAtPosition(uint32_t p) {
	   /*
	    auto rsafe = rankDictSafe.Rank(p, 1);
	    auto r = rankDict->rank(p);
	    if (r != rsafe) {
		std::cerr << "RANK VECTOR IMPLEMENTATIONS DISAGREED!\n";
		std::cerr << "RSDic said " << rsafe << ", uncompressed said " << r << '\n';
	    }
	    */
	    return rankDict->rank(p);
	}

        bool load(const std::string& indDir) {

            auto logger = spdlog::get("stderrLog");
            size_t n{0};
            std::ifstream saStream(indDir + "sa.bin");
            {
                logger->info("Loading Suffix Array ");
                cereal::BinaryInputArchive saArchive(saStream);
                saArchive(SA);
                //saArchive(LCP);
            }
            saStream.close();

            khash.set_empty_key(std::numeric_limits<uint64_t>::max());
            uint32_t k;
            std::ifstream hashStream(indDir + "hash.bin");
            {
                logger->info("Loading Position Hash");
                cereal::BinaryInputArchive hashArchive(hashStream);
                hashArchive(k);
                khash.unserialize(google::dense_hash_map<uint64_t,
                        rapmap::utils::SAInterval,
                        rapmap::utils::KmerKeyHasher>::NopointerSerializer(), &hashStream);

                //hashArchive(khash);
            }
            hashStream.close();
            rapmap::utils::my_mer::k(k);

            std::ifstream seqStream(indDir + "txpInfo.bin");
            {
                logger->info("Loading Transcript Info ");
                cereal::BinaryInputArchive seqArchive(seqStream);
                seqArchive(txpNames);
                seqArchive(txpOffsets);
                //seqArchive(positionIDs);
                seqArchive(seq);
            }
            seqStream.close();

	    /*
	    std::ifstream rsStream(indDir + "rsdSafe.bin", std::ios::binary);
	    {
		    logger->info("Loading Rank-Select Data");
		    rankDictSafe.Load(rsStream);
	    }
	    rsStream.close();
	    */
	    std::string rsFileName = indDir + "rsd.bin";
	    FILE* rsFile = fopen(rsFileName.c_str(), "r");
	    {
		logger->info("Loading Rank-Select Bit Array");
		bitArray.reset(bit_array_create(0));
		if (!bit_array_load(bitArray.get(), rsFile)) {
			logger->error("Couldn't load bit array from {}!", rsFileName);
			std::exit(1);
		}
		logger->info("There were {} set bits in the bit array", bit_array_num_bits_set(bitArray.get()));
		rankDict.reset(new rank9b(bitArray->words, bitArray->num_of_bits));
	    }
	    fclose(rsFile);

            {
                logger->info("Computing transcript lengths");
                txpLens.resize(txpOffsets.size());
                if (txpOffsets.size() > 1) {
                    for(size_t i = 0; i < txpOffsets.size() - 1; ++i) {
                        auto nextOffset = txpOffsets[i+1];
                        auto currentOffset = txpOffsets[i];
                        txpLens[i] = (nextOffset - 1) - currentOffset;
                    }
                }
                // The last length is just the length of the suffix array - the last offset
                txpLens[txpOffsets.size()-1] = (SA.size() - 1) - txpOffsets[txpOffsets.size() - 1];
            }

            logger->info("Done loading index");
            return true;
        }

    std::vector<int> SA;
    //rsdic::RSDic rankDictSafe;
	
    BitArrayPointer bitArray{nullptr};
    std::unique_ptr<rank9b> rankDict{nullptr};

    std::string seq;
    std::vector<std::string> txpNames;
    std::vector<uint32_t> txpOffsets;
    std::vector<uint32_t> txpLens;
    std::vector<uint32_t> positionIDs;
    google::dense_hash_map<uint64_t,
                        rapmap::utils::SAInterval,
                        rapmap::utils::KmerKeyHasher> khash;
                       // ::NopointerSerializer(), &hashStream);
        /*
    std::unordered_map<uint64_t,
                       rapmap::utils::SAInterval,
                       rapmap::utils::KmerKeyHasher> khash;
                       */
};
#endif //__RAPMAP_SA_INDEX_HPP__
