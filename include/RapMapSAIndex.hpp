#ifndef __RAPMAP_SA_INDEX_HPP__
#define __RAPMAP_SA_INDEX_HPP__

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/details/format.h"
//#include "RSDic.hpp"
#include "google/dense_hash_map"

#include <vector>

class RapMapSAIndex {
    public:
        RapMapSAIndex() {}

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
                seqArchive(positionIDs);
                seqArchive(seq);
            }
            seqStream.close();

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
    //std::vector<int> LCP;
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
