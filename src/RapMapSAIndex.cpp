#include "RapMapSAIndex.hpp"

RapMapSAIndex::RapMapSAIndex() {}

// Given a position, p, in the concatenated text,
// return the corresponding transcript
uint32_t RapMapSAIndex::transcriptAtPosition(uint32_t p) {
    return rankDict->rank(p);
}

bool RapMapSAIndex::load(const std::string& indDir) {

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

