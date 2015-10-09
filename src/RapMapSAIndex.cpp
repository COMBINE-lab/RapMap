#include "RapMapSAIndex.hpp"
#include <future>
#include <thread>

template <typename IndexT>
RapMapSAIndex<IndexT>::RapMapSAIndex() {}

// Given a position, p, in the concatenated text,
// return the corresponding transcript
template <typename IndexT>
IndexT RapMapSAIndex<IndexT>::transcriptAtPosition(IndexT p) {
    return rankDict->rank(p);
}

template <typename IndexT>
bool RapMapSAIndex<IndexT>::load(const std::string& indDir) {

    auto logger = spdlog::get("stderrLog");
    size_t n{0};

    // This part takes the longest, so do it in it's own asynchronous task
    std::future<std::pair<bool, uint32_t>> loadingHash = std::async(std::launch::async, [this, logger, indDir]() -> std::pair<bool, uint32_t> {
        this->khash.set_empty_key(std::numeric_limits<uint64_t>::max());
        uint32_t k;
        std::ifstream hashStream(indDir + "hash.bin");
        {
            logger->info("Loading Position Hash");
            cereal::BinaryInputArchive hashArchive(hashStream);
            hashArchive(k);
            khash.unserialize(typename google::dense_hash_map<uint64_t,
                    rapmap::utils::SAInterval<IndexT>,
                    rapmap::utils::KmerKeyHasher>::NopointerSerializer(), &hashStream);

            //hashArchive(khash);
        }
        hashStream.close();
        return std::make_pair(true, k);
    });



    std::ifstream saStream(indDir + "sa.bin");
    {
        logger->info("Loading Suffix Array ");
        cereal::BinaryInputArchive saArchive(saStream);
        saArchive(SA);
        //saArchive(LCP);
    }
    saStream.close();

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

    logger->info("Waiting to finish loading hash");
    loadingHash.wait();
    auto hashLoadRes = loadingHash.get();
    if (!hashLoadRes.first) {
        logger->error("Failed to load hash!");
        std::exit(1);
    }
    rapmap::utils::my_mer::k(hashLoadRes.second);



    logger->info("Done loading index");
    return true;
}

template class RapMapSAIndex<int32_t>;
template class RapMapSAIndex<int64_t>;
