#include "RapMapSAIndex.hpp"
#include "ScopedTimer.hpp"
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


// Given a key (k-mer), return a pointer to the SA
// interval to which it corresponds.
// NOTE: Returns nullptr if the key is not found.
// template <typename IndexT>
// inline rapmap::utils::SAInterval<IndexT>*
// RapMapSAIndex<IndexT>::intervalForKmer(rapmap::utils::my_mer& key);

template <typename IndexT>
bool RapMapSAIndex<IndexT>::load(const std::string& indDir) {

    auto logger = spdlog::get("stderrLog");
    size_t n{0};

    // This part takes the longest, so do it in it's own asynchronous task
    std::future<std::pair<bool, uint32_t>> loadingHash = std::async(std::launch::async, [this, logger, indDir]() -> std::pair<bool, uint32_t> {
            std::string jfFileName = indDir + "kmers.jfhash";
            std::string saIntervalFileName = indDir + "sa_intervals.bin";

            // Load the kmer info list first --- this will
            // give us the # of unique k-mers
            std::ifstream kmerInfoStream(saIntervalFileName, std::ios::binary);
            {
            logger->info("loading k-mer info list . . .");
            ScopedTimer timer;
            cereal::BinaryInputArchive kmerInfoArchive(kmerInfoStream);
            kmerInfoArchive(this->saIntervals);
            logger->info("done\n");
            }
            kmerInfoStream.close();

            // Now load the jellyfish hash
            size_t numDistinctKmers = saIntervals.size();
            uint32_t k{0};
            {
            ScopedTimer timer;
            logger->info("loading k-mer => id hash . . . ");
            std::ifstream bis(jfFileName);
            const SpecialHeader bh(bis);
            // mapFile.reset(new jellyfish::mapped_file(jfFileName.c_str()));
            const size_t sizeInBytes = bh.size_bytes();
            // Load the hash from file
            logger->info("\theader format = {}"      , bh.format());
            logger->info("\t# distinct k-mers = {}"  , numDistinctKmers);
            logger->info("\thash key len = {}"       , bh.key_len());
            logger->info("\tcounter len = {}"        , bh.counter_len());
            logger->info("\tmax reprobe offset = {}" , bh.max_reprobe());
            logger->info("\tsize in bytes = {}"      , sizeInBytes);

            // Allocate the actual storage
            rawHashMem.reset(new char[sizeInBytes]);
            bis.read(this->rawHashMem.get(), sizeInBytes);
            // We can close the file now
            bis.close();

            khash.reset( new FileMerArray(this->rawHashMem.get(),//mapFile->base() + bh.offset(),
                        sizeInBytes,
                        bh.size(),
                        bh.key_len(),
                        bh.counter_len(),
                        bh.max_reprobe(),
                        bh.matrix()));
            // Set the key size
            rapmap::utils::my_mer::k(bh.key_len() / 2);
            k = bh.key_len() / 2;
            logger->info("done");
            }
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
