#include "BooMap.hpp"
#include "RapMapSAIndex.hpp"
#include "IndexHeader.hpp"
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>


#include <future>
#include <thread>

// These are **free** functions that are used for loading the
// appropriate type of hash.
template <typename IndexT>
bool loadHashFromIndex(const std::string& indexDir,
                       google::dense_hash_map<uint64_t,
                       rapmap::utils::SAInterval<IndexT>,
                       rapmap::utils::KmerKeyHasher>& khash) {
      khash.set_empty_key(std::numeric_limits<uint64_t>::max());
      std::ifstream hashStream(indexDir + "hash.bin");
      khash.unserialize(typename google::dense_hash_map<uint64_t,
                      rapmap::utils::SAInterval<IndexT>,
                      rapmap::utils::KmerKeyHasher>::NopointerSerializer(), &hashStream);
      return true;
}

template <typename IndexT>
bool loadHashFromIndex(const std::string& indexDir,
		       BooMap<uint64_t, rapmap::utils::SAInterval<IndexT>> & h) {
    std::string hashBase = indexDir + "hash_info";
    h.load(hashBase);
    return true;
}

template <typename IndexT, typename HashT>
RapMapSAIndex<IndexT, HashT>::RapMapSAIndex() {}

// Given a position, p, in the concatenated text,
// return the corresponding transcript
template <typename IndexT, typename HashT>
IndexT RapMapSAIndex<IndexT, HashT>::transcriptAtPosition(IndexT p) {
    return rankDict->rank(p);
}

template <typename IndexT, typename HashT>
bool RapMapSAIndex<IndexT, HashT>::load(const std::string& indDir) {

    auto logger = spdlog::get("stderrLog");
    size_t n{0};

    IndexHeader h;
    std::ifstream indexStream(indDir + "header.json");
    {
      cereal::JSONInputArchive ar(indexStream);
      ar(h);
    }
    indexStream.close();
    uint32_t idxK = h.kmerLen();

    // This part takes the longest, so do it in it's own asynchronous task
    std::future<bool> loadingHash = std::async(std::launch::async, [this, logger, indDir]() -> bool {
	   if (loadHashFromIndex(indDir, khash)) {
                logger->info("Successfully loaded position hash");
                return true;
            } else {
                logger->error("Failed to load position hash!");
                return false;
            }
	// If using a google dense hash
        //this->khash.set_empty_key(std::numeric_limits<uint64_t>::max());
        //uint32_t k = 31;
        //std::ifstream hashStream(indDir + "hash.bin");
        //{

	  //logger->info("Loading Position Hash");
            //khash.load(hashStream);
            //cereal::BinaryInputArchive hashArchive(hashStream);
            //hashArchive(k);
            //khash.unserialize(typename google::dense_hash_map<uint64_t,
            //        rapmap::utils::SAInterval<IndexT>,
            //        rapmap::utils::KmerKeyHasher>::NopointerSerializer(), &hashStream);
            //hashArchive(khash);
	   //}
        //hashStream.close();
        //std::cerr << "had " << khash.size() << " entries\n";
        //return true;
    });

    /*
    std::ifstream intervalStream(indDir + "kintervals.bin");
    {
        logger->info("Loading k-mer intervals");
        cereal::BinaryInputArchive intervalArchive(intervalStream);
        intervalArchive(kintervals);
    }
    intervalStream.close();
    */

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
    if (!hashLoadRes) {
        logger->error("Failed to load hash!");
        std::exit(1);
    }
    rapmap::utils::my_mer::k(idxK);

    logger->info("Done loading index");
    return true;
}

template class RapMapSAIndex<int32_t,  google::dense_hash_map<uint64_t,
                      rapmap::utils::SAInterval<int32_t>,
                      rapmap::utils::KmerKeyHasher>>;
template class RapMapSAIndex<int64_t,  google::dense_hash_map<uint64_t,
                      rapmap::utils::SAInterval<int64_t>,
                      rapmap::utils::KmerKeyHasher>>;
template class RapMapSAIndex<int32_t, BooMap<uint64_t, rapmap::utils::SAInterval<int32_t>>>;
template class RapMapSAIndex<int64_t, BooMap<uint64_t, rapmap::utils::SAInterval<int64_t>>>;
