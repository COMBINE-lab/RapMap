#include "RapMapIndex.hpp"

RapMapIndex::RapMapIndex() {}

bool RapMapIndex::load(std::string& indexPrefix) {
    auto logger = spdlog::get("stderrLog");
    std::string kmerInfosName = indexPrefix + "kinfo.bin";
    std::string eqClassListName = indexPrefix + "eqclass.bin";
    std::string eqLabelListName = indexPrefix + "eqlab.bin";
    std::string posListName = indexPrefix + "pos.bin";
    std::string jfFileName = indexPrefix + "rapidx.jfhash";
    std::string txpNameFile = indexPrefix + "txpnames.bin";
    std::string txpLenFile = indexPrefix + "txplens.bin";
    std::string fwdJumpFile = indexPrefix + "fwdjump.bin";
    std::string revJumpFile = indexPrefix + "revjump.bin";

    // Load the kmer info list first --- this will
    // give us the # of unique k-mers
    std::ifstream kmerInfoStream(kmerInfosName, std::ios::binary);
    {
        logger->info("loading k-mer info list . . .");
        ScopedTimer timer;
        cereal::BinaryInputArchive kmerInfoArchive(kmerInfoStream);
        kmerInfoArchive(kmerInfos);
        logger->info("done\n");
    }
    kmerInfoStream.close();

    size_t numDistinctKmers = kmerInfos.size();
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
        bis.read(rawHashMem.get(), sizeInBytes);
        // We can close the file now
        bis.close();

        merHash.reset( new FileMerArray(rawHashMem.get(),//mapFile->base() + bh.offset(),
                    sizeInBytes,
                    bh.size(),
                    bh.key_len(),
                    bh.counter_len(),
                    bh.max_reprobe(),
                    bh.matrix()));
        // Set the key size
        rapmap::utils::my_mer::k(bh.key_len() / 2);
        logger->info("done");
    }


    std::ifstream eqClassStream(eqClassListName, std::ios::binary);
    {
        logger->info("loading eq classes . . . ");
        ScopedTimer timer;
        cereal::BinaryInputArchive eqClassArchive(eqClassStream);
        eqClassArchive(eqClassList);
        logger->info("[{}] classes", eqClassList.size());
        logger->info("done");
    }
    eqClassStream.close();
    std::ifstream eqLabelStream(eqLabelListName, std::ios::binary);
    {
        logger->info("loading eq class labels . . . ");
        ScopedTimer timer;
        cereal::BinaryInputArchive eqLabelArchive(eqLabelStream);
        eqLabelArchive(eqLabelList);
        logger->info("[{}] labels", eqLabelList.size());
        logger->info("done");
    }
    eqLabelStream.close();
    std::ifstream posStream(posListName, std::ios::binary);
    {
        logger->info("loading position list . . . ");
        ScopedTimer timer;
        cereal::BinaryInputArchive posArchive(posStream);
        posArchive(posList);
        logger->info("[{}] total k-mer positions", posList.size());
        logger->info("done");
    }
    posStream.close();
    std::ifstream txpNameStream(txpNameFile, std::ios::binary);
    {
        logger->info("loading transcript names ");
        ScopedTimer timer;
        cereal::BinaryInputArchive txpNameArchive(txpNameStream);
        txpNameArchive(txpNames);
        logger->info("[{}] transcripts in index ", txpNames.size());
        logger->info("done ");
    }
    txpNameStream.close();

    std::ifstream txpLenStream(txpLenFile, std::ios::binary);
    {
        logger->info("loading transcript lengths");
        ScopedTimer timer;
        cereal::BinaryInputArchive txpLenArchive(txpLenStream);
        txpLenArchive(txpLens);
        logger->info("[{}] transcripts in index", txpLens.size());
        logger->info("done ");
    }
    txpLenStream.close();

    std::ifstream fwdJumpStream(fwdJumpFile, std::ios::binary);
    {
        logger->info("loading forward jumps");
        ScopedTimer timer;
        cereal::BinaryInputArchive fwdJumpArchive(fwdJumpStream);
        fwdJumpArchive(fwdJumpTable);
        logger->info("[{}] forward jumps", fwdJumpTable.size());
        logger->info("done ");
    }
    fwdJumpStream.close();

    std::ifstream revJumpStream(revJumpFile, std::ios::binary);
    {
        logger->info("loading reverse jumps");
        ScopedTimer timer;
        cereal::BinaryInputArchive revJumpArchive(revJumpStream);
        revJumpArchive(revJumpTable);
        logger->info("[{}] reverse jumps", revJumpTable.size());
        logger->info("done ");
    }
    revJumpStream.close();
    return true;
}

