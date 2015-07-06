#ifndef __RAP_MAP_INDEX_HPP__
#define __RAP_MAP_INDEX_HPP__

#include <fstream>

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "ScopedTimer.hpp"
#include "RapMapUtils.hpp"

class RapMapIndex {
    using TranscriptList = std::vector<uint32_t>;
    using PositionList = std::vector<uint32_t>;
    using KmerIndex = std::unordered_map<uint64_t, TranscriptList, rapmap::utils::KmerKeyHasher>;
    using IntervalIndex = std::unordered_map<uint64_t, rapmap::utils::KmerInterval, rapmap::utils::KmerKeyHasher>;

    public:
        RapMapIndex() {}

        bool load(std::string& indexPrefix) {

            std::string intervalIdxName = indexPrefix + ".imap.bin";
            std::string tidListName = indexPrefix + ".tid.bin";
            std::string posListName = indexPrefix + ".pos.bin";
            std::string txpNameFile = indexPrefix + ".txps.bin";

            std::ifstream imapStream(intervalIdxName, std::ios::binary);
            {
                std::cerr << "loading interval index ";
                ScopedTimer timer;
                cereal::BinaryInputArchive idxArchive(imapStream);
                idxArchive(idx);
                std::cerr << "done ";
            }
            imapStream.close();
            std::ifstream tidStream(tidListName, std::ios::binary); 
            {
                std::cerr << "loading transcript list ";
                ScopedTimer timer;
                cereal::BinaryInputArchive tidArchive(tidStream);
                tidArchive(tidList);
                std::cerr << "done ";
            }
            tidStream.close();
            std::ifstream posStream(posListName, std::ios::binary); 
            {
                std::cerr << "loading position list ";
                ScopedTimer timer;
                cereal::BinaryInputArchive posArchive(posStream);
                posArchive(posList);
                std::cerr << "done ";
            }
            posStream.close();

            std::ifstream txpNameStream(txpNameFile, std::ios::binary); 
            {
                std::cerr << "loading position list ";
                ScopedTimer timer;
                cereal::BinaryInputArchive txpNameArchive(txpNameStream);
                txpNameArchive(txpNames);
                std::cerr << "done ";
            }
            txpNameStream.close();
            return true;
        }

    public:
    IntervalIndex idx;
    TranscriptList tidList;
    PositionList posList;
    std::vector<std::string> txpNames;
};
 
#endif //__RAP_MAP_INDEX_HPP__
