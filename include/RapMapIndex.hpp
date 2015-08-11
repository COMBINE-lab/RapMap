#ifndef __RAP_MAP_INDEX_HPP__
#define __RAP_MAP_INDEX_HPP__

#include <fstream>
#include <memory>

//#include "jellyfish/jellyfish.hpp"
#include "jellyfish/file_header.hpp"
#include "jellyfish/binary_dumper.hpp"
#include "jellyfish/hash_counter.hpp"
#include "jellyfish/mapped_file.hpp"
#include "JFRaw.hpp"

#include "spdlog/spdlog.h"

#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>

#include "RapMapUtils.hpp"
#include "ScopedTimer.hpp"

class RapMapIndex {
    using PositionList = std::vector<uint32_t>;
    using KmerInfoList = std::vector<rapmap::utils::KmerInfo>;
    using EqClassList = std::vector<rapmap::utils::EqClass>;
    //using MerMapT = jellyfish::cooperative::hash_counter<rapmap::utils::my_mer>;
    using FileMerArray = jellyfish::large_hash::array_raw<rapmap::utils::my_mer>;
    using EqClassLabelVec = std::vector<uint32_t>;

    //using KmerIndex = std::unordered_map<uint64_t, TranscriptList, rapmap::utils::KmerKeyHasher>;
    //using IntervalIndex = std::unordered_map<uint64_t, rapmap::utils::KmerInterval, rapmap::utils::KmerKeyHasher>;

    public:
    RapMapIndex();

    bool load(std::string& indexPrefix);

    KmerInfoList kmerInfos;
    std::unique_ptr<char> rawHashMem{nullptr};
    std::unique_ptr<FileMerArray> merHash{nullptr};
    EqClassList eqClassList;
    EqClassLabelVec eqLabelList;
    PositionList posList;
    std::vector<std::string> txpNames;
    std::vector<uint32_t> txpLens;
    std::vector<uint8_t> fwdJumpTable;
    std::vector<uint8_t> revJumpTable;
};

#endif //__RAP_MAP_INDEX_HPP__
