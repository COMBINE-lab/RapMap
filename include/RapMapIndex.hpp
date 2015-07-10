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
        RapMapIndex() {}

        bool load(std::string& indexPrefix) {

            std::string kmerInfosName = indexPrefix + "kinfo.bin";
            std::string eqClassListName = indexPrefix + "eqclass.bin";
            std::string eqLabelListName = indexPrefix + "eqlab.bin";
            std::string posListName = indexPrefix + "pos.bin";
            std::string jfFileName = indexPrefix + "rapidx.jfhash";
            std::string txpNameFile = indexPrefix + "txpnames.bin";

            // Load the kmer info list first --- this will
            // give us the # of unique k-mers
            std::ifstream kmerInfoStream(kmerInfosName, std::ios::binary);
            {
                std::cerr << "loading k-mer info list . . .";
                ScopedTimer timer;
                cereal::BinaryInputArchive kmerInfoArchive(kmerInfoStream);
                kmerInfoArchive(kmerInfos);
                std::cerr << "done\n";
            }
            kmerInfoStream.close();

            size_t numDistinctKmers = kmerInfos.size();

            {
            std::ifstream bis(jfFileName);
            const SpecialHeader bh(bis);
            mapFile.reset(new jellyfish::mapped_file(jfFileName.c_str()));
            const size_t sizeInBytes = bh.size_bytes();

            // Load the hash from file
            std::cerr << "Header format = " << bh.format() << "\n";
            std::cerr << "# distinct k-mers = " << numDistinctKmers << "\n";
            std::cerr << "jf hash key len = " << bh.key_len() << "\n";
            std::cerr << "jf counter len = " << bh.counter_len() << "\n";
            std::cerr << "jf max reprobe offset = " << bh.max_reprobe() << "\n";
            std::cerr << "Size in bytes = " << sizeInBytes << "\n";
            std::cerr << std::endl;

            merHash.reset( new FileMerArray(mapFile->base() + bh.offset(),
                                            sizeInBytes,
                                            bh.size(),
                                            bh.key_len(),
                                            bh.counter_len(),
                                            bh.max_reprobe(),
                                            bh.matrix()));
            // Set the key size
            rapmap::utils::my_mer::k(bh.key_len() / 2);

            /*
            // Yes --- this is a leak.  For some reason we get a segfault
            // when this goes away?!?
            jellyfish::binary_reader<rapmap::utils::my_mer, uint32_t>* br =
                new jellyfish::binary_reader<rapmap::utils::my_mer, uint32_t>(bis, &bh);

            size_t hashSize = 2 << (static_cast<uint32_t>(
                                    std::round(std::log2(numDistinctKmers)) + 1));
            merHash.reset( new MerMapT(hashSize, bh.key_len(),
                             bh.counter_len(), 1, 126));

            // Iterate through the reader, and copy the
            // k-mers over to the new merHash.
            size_t i = 0;
            while (br->next()) {
                merHash->add(br->key(), br->val());
                if (i % 1000000 == 0) {
                    std::cerr << "\r\rpopulated " << i << " k-mers";
                }
                ++i;
            }
            std::cerr << "\n";
            */
            //std::swap(kmerHash, merHash);
            //bis.close();
            // Done loading JF hash?
            }


            std::ifstream eqClassStream(eqClassListName, std::ios::binary);
            {
                std::cerr << "loading eq classes . . . ";
                ScopedTimer timer;
                cereal::BinaryInputArchive eqClassArchive(eqClassStream);
                eqClassArchive(eqClassList);
                std::cerr << "done\n";
            }
            eqClassStream.close();
            std::ifstream eqLabelStream(eqLabelListName, std::ios::binary);
            {
                std::cerr << "loading eq class labels . . . ";
                ScopedTimer timer;
                cereal::BinaryInputArchive eqLabelArchive(eqLabelStream);
                eqLabelArchive(eqLabelList);
                std::cerr << "done\n";
            }
            eqLabelStream.close();
            std::ifstream posStream(posListName, std::ios::binary);
            {
                std::cerr << "loading position list . . . ";
                ScopedTimer timer;
                cereal::BinaryInputArchive posArchive(posStream);
                posArchive(posList);
                std::cerr << "done\n";
            }
            posStream.close();
            std::ifstream txpNameStream(txpNameFile, std::ios::binary);
            {
                std::cerr << "loading transcript names ";
                ScopedTimer timer;
                cereal::BinaryInputArchive txpNameArchive(txpNameStream);
                txpNameArchive(txpNames);
                std::cerr << "done ";
            }
            txpNameStream.close();
            return true;
        }

    public:
    KmerInfoList kmerInfos;
    std::unique_ptr<jellyfish::mapped_file> mapFile{nullptr};
    std::unique_ptr<FileMerArray> merHash{nullptr};
    EqClassList eqClassList;
    EqClassLabelVec eqLabelList;
    PositionList posList;
    std::vector<std::string> txpNames;
};

#endif //__RAP_MAP_INDEX_HPP__
