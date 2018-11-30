//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __INDEX_HEADER_HPP__
#define __INDEX_HEADER_HPP__

#include "spdlog/spdlog.h"
#include <cereal/types/string.hpp>

// The different types of indices supported
enum class IndexType : uint8_t {
    PSEUDO = 0,
    QUASI,
    INVALID
};

class IndexHeader {
    public:
        IndexHeader () : type_(IndexType::INVALID), versionString_("invalid"), usesKmers_(false), kmerLen_(0), perfectHash_(false) {}

        IndexHeader(IndexType typeIn, const std::string& versionStringIn,
                    bool usesKmersIn, uint32_t kmerLenIn, bool bigSA = false, bool perfectHash = false):
                    type_(typeIn), versionString_(versionStringIn),
                    usesKmers_(usesKmersIn), kmerLen_(kmerLenIn), bigSA_(bigSA),
                    perfectHash_(perfectHash) {}

        template <typename Archive>
            void save(Archive& ar) const {
                ar( cereal::make_nvp("IndexType", type_) );
                ar( cereal::make_nvp("IndexVersion", versionString_) );
                ar( cereal::make_nvp("UsesKmers", usesKmers_) );
                ar( cereal::make_nvp("KmerLen", kmerLen_) );
                ar( cereal::make_nvp("BigSA", bigSA_) );
                ar( cereal::make_nvp("PerfectHash", perfectHash_) );
                ar( cereal::make_nvp("SeqHash", seqHash256_) );
                ar( cereal::make_nvp("NameHash", nameHash256_) );
                ar( cereal::make_nvp("SeqHash512", seqHash512_) );
                ar( cereal::make_nvp("NameHash512", nameHash512_) );
            }

        template <typename Archive>
        void load(Archive& ar) {
            try {
                ar( cereal::make_nvp("IndexType", type_) );
                ar( cereal::make_nvp("IndexVersion", versionString_) );
                ar( cereal::make_nvp("UsesKmers", usesKmers_) );
                ar( cereal::make_nvp("KmerLen", kmerLen_) );
                ar( cereal::make_nvp("BigSA", bigSA_) );
                ar( cereal::make_nvp("PerfectHash", perfectHash_) );
                ar( cereal::make_nvp("SeqHash", seqHash256_) );
                ar( cereal::make_nvp("NameHash", nameHash256_) );
                ar( cereal::make_nvp("SeqHash512", seqHash512_) );
                ar( cereal::make_nvp("NameHash512", nameHash512_) );
            } catch (const cereal::Exception& e) {
                auto cerrLog = spdlog::get("stderrLog");
                cerrLog->error("Encountered exception [{}] when loading index.", e.what());
                cerrLog->error("The index was likely build with an older (and incompatible) "
                               "version of RapMap.  Please re-build the index with a compatible version.");
                cerrLog->flush(); 
                std::exit(1);
            }
        }

        IndexType indexType() const { return type_; }
        std::string version() const { return versionString_; }
        bool usesKmers() const { return usesKmers_; }
        uint32_t kmerLen() const { return kmerLen_; }
        bool bigSA() const { return bigSA_; }
        bool perfectHash() const { return perfectHash_; }

        void setSeqHash256(const std::string& seqHash) { seqHash256_ = seqHash; }
        void setNameHash256(const std::string& nameHash) { nameHash256_ = nameHash; }
        std::string seqHash256() const { return seqHash256_; }
        std::string nameHash256() const { return nameHash256_; }

        void setSeqHash512(const std::string& seqHash) { seqHash512_ = seqHash; }
        void setNameHash512(const std::string& nameHash) { nameHash512_ = nameHash; }
        std::string seqHash512() const { return seqHash512_; }
        std::string nameHash512() const { return nameHash512_; }
    private:
        // The type of index we have
        IndexType type_;
        // The version string for the index
        std::string versionString_;
        // True if this index makes use of k-mers false otherwise
        // (currently, all supported indices use k-mers in some form)
        bool usesKmers_;
        // The length of k-mer used by the index
        uint32_t kmerLen_;
        // Do we have a 64-bit suffix array or not
        bool bigSA_;
        // Are we using a perfect hash in the index or not?
        bool perfectHash_;
        // Hash of sequences in txome
        std::string seqHash256_;
        // Hash of sequence names in txome
        std::string nameHash256_;
        // Hash of sequences in txome
        std::string seqHash512_;
        // Hash of sequence names in txome
        std::string nameHash512_;
};


#endif // __INDEX_HEADER_HPP__
