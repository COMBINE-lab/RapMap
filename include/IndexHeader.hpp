#ifndef __INDEX_HEADER_HPP__
#define __INDEX_HEADER_HPP__

#include <cereal/types/string.hpp>

// The different types of indices supported
enum class IndexType : uint8_t {
    PSEUDO = 0,
    QUASI,
    INVALID
};

class IndexHeader {
    public:
        IndexHeader () : type_(IndexType::INVALID), versionString_("invalid"), usesKmers_(false), kmerLen_(0) {}

        IndexHeader(IndexType typeIn, const std::string& versionStringIn,
                    bool usesKmersIn, uint32_t kmerLenIn, bool bigSA = false):
                    type_(typeIn), versionString_(versionStringIn),
                    usesKmers_(usesKmersIn), kmerLen_(kmerLenIn), bigSA_(bigSA) {}

        template <typename Archive>
            void save(Archive& ar) const {
                ar( cereal::make_nvp("IndexType", type_) );
                ar( cereal::make_nvp("IndexVersion", versionString_) );
                ar( cereal::make_nvp("UsesKmers", usesKmers_) );
                ar( cereal::make_nvp("KmerLen", kmerLen_) );
                ar( cereal::make_nvp("BigSA", bigSA_) );
            }

        template <typename Archive>
            void load(Archive& ar) {
                ar( cereal::make_nvp("IndexType", type_) );
                ar( cereal::make_nvp("IndexVersion", versionString_) );
                ar( cereal::make_nvp("UsesKmers", usesKmers_) );
                ar( cereal::make_nvp("KmerLen", kmerLen_) );
                ar( cereal::make_nvp("BigSA", bigSA_) );
            }

        IndexType indexType() const { return type_; }
        std::string version() const { return versionString_; }
        bool usesKmers() const { return usesKmers_; }
        uint32_t kmerLen() const { return kmerLen_; }
        bool bigSA() const { return bigSA_; }

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
};


#endif // __INDEX_HEADER_HPP__
