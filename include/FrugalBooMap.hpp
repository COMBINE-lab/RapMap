#ifndef __BOO_MAP_FRUGAL__
#define __BOO_MAP_FRUGAL__

#include "BooPHF.hpp"
#include "RapMapUtils.hpp"

#include "cereal/types/vector.hpp"
#include "cereal/types/utility.hpp"
#include "cereal/archives/binary.hpp"

#include <fstream>
#include <vector>
#include <iterator>
#include <type_traits>

#include <sys/stat.h>

// adapted from :
// http://stackoverflow.com/questions/34875315/implementation-my-own-list-and-iterator-stl-c
template <typename Iter, typename IndexT>
class KeyProxyIterator {
public:
    typedef KeyProxyIterator<Iter, IndexT> self_type;
    typedef uint64_t value_type;//std::iterator_traits<Iter>::value_type::first_type value_type;
    typedef value_type& reference;
    typedef value_type* pointer;
    typedef std::forward_iterator_tag iterator_category;
    typedef int64_t difference_type;

    KeyProxyIterator(Iter first, std::vector<IndexT>* saPtr, const char* txtPtr, unsigned int k) : 
        curr_(first), saPtr_(saPtr), txtPtr_(txtPtr) {}
    KeyProxyIterator operator++() { KeyProxyIterator i = *this; curr_++; return i; }
    KeyProxyIterator operator++(int) { ++curr_; return *this; }
    reference operator*() { 
        // Does not validate the argument!
        mer_.from_chars(txtPtr_ + (*saPtr_)[curr_->begin()]);
        intRep_ = mer_.get_bits(0, 2*k_);
        return intRep_; 
    }
    /*
    pointer operator->() { 
        return &(curr_->first); 
    }
    */
    bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }
    bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }
    bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }
    bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }
    
private:
    rapmap::utils::my_mer mer_;
    uint64_t intRep_;
    Iter curr_;
    unsigned int k_;
    std::vector<IndexT>* saPtr_{nullptr};
    const char* txtPtr_{nullptr};
};

template <typename IterT>
class KVProxy {
public:
    typedef KVProxy<IterT> self_type;
    typedef typename std::iterator_traits<IterT>::value_type ValueT;
    typedef std::pair<uint64_t, ValueT> value_type;
    typedef value_type& reference;
    typedef value_type* pointer;

    KVProxy(uint64_t mer, IterT it, bool isEnd = false) : curr_(it) {
        if(!isEnd) {pair_ = std::make_pair(mer, *it);}
    }
    reference operator*() { return pair_; }
    pointer operator->() { return &pair_; }
    bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }
    bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }
    bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }
    bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }
private:
    IterT curr_;
    std::pair<uint64_t, ValueT> pair_;
};


// Unlike the standard "generic" BooMap, the frugal variant
// *does not* store the key.  Rather, it assumes we have a
// pointer to the suffix array, and it "spot checks" the index 
// returned by the perfect hash by ensuring that the suffix at
// the corresponding offset starts with the query k-mer.
template <typename KeyT, typename ValueT>
class FrugalBooMap {
public:
    using HasherT = boomphf::SingleHashFunctor<KeyT>;
    using BooPHFT = boomphf::mphf<KeyT, HasherT>;
    typedef typename ValueT::index_type IndexT;
    using IteratorT = KVProxy<typename std::vector<ValueT>::iterator>;

    //using IteratorT = typename std::vector<std::pair<KeyT, ValueT>>::iterator;

    FrugalBooMap() : built_(false) {}
    void setSAPtr(std::vector<IndexT>* saPtr) { saPtr_ = saPtr; }
    void setTextPtr(const char* txtPtr, size_t textLen) { txtPtr_ = txtPtr; textLen_ = textLen; }

    void add(KeyT&& k, ValueT&& v) {
        // In the frugal map, we don't even keep the key!
        data_.emplace_back(v);
    }


    bool validate_hash(){
        auto k_ = rapmap::utils::my_mer::k();
        rapmap::utils::my_mer mer_;
        for( auto& e : data_ ) {
            auto mer = getKmerFromInterval_(e);
            if (mer != getKmerFromInterval_(data_[boophf_->lookup(mer)]) ) {
                std::cerr << "lookup of " << mer << " failed!\n";
            }
        }
        return true;
    }


    bool build(int nthreads=1) {
        size_t numElem = data_.size();
        k_ = rapmap::utils::my_mer::k();
        twok_ = 2 * k_; 
        KeyProxyIterator<decltype(data_.begin()), IndexT> kb(data_.begin(), saPtr_, txtPtr_, k_);
        KeyProxyIterator<decltype(data_.begin()), IndexT> ke(data_.end(), saPtr_, txtPtr_, k_);
        auto keyIt = boomphf::range(kb, ke);
        BooPHFT* ph = new BooPHFT(numElem, keyIt, nthreads);
        boophf_.reset(ph);
        std::cerr << "reordering keys and values to coincide with phf ... ";
        /*
        std::vector<size_t> inds; inds.reserve(data_.size());
        for (size_t i = 0; i < data_.size(); ++i) {
            inds.push_back(ph->lookup(data_[i].first));
        }
        reorder_destructive_(inds.begin(), inds.end(), data_.begin());
        */
        reorder_fn_();
        validate_hash();
        std::cerr << "done\n";
        built_ = true;
        return built_;
    }

    inline IteratorT find(const KeyT& k) {
        auto intervalIndex = boophf_->lookup(k);
        auto ind = data_[intervalIndex].begin();

        if (ind < saPtr_->size()) {
            auto textInd = (*saPtr_)[ind];
            if (textInd + k < textLen_) {
                rapmap::utils::my_mer m(txtPtr_ + textInd);
                // If what we find matches the key, return the iterator
                // otherwise we don't have the key (it must have been here if it
                // existed).
                return (m == k) ? IteratorT(m, data_.begin() + intervalIndex) : end();
            } 
        } 
        // The search sould have been invalid, the key must not exist.
        return end();
    }
    
    /**
     * NOTE: This function *assumes* that the key is in the hash.
     * If it isn't, you'll get back a random element!
     */
    /*
    inline ValueT& operator[](const KeyT& k) {
        auto ind = boophf_->lookup(k);
        return (ind < data_.size() ? data_[ind].second : data_[0].second);
    }
    */
    
    inline IteratorT begin() { return IteratorT(0, data_.begin()); }
    inline IteratorT end() { return IteratorT(0, data_.end(), true); }
    inline IteratorT cend() const { return IteratorT(0, data_.cend(), true); }
    inline IteratorT cbegin() const { return IteratorT(0, data_.cbegin()); }
    
    void save(const std::string& ofileBase) {
        if (built_) {
            std::string hashFN = ofileBase + ".bph";
            // save the perfect hash function
            {
                std::ofstream os(hashFN, std::ios::binary);
                if (!os.is_open()) {
                    std::cerr << "BooM: unable to open output file [" << hashFN << "]; exiting!\n";
                    std::exit(1);
                }
                boophf_->save(os);
                os.close();
            }
            // and the values
            std::string dataFN = ofileBase + ".val";
            {
                std::ofstream valStream(dataFN, std::ios::binary);
                if (!valStream.is_open()) {
                    std::cerr << "BooM: unable to open output file [" << dataFN << "]; exiting!\n";
                    std::exit(1);
                }
                {
                    cereal::BinaryOutputArchive outArchive(valStream);
                    outArchive(data_);
                }
                valStream.close();
            }
        }
    }
    
    void load(const std::string& ofileBase) {
        std::string hashFN = ofileBase + ".bph";
        std::string dataFN = ofileBase + ".val";

        if ( !FileExists_(hashFN.c_str()) ) {
            std::cerr << "BooM: Looking for perfect hash function file [" << hashFN << "], which doesn't exist! exiting.\n";
            std::exit(1);
        }
        if ( !FileExists_(dataFN.c_str()) ) {
            std::cerr << "BooM: Looking for key-value file [" << dataFN << "], which doesn't exist! exiting.\n";
            std::exit(1);
        }

        // load the perfect hash function
        {
            boophf_.reset(new BooPHFT);
            std::ifstream is(hashFN, std::ios::binary);
            boophf_->load(is);
            is.close();
        }
        // and the values
        {
            std::ifstream dataStream(dataFN, std::ios::binary);
            {
                cereal::BinaryInputArchive inArchive(dataStream);
                inArchive(data_);
            }
            dataStream.close();
        }

        k_ = rapmap::utils::my_mer::k();
        twok_ = 2*k_;
        built_ = true;
    }

private:
    // Taken from http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
    bool FileExists_(const char *path) {
        struct stat fileStat;
        if ( stat(path, &fileStat) ) {
            return false;
        }
        if ( !S_ISREG(fileStat.st_mode) ) {
            return false;
        }
        return true;
    }

    inline KeyT getKmerFromInterval_(ValueT& ival) {
        auto m = mer_; // copy the global mer to get k-mer object
        mer_.from_chars(txtPtr_ + (*saPtr_)[ival.begin()]);
        return mer_.get_bits(0, twok_);
    }

    void reorder_fn_()  {
        /* Adapted from code at: http://blog.merovius.de/2014/08/12/applying-permutation-in-constant.html */
        // Note, we can actually do this with out the bitvector by using the high-order bit 
        // of the start of the suffix array intervals (since they are signed integers and negative
        // positions are forbidden). 
        std::vector<bool> bits(data_.size(), false);
        for ( size_t i = 0; i < data_.size(); ++i ) {
            if (!bits[i]) {
                decltype(data_.front()) v = data_[i];
                auto j = boophf_->lookup(getKmerFromInterval_(data_[i]));
                while (i != j) {
                    auto pj = boophf_->lookup(getKmerFromInterval_(data_[j]));
                    std::swap(data_[j], v);
                    bits[j] = 1;
                    j = pj; 
                }
                data_[i] = v;
            }
        }
    }

    rapmap::utils::my_mer mer_;
    std::vector<IndexT>* saPtr_;
    const char* txtPtr_; 
    size_t textLen_;
    bool built_;
    std::vector<ValueT> data_;
    std::unique_ptr<BooPHFT> boophf_{nullptr};
    unsigned int k_;
    unsigned int twok_;
};
#endif // __BOO_MAP_FRUGAL__
