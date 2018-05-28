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
template <typename Iter, typename IndexT, typename HashMapT>
class KeyProxyIterator {
public:
    typedef KeyProxyIterator<Iter, IndexT, HashMapT> self_type;
    typedef uint64_t value_type;//std::iterator_traits<Iter>::value_type::first_type value_type;
    typedef value_type& reference;
    typedef value_type* pointer;
    typedef std::forward_iterator_tag iterator_category;
    typedef int64_t difference_type;

    KeyProxyIterator(Iter first, HashMapT* hm) : 
        curr_(first), hm_(hm) {}
    KeyProxyIterator operator++() { KeyProxyIterator i = *this; curr_++; return i; }
    KeyProxyIterator operator++(int) { ++curr_; return *this; }
    reference operator*() { 
        intRep_ = hm_->getKmerFromPos_(*curr_, mer_);
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
    HashMapT* hm_{nullptr}; 
};


template <typename IterT, typename ProxyValueT>
class KVProxy {
public:
    typedef KVProxy<IterT, ProxyValueT> self_type;
    typedef typename std::iterator_traits<IterT>::value_type ValueT;
    typedef std::pair<uint64_t, ValueT> value_type;
    typedef std::pair<uint64_t, ProxyValueT>& reference;
    typedef std::pair<uint64_t, ProxyValueT>* pointer;

    KVProxy(uint64_t mer, IterT it, ValueT len, bool isEnd = false) : curr_(it) {
        if(!isEnd) {ProxyValueT x{*it, len}; pair_ = std::make_pair(mer, x);}
    }
    reference operator*() { return pair_; }
    pointer operator->() { return &pair_; }
    bool operator==(const self_type& rhs) { return curr_ == rhs.curr_; }
    bool operator!=(const self_type& rhs) { return curr_ != rhs.curr_; }
    bool operator<(const self_type& rhs) { return curr_ < rhs.curr_; }
    bool operator<=(const self_type& rhs) { return curr_ <= rhs.curr_; }
private:
    IterT curr_;
    std::pair<uint64_t, ProxyValueT> pair_;
};


// Unlike the standard "generic" BooMap, the frugal variant
// *does not* store the key.  Rather, it assumes we have a
// pointer to the suffix array, and it "spot checks" the index 
// returned by the perfect hash by ensuring that the suffix at
// the corresponding offset starts with the query k-mer.
template <typename KeyT, typename ValueT>
class FrugalBooMap {
public:
    using self_type = FrugalBooMap<KeyT, ValueT>;
    using HasherT = boomphf::SingleHashFunctor<KeyT>;
    using BooPHFT = boomphf::mphf<KeyT, HasherT>;
    typedef typename ValueT::index_type IndexT;
    using IteratorT = KVProxy<typename std::vector<IndexT>::iterator, ValueT>;

    //using IteratorT = typename std::vector<std::pair<KeyT, ValueT>>::iterator;

    FrugalBooMap() : built_(false) {}
    void setSAPtr(std::vector<IndexT>* saPtr) { saPtr_ = saPtr; }
    void setTextPtr(const char* txtPtr, size_t textLen) { txtPtr_ = txtPtr; textLen_ = textLen; }

    void add(KeyT&& k, ValueT&& v) {
        // In the frugal map, we don't even keep the key!
        data_.emplace_back(v.begin());
        IndexT l = v.end() - v.begin();
        if (l >= std::numeric_limits<uint8_t>::max()) {
            overflow_[v.begin()] = l;
            lens_.emplace_back(std::numeric_limits<uint8_t>::max());
        } else {
            lens_.emplace_back(static_cast<uint8_t>(l));
        }
    }


    bool validate_hash(){
        for( auto& e : data_ ) {
            rapmap::utils::my_mer kmer(txtPtr_ + (*saPtr_)[e]);
            auto ind = boophf_->lookup(kmer.word(0));
            if (ind >= data_.size()) { 
                rapmap::utils::my_mer km(txtPtr_ + (*saPtr_)[e]);
                std::cerr << "index for " << km << " was " << ind << ", outside bounds of data_ (" << data_.size() << ")\n";
                return false;
            }
            auto mer = getKmerFromInterval_(e);
            if (mer != getKmerFromInterval_(data_[ind]) ) {
                std::cerr << "lookup of " << mer << " failed!\n";
            }
        }
        return true;
    }


    bool build(int nthreads=1) {
        size_t numElem = data_.size();
        KeyProxyIterator<decltype(data_.begin()), IndexT, self_type> kb(data_.begin(), this);
        KeyProxyIterator<decltype(data_.begin()), IndexT, self_type> ke(data_.end(), this);
        auto keyIt = boomphf::range(kb, ke);
        BooPHFT* ph = new BooPHFT(numElem, keyIt, nthreads);
        boophf_.reset(ph);
        std::cerr << "reordering keys and values to coincide with phf ... ";
        reorder_fn_();
        //std::cerr << "validating hash\n";
        //validate_hash();
        std::cerr << "done\n";
        std::cerr << "size of overflow table is " << overflow_.size() << '\n';
        built_ = true;
        return built_;
    }

    inline IteratorT find(const KeyT& k) {
        auto intervalIndex = boophf_->lookup(k);
        if (intervalIndex >= data_.size()) return end();
        auto ind = data_[intervalIndex];
        auto textInd = (*saPtr_)[ind];
        rapmap::utils::my_mer m(txtPtr_ + textInd);

        // If what we find matches the key, return the iterator
        // otherwise we don't have the key (it must have been here if it
        // existed).
        if (m.word(0) == k) {
            IndexT l = *(lens_.begin() + intervalIndex);
            if (l == std::numeric_limits<uint8_t>::max()) {
                l = overflow_[ind];
            }
            return IteratorT(m.word(0), data_.begin() + intervalIndex, ind + l);
        }
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
    
    inline IteratorT begin() { return IteratorT(0, data_.begin(), lens_.front()); }
    inline IteratorT end() { return IteratorT(0, data_.end(), 0, true); }
    inline IteratorT cend() const { return IteratorT(0, data_.cend(), 0, true); }
    inline IteratorT cbegin() const { return IteratorT(0, data_.cbegin(), lens_.front()); }
    
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
                    outArchive(lens_);
                    overflow_.serialize(typename spp_utils::pod_hash_serializer<IndexT, IndexT>(), &valStream);
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
                inArchive(lens_);
                overflow_.unserialize(typename spp_utils::pod_hash_serializer<IndexT, IndexT>(), &dataStream);
            }
            dataStream.close();
        }

        built_ = true;
    }

    inline KeyT getKmerFromInterval_(ValueT& ival) {
        rapmap::utils::my_mer m;// copy the global mer to get k-mer object
        m.fromChars(txtPtr_ + (*saPtr_)[ival.begin()]);
        return m.word(0);
    }

    // variant where we provide an existing mer object
    inline KeyT getKmerFromInterval_(ValueT& ival, rapmap::utils::my_mer& m) {
        m.fromChars(txtPtr_ + (*saPtr_)[ival.begin()]);
        return m.word(0);
    }

    // variant where we provide an existing mer object
    inline KeyT getKmerFromPos_(IndexT pos, rapmap::utils::my_mer& m) {
        m.fromChars(static_cast<decltype(txtPtr_)>(txtPtr_ + (*saPtr_)[pos]));
        return m.word(0);
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

    void reorder_fn_()  {
        /* Adapted from code at: http://blog.merovius.de/2014/08/12/applying-permutation-in-constant.html */
        // Note, we can actually do this with out the bitvector by using the high-order bit 
        // of the start of the suffix array intervals (since they are signed integers and negative
        // positions are forbidden). 
        rapmap::utils::my_mer mer;
        std::vector<bool> bits(data_.size(), false);
        for ( size_t i = 0; i < data_.size(); ++i ) {
            if (!bits[i]) {
                decltype(data_.front()) v = data_[i];
                decltype(lens_.front()) v2 = lens_[i];
                auto j = boophf_->lookup(getKmerFromPos_(data_[i], mer));
                while (i != j) {
                    auto pj = boophf_->lookup(getKmerFromPos_(data_[j], mer));
                    std::swap(data_[j], v);
                    std::swap(lens_[j], v2);
                    bits[j] = 1;
                    j = pj; 
                }
                data_[i] = v;
                lens_[i] = v2;
            }
        }
    }

    std::vector<IndexT>* saPtr_;
    const char* txtPtr_; 
    size_t textLen_;
    rapmap::utils::my_mer mer_;
    bool built_;
    // Starting offset in the suffix array
    std::vector<IndexT> data_;
    // Length of the interval
    std::vector<uint8_t> lens_;
    // Overflow table if interval is >= std::numeric_limits<uint8_t>::max()
    spp::sparse_hash_map<IndexT, IndexT> overflow_;
    std::unique_ptr<BooPHFT> boophf_{nullptr};
};


#endif // __BOO_MAP_FRUGAL__
