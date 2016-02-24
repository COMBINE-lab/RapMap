#ifndef __MPH_MAP_HPP__
#define __MPH_MAP_HPP__

#include <fstream>

#include "cereal/types/vector.hpp"
#include "cereal/types/utility.hpp"
#include "cereal/archives/binary.hpp"

#include "RapMapFileSystem.hpp"
#include "RapMapUtils.hpp"

#include "emphf/common.hpp"
#include "emphf/mphf.hpp"
#include "emphf/base_hash.hpp"
#include "emphf/perfutils.hpp"
#include "emphf/mmap_memory_model.hpp"
#include "emphf/hypergraph_sorter_scan.hpp"

template <typename K, typename V>
class MPHMap {
public:
    using iterator_t = typename std::vector<V>::iterator;
    MPHMap() = default;

    
  inline void load(const std::string& indexDir) { 
    std::string hashFname = indexDir + "hash.bin";
    std::string intervalFname = indexDir + "kintervals.bin";
    if ( !rapmap::fs::FileExists(hashFname.c_str()) ) {
      std::cerr << "Looking for index-related file " << hashFname << ", which doesn't exist! exiting.\n";
      std::exit(1);
    }
    if ( !rapmap::fs::FileExists(intervalFname.c_str()) ) {
      std::cerr << "Looking for index-related file " << hashFname << ", which doesn't exist! exiting.\n";
      std::exit(1);
    }
    // Load the hash
    {
      std::ifstream hashStream(hashFname, std::ios::binary);
      hash_.load(hashStream);
      hashStream.close();
    }
	
    // Load the k-mer intervals
    std::ifstream intervalStream(intervalFname, std::ios::binary);
    {
      cereal::BinaryInputArchive inArchive(intervalStream);
      inArchive(varray_);
    }
    intervalStream.close();
    numVals_ = varray_.size();
  }

      /**
       * Given a k-mer key (key), return an iterator 
       * pointing to the interval of this k-mer's occurences
       * in the suffix array.  If the k-mer isn't present in 
       * the suffix array, returns end().
       */
    inline iterator_t find(const K& key) {
      // Get the index from the perfect hash
      auto idx = hash_.lookup(key, adaptor_);
      // The index is valid if it's within the range, and
      // it maps to an interval that matches the requested key.
      // Otherwise, it's an invalid index and we return varray_.end(); 
      return (idx < numVals_) ? 
	     ((varray_[idx].first == key) ? (begin() + idx) : varray_.end()) : varray_.end();
    }

    inline iterator_t cend() const { return varray_.cend(); }
    inline iterator_t cbegin() const { return varray_.cbegin(); }
    inline iterator_t end() { return varray_.end(); }
    inline iterator_t begin() { return varray_.begin(); }

  private:

  struct MerAdaptor {
    emphf::byte_range_t operator()(const uint64_t& ival) const {
      const uint8_t* buf = reinterpret_cast<uint8_t const*>(&ival);
      const uint8_t* end = buf + sizeof(K);
      return emphf::byte_range_t(buf, end);
    }
  };


    using mphf_t = emphf::mphf<emphf::jenkins64_hasher>;
    mphf_t hash_;
    MerAdaptor adaptor_;
    std::vector<V> varray_;
    size_t numVals_;
};

#endif // __MPH_MAP_HPP__
