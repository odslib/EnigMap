#ifndef NDEBUG
#ifndef ENCLAVE_MODE
#define OCI_ASSERTIONS
#endif
#endif

#ifdef OCI_ASSERTIONS
#include <boost/functional/hash.hpp>
#endif

namespace _ORAM::ORAMClientInterface {
#ifdef OCI_ASSERTIONS
  template<typename Block, typename Bucket>
  struct ORAMInterfaceVerifier {
    using Block_t = Block;
    using BucketMetadata_t = typename Bucket::BucketMetadata_t;
    typedef std::tuple<Position, Index, Index> Triplet;
    typedef std::tuple<Position, Index> Duplet;
  
    template<typename K>
    struct KeyHash {
        std::size_t operator()(const K &key) const
        {
            return boost::hash_value(key);
        }
    };
    uint64_t L_;
    std::unordered_map<Triplet, Block_t, KeyHash<Triplet> > checkdata_block;
    std::unordered_map<Duplet, BucketMetadata_t, KeyHash<Duplet> > checkdata_metadata;

    ORAMInterfaceVerifier(uint64_t L) : L_(L) { }

    void AssertSameMetadata(const Index& pos, const Index& depth, const BucketMetadata_t& ret) {
      Position realPos = Indexers::GetArrIndex(L_, pos, depth);
    
      Duplet location = Duplet{realPos, depth};
      if (checkdata_metadata.contains(location)) {
        BucketMetadata_t& expected = checkdata_metadata[location]; 
        Assert(expected == ret, "\n", expected, "\n", ret, pos, realPos, depth, ret, ret);
      } else {
        Assert(ret == BucketMetadata_t::DUMMY(), "\n", ret, pos, realPos, depth, ret, ret);
      }
    }

    void AssertSameBlock(const Index& pos, const Index& depth, const Index& offset, Block_t& ret) {
       Position realPos = Indexers::GetArrIndex(L_, pos, depth);

      // if (ret.data.data.k != uint64_t(-1)) {
      //   TRACE_FUNCTION(pos, realPos, depth, offset, ret);
      // }

      Triplet location = Triplet{realPos, depth, offset};
      if (checkdata_block.contains(location)) {
        auto& expected = checkdata_block[location];
        Assert(expected == ret, "\n", expected, "\n", ret);
      } else {
        Assert(ret == Block_t::DUMMY(), "\n", ret);
      }
    }

    void UpdateMetadata(const Index& pos, const Index& depth, const BucketMetadata_t& val) {
      Position realPos = Indexers::GetArrIndex(L_, pos, depth);
      // TRACE_FUNCTION(pos, realPos, depth, val, &bucket);

      Duplet location = Duplet{realPos, depth};
      checkdata_metadata[location] = val;
    }

    void UpdateBlock(const Index& pos, const Index& depth, const Index& offset, const Block_t& val) {
      Assert(offset < Bucket::BUCKET_SIZE);
      Position realPos = Indexers::GetArrIndex(L_, pos, depth);

      Triplet location = Triplet{realPos, depth, offset};
      checkdata_block[location] = val;
    }
  };

  #define OCI_ONLY(...) __VA_ARGS__
#else
  #define OCI_ONLY(...) 
#endif
};