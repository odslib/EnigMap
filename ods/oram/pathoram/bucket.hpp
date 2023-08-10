#pragma once
#include "oram/common/bucket.hpp"

namespace _ORAM::PathORAM::Bucket {
  template <unsigned int Z=ORAM__Z>
  struct BucketMetadata
  {
    static inline constexpr unsigned int BUCKET_SIZE = Z;
    static inline constexpr unsigned int _Z = Z;
    
    _ORAM::ORAMAddress addresses[BUCKET_SIZE];

    #ifndef ENCLAVE_MODE
    friend std::ostream &operator<<(std::ostream &o, const BucketMetadata &x)
    {
      o << "{" << std::endl;
      for (int i=0; i<BUCKET_SIZE; i++) {
        if (x.addresses[i] == ORAMAddress::DUMMY()) {
          continue;
        }
        o << i << ": " << x.addresses[i] << std::endl;
      }
      o << "}";
      return o;
    }
    #endif

    
    static INLINE CLANG_OR_GCC(constexpr, consteval) BucketMetadata DUMMY() {
      BucketMetadata ret;
      for (int i=0; i<BUCKET_SIZE;i++) {
        ret.addresses[i] = ORAMAddress::DUMMY();
      }
      return ret;
    }

    bool operator==(BucketMetadata const&) const = default;

    using Encrypted_t = typename std::conditional_t<false, NonEncrypted<BucketMetadata>, NonEncrypted<BucketMetadata> >;
  };

  template<
    typename Block=_ORAM::Block::Block<> 
  , unsigned int Z=ORAM__Z 
  , typename BucketMetadata=BucketMetadata<Z> >
  using Bucket = _ORAM::Bucket::Bucket<Block,Z,0,false,BucketMetadata,false>;
  
} // namespace _ORAM::PathORAM::Bucket