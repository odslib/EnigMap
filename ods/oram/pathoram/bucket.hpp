#pragma once
#include "oram/common/bucket.hpp"

namespace _ORAM::PathORAM::Bucket {
  struct PublicMetadata {
    // UNDONE(): overload == correctly for obliviousness
    bool operator==(PublicMetadata const&) const = default;
  };
  static_assert(IS_POD<PublicMetadata>());


  template <unsigned int Z=ORAM__Z>
  struct PrivateMetadata {
    static inline constexpr unsigned int BUCKET_SIZE = Z;
    static inline constexpr unsigned int _Z = Z;
    
    _ORAM::ORAMAddress addresses[BUCKET_SIZE];

    // UNDONE(): overload correctly == for obliviousness
    bool operator==(PrivateMetadata const&) const = default;

    using Encrypted_t = Encrypted<PrivateMetadata>;
  };
  static_assert(IS_POD<PrivateMetadata<>>());

  template <unsigned int Z=ORAM__Z, bool ENCRYPT_BLOCKS=ORAM__ENCRYPT_BLOCKS>
  struct BucketMetadata : MixedEncryptable<PublicMetadata, PrivateMetadata<Z> > 
  {
    using Base = MixedEncryptable<PublicMetadata, PrivateMetadata<Z>>;
    using PublicData_t = typename Base::PublicData_t;
    using PrivateData_t = typename Base::PrivateData_t; 
    static inline constexpr unsigned int BUCKET_SIZE = Z;
    static inline constexpr unsigned int _Z = Z;

    friend std::ostream &operator<<(std::ostream &o, const BucketMetadata &x)
    {
      o << "{" << std::endl;
      for (int i=0; i<BUCKET_SIZE; i++) {
        if (x.priv.addresses[i] == ORAMAddress::DUMMY()) {
          continue;
        }
        o << i << ": " << x.priv.addresses[i] << std::endl;
      }
      o << "}";
      return o;
    }

    
    static inline CLANG_OR_GCC(constexpr, consteval) BucketMetadata DUMMY() {
      BucketMetadata ret;
      for (int i=0; i<BUCKET_SIZE;i++) {
        ret.priv.addresses[i] = ORAMAddress::DUMMY();
      }
      return ret;
    }

    using Encrypted_t = typename std::conditional_t<ENCRYPT_BLOCKS, MixedEncrypted<BucketMetadata>, NonEncrypted<BucketMetadata> >;
  };

  template<
    typename Block=_ORAM::Block::Block<> 
  , unsigned int Z=ORAM__Z 
  , bool ENCRYPT_BUCKETS=ORAM__ENCRYPT_BLOCKS
  , typename BucketMetadata=BucketMetadata<Z,false> >
  using Bucket = _ORAM::Bucket::Bucket<Block,Z,0,false,BucketMetadata,ENCRYPT_BUCKETS>;

} // namespace _ORAM::PathORAM::Bucket