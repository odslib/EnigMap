#pragma once
#include "oram/common/bucket.hpp"

namespace _ORAM::RingORAM::Bucket {
  template <
      unsigned int Z=ORAM__Z 
    , unsigned int S=ORAM__S>
  struct PublicMetadata {
    static inline constexpr unsigned int BUCKET_SIZE = S + Z;
    static inline constexpr unsigned int _S = S;
    static inline constexpr unsigned int _Z = Z;

    uint64_t counter;
    bool valid[BUCKET_SIZE];
    bool invalidated;
    // UNDONE(): overload == correctly for obliviousness
    bool operator==(PublicMetadata const&) const = default;
  };
  static_assert(IS_POD<PublicMetadata<>>());


  template <
      unsigned int Z=ORAM__Z 
    , unsigned int S=ORAM__S>
  struct PrivateMetadata {
    static inline constexpr unsigned int BUCKET_SIZE = S + Z;
    static inline constexpr unsigned int _S = S;
    static inline constexpr unsigned int _Z = Z;
    
    // UNDONE(): this could be just Z, investigate what is better
    //
    _ORAM::ORAMAddress addresses[BUCKET_SIZE];

    // UNDONE(): overload correctly == for obliviousness
    bool operator==(PrivateMetadata const&) const = default;

    using Encrypted_t = Encrypted<PrivateMetadata>;
  };
  static_assert(IS_POD<PrivateMetadata<>>());

  template <unsigned int Z=ORAM__Z, unsigned int S=ORAM__S,bool ENCRYPT_BLOCKS=ORAM__ENCRYPT_BLOCKS>
  struct BucketMetadata : MixedEncryptable<PublicMetadata<Z,S>, PrivateMetadata<Z,S> > 
  {
    using Base = MixedEncryptable<PublicMetadata<Z,S>, PrivateMetadata<Z,S>>;
    using PublicData_t = typename Base::PublicData_t;
    using PrivateData_t = typename Base::PrivateData_t; 
    static inline constexpr unsigned int BUCKET_SIZE = S + Z;
    static inline constexpr unsigned int _S = S;
    static inline constexpr unsigned int _Z = Z;

    friend std::ostream &operator<<(std::ostream &o, const BucketMetadata &x)
    {
      if (x.pub.invalidated) {
        o << "{INV}";
        return o;
      }
      o << "{" << x.pub.counter << "," << std::endl;
      for (int i=0; i<BUCKET_SIZE; i++) {
        if (x.priv.addresses[i] == ORAMAddress::DUMMY()) {
          continue;
        }
        if (x.pub.valid[i]) {
          o << "   ";
        } else {
          o << " X ";
        }
        o << i << ": " << x.priv.addresses[i] << std::endl;
      }
      o << "}";
      return o;
    }

    
    static INLINE CLANG_OR_GCC(constexpr, consteval) BucketMetadata DUMMY() {
      BucketMetadata ret;
      ret.pub.counter = 0;
      ret.pub.invalidated = false;
      for (int i=0; i<BUCKET_SIZE;i++) {
        ret.pub.valid[i] = true;
        ret.priv.addresses[i] = ORAMAddress::DUMMY();
      }
      return ret;
    }

    using Encrypted_t = typename std::conditional_t<ENCRYPT_BLOCKS, MixedEncrypted<BucketMetadata>, NonEncrypted<BucketMetadata> >;
  };

  template<
    typename Block=_ORAM::Block::Block<> 
  , unsigned int Z=ORAM__Z 
  , unsigned int S=ORAM__S
  , bool ENCRYPT_BLOCKS=ORAM__ENCRYPT_BLOCKS
  , typename BucketMetadata=BucketMetadata<Z,S,ENCRYPT_BLOCKS> >
  using Bucket = _ORAM::Bucket::Bucket<Block,Z,S,ENCRYPT_BLOCKS,BucketMetadata>;

} // namespace _ORAM::RingORAM::Bucket