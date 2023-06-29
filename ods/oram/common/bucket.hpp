#pragma once
#include <cinttypes>


#pragma once
#include "common/defs.hpp"
#include "common/encrypted.hpp"
#include "oram/common/block.hpp"
#include "oram/common/concepts.hpp"

namespace _ORAM::Bucket {
// These DefaultXMetadata are just examples.
// Each oram overrides them with whatever metadata a bucket should hold.
//
template <bool make_template=true>
struct DefaultPublicMetadata {
  bool operator==(DefaultPublicMetadata const&) const = default;
};
static_assert(IS_POD<DefaultPublicMetadata<>>());


template <unsigned int Z=ORAM__Z>
struct DefaultPrivateMetadata {
  _ORAM::ORAMAddress addresses[Z];

  // UNDONE(): overload correctly == for obliviousness
  bool operator==(DefaultPrivateMetadata const&) const = default;
  using Encrypted_t = Encrypted<DefaultPrivateMetadata>;
};
static_assert(IS_POD<DefaultPrivateMetadata<>>());


template <unsigned int Z=ORAM__Z, unsigned int S=ORAM__S,bool ENCRYPT_BLOCKS=ORAM__ENCRYPT_BLOCKS>
requires (S==0)
struct DefaultBucketMetadata : MixedEncryptable<DefaultPublicMetadata<>, DefaultPrivateMetadata<Z> > 
{
  using Base = MixedEncryptable<DefaultPublicMetadata<>, DefaultPrivateMetadata<Z>>;
  using PublicData_t = typename Base::PublicData_t;
  using PrivateData_t = typename Base::PrivateData_t; 
  static inline constexpr unsigned int _Z = Z;
  static inline constexpr unsigned int BUCKET_SIZE = _Z;

  #ifndef ENCLAVE_MODE
  friend std::ostream &operator<<(std::ostream &o, const DefaultBucketMetadata &x)
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
  #endif

  static INLINE CLANG_OR_GCC(constexpr, consteval) DefaultBucketMetadata DUMMY() {
    DefaultBucketMetadata ret;
    for (int i=0; i<Z-1; i++) {
      ret.priv.addresses[i] = ORAMAddress::DUMMY();
    }
    return ret;
  }

  using Encrypted_t = typename std::conditional_t<ENCRYPT_BLOCKS, MixedEncrypted<DefaultBucketMetadata>, NonEncrypted<DefaultBucketMetadata> >;
};

template<
    typename Block=_ORAM::Block::Block<> 
  , unsigned int Z=ORAM__Z 
  , unsigned int S=0
  , bool ENCRYPT_BLOCKS=ORAM__ENCRYPT_BLOCKS
  , typename BucketMetadata=DefaultBucketMetadata<Z,S,ENCRYPT_BLOCKS>
  , bool ENCRYPT_BUCKETS=ORAM__ENCRYPT_BUCKETS>
struct Bucket {
  static inline constexpr unsigned int BUCKET_SIZE = S + Z;
  static inline constexpr unsigned int _S = S;
  static inline constexpr unsigned int _Z = Z;
  using Block_t = Block;
  using BucketMetadata_t = BucketMetadata;
  // </struct BucketMetadata>
  
  static_assert(IS_POD<Block>());
  static_assert(IS_POD<BucketMetadata>());

  // UNDONE(): We might also want to include in bucket information about its
  // children signatures (this way we have a merkle tree for free).
  //
  BucketMetadata_t md;
  Block_t blocks[BUCKET_SIZE];

  #ifndef ENCLAVE_MODE
  friend std::ostream &operator<<(std::ostream &o, const Bucket &x)
  {
    if (x.md.pub.invalidated) {
      o << "{INV}";
      return o;
    }
    o << "{" << x.md.pub.counter << "," << std::endl;
    for (int i=0; i<BUCKET_SIZE; i++) {
      if (x.md.priv.addresses[i] == ORAMAddress::DUMMY()) {
        continue;
      }
      if (x.md.pub.valid[i]) {
        o << "   ";
      } else {
        o << " X ";
      }
      o << i << ": " << x.md.priv.addresses[i];
      o << x.blocks[i];

      o << std::endl;
    }
    o << "}";
    return o;
  }
  #endif

  // UNDONE(): We want the bucket to forward encryptions, as, for instance, 
  // Metadata is only partially encrypted.
  //
  using Encrypted_t = std::conditional_t<ENCRYPT_BUCKETS, Encrypted<Bucket>, NonEncrypted<Bucket> >; 
};
// </struct Bucket>

static_assert(IS_POD<Bucket<>>());
static_assert(IS_POD<Bucket<>::Encrypted_t>());

} // namespace _ORAM::Bucket




namespace _ORAM::LargeBucket {
  template<typename Bucket=_ORAM::Bucket::Bucket<>
    , bool ENCRYPT_LARGE_BUCKETS=ORAM_SERVER__ENCRYPT_LARGE_BUCKETS
    , unsigned int LEVELS_PER_PACK=ORAM_SERVER__LEVELS_PER_PACK>
  struct LargeBucket {
    static inline constexpr uint64_t BUCKETS_PER_PACK = (1<<LEVELS_PER_PACK)-1; // 1 + 2 + 4 + 8

    static auto assertions() {
      static_assert( sizeof(Bucket) * (BUCKETS_PER_PACK) == sizeof(LargeBucket));
      static_assert(BUCKETS_PER_PACK == -1 + (1<<(LEVELS_PER_PACK)));
    }

    Bucket buckets[BUCKETS_PER_PACK];
    using Encrypted_t = std::conditional_t<ENCRYPT_LARGE_BUCKETS, Encrypted<LargeBucket>, NonEncrypted<LargeBucket>>;
  }; // struct LargeBucket
  
  
static_assert(IS_POD<LargeBucket<>>());
static_assert(IS_POD<LargeBucket<>::Encrypted_t>());
} // namespace _ORAM::LargeBucket