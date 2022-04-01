#pragma once
#include <inttypes.h>
#include <utility>

#include "common/utils.hpp"
#include "common/encutils.hpp"
#include "common/tracing/tracer.hpp"

namespace Concepts {
  template<typename T>
  concept Encryptable = requires(typename T::Encrypted_t et) {
    { et };
  };
}

// "random" (:P) AES key:
// UNDONE(): move key to oram_client_interface.hpp
inline constexpr uint8_t KEY[16] = {0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41};

template<typename T>
// requires (IS_POD<T>())
struct Encrypted {
  static_assert(IS_POD<T>());


  static constexpr uint64_t SIZE = sizeof(T);
  uint8_t data[SIZE]; //We don't need to adjust this because CTR modes don't need padding.
  uint8_t iv[AES_BLOCK_SIZE];
  uint8_t tag[AES_BLOCK_SIZE];

  // Encrypted& operator=(const Encrypted&) = delete;
  // Encrypted(const Encrypted&) = delete;

  inline void Encrypt(const T& in) {
    // PROFILE_F();
    GetRand16(iv);
    aes_256_gcm_encrypt(SIZE, const_cast<uint8_t*>(reinterpret_cast<const uint8_t*>(&in)), KEY, iv, tag, data);
  }

  inline void Decrypt(T& out) /*const*/ {
    // PROFILE_F();
    bool r = aes_256_gcm_decrypt(SIZE, data, KEY, iv, tag, reinterpret_cast<uint8_t*>(&out));
    Assert(r);
    IGNORE_UNUSED(r);
  }

  
  friend std::ostream &operator<<(std::ostream &o, const Encrypted &x)
  {
    T v;
    const_cast<Encrypted&>(x).Decrypt(v);
    o << "E{" << v << "}";
    return o;
  }

  bool operator==(Encrypted& o) {
    T d1;
    Decrypt(d1);
    T d2;
    o.Decrypt(d2);
    return d1 == d2;
  }
};
static_assert(IS_POD<Encrypted<int>>());


// Structure with same interfaces as Encrypted, but with
// everything as plaintext for debugging.
//
template<typename T>
struct NonEncrypted {
  static_assert(IS_POD<T>());

  static constexpr uint64_t SIZE = sizeof(T);
  T data;

  inline void Encrypt(const T& in) {
    data = in;
  }

  inline void Decrypt(T& out) /*const*/ {
    out = data;
  }

  friend std::ostream &operator<<(std::ostream &o, const NonEncrypted &x)
  {
    T v;
    const_cast<NonEncrypted&>(x).Decrypt(v);
    o << "E{" << v << "}";
    return o;
  }

  bool operator==(const NonEncrypted& o) const {
    return data == o.data;
  }
};
static_assert(IS_POD<NonEncrypted<int>>());


template<typename PublicData, typename PrivateData>
requires (IS_POD<PublicData>())
      && (IS_POD<PrivateData>())
      && ::Concepts::Encryptable<PrivateData>
struct MixedEncryptable {
  using PublicData_t = PublicData;
  using PrivateData_t = PrivateData;
  PublicData pub;
  PrivateData priv;

  friend std::ostream &operator<<(std::ostream &o, const MixedEncryptable& x)
  {
    o << "ME{"
      << x.pub
      << ", "
      << x.priv
      << "}";
    return o;
  }

  bool operator==(const MixedEncryptable& o) const {
    return (pub == o.pub) * (priv == o.priv);
  }

  // Classes that extend this need to declare Encrypted_t.
  //
};

template <typename PublicData, typename PrivateData, typename T>
concept ExtendsMixedEncrytable = std::is_base_of<MixedEncryptable<PublicData, PrivateData>, T>::value;

template<typename T>
struct MixedEncrypted
{
  using PublicData_t = typename T::PublicData_t;
  using PrivateData_t = typename T::PrivateData_t;
  typename PrivateData_t::Encrypted_t _priv;
  PublicData_t _pub;

  static_assert(std::is_base_of<MixedEncryptable<PublicData_t, PrivateData_t>, T>::value);

  inline void Encrypt(const T& in) {
    _pub = in.pub;
    _priv.Encrypt(in.priv);
  }

  inline void Decrypt(T& out) /*const*/ {
    out.pub = _pub;
    _priv.Decrypt(out.priv);
  }

  friend std::ostream &operator<<(std::ostream &o, const MixedEncrypted& x)
  {
    o << "ME{"
      << x.pub
      << ", "
      << x.priv
      << "}";
    return o;
  }

  bool operator==(const MixedEncrypted& o) const {
    return (_pub == o._pub) * (_priv == o._priv);
  }
};