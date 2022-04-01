#pragma once
#include <cinttypes>
#include <vector>

#include "common/lrucache.hpp"

namespace EM {
namespace MemServer {
// T is typically a LargeBucket
//
template<typename T, uint64_t CACHE_SIZE=ORAM_SERVER__CACHE_SIZE>
struct MemServer {
  static inline constexpr auto sizeOfT = sizeof(T);
  typedef uint64_t IndexType;
  IndexType N;
  std::vector<typename T::Encrypted_t> data;
  _ORAM::Cache<T,100000000> cache;

  MemServer(IndexType _N) : N(_N), data(_N) {
    PROFILE_F();

    uint64_t byteCount = ((uint64_t)sizeof(typename T::Encrypted_t))*_N;
    // std::cerr << "Lbsize: " << sizeof(typename T::Encrypted_t) << std::endl;
    // std::cerr << "Datasize: " << _N << std::endl;
    // std::cerr << "ORAMClientInterface: Creating in memory file (" << byteCount << " bytes)" << std::endl;

    memset(&data[0], 0, byteCount);
    X_LOG("ORAMClientInterface: Done creating in memory file (", byteCount , " bytes)");
  }

  T& Access(const IndexType i) {
    Assert(i < N, i, N);
    if (!cache.CheckContains(i)) {
      if (cache.IsFull()) {
        IndexType evictedIndex;
        {
          T& evicted = cache.GetNextToEvict(evictedIndex);
          typename T::Encrypted_t evictedEnc;
          evictedEnc.Encrypt(evicted);
          Write(evictedIndex, evictedEnc);
        }
        cache.EvictLRU(evictedIndex);
      }
      typename T::Encrypted_t retEnc;
      T ret;
      Read(i, retEnc);
      retEnc.Decrypt(ret);
      std::ignore = cache.Insert(i, ret);
    }

    return cache.Access(i);
  }

  void Write(const IndexType i, const typename T::Encrypted_t & in) {
    Assert(i < N);
    data[i] = in;
  }

  void Read(const IndexType i, typename T::Encrypted_t & out) {
    Assert(i < N);
    out = data[i];
  }
}; // struct MemServer
} // namespace MemServer
} // namespace _ORAM