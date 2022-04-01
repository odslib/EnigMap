#pragma once
#include <cinttypes>

#include "common/lrucache.hpp"

// Similar to fileServer.hpp, but separated to be used within an enclave.
//
namespace EM {
namespace EnclaveFileServer {
// T is typically a LargeBucket
//
template<typename T, uint64_t CACHE_SIZE=ORAM_SERVER__CACHE_SIZE>
struct EnclaveFileServer {
  static inline constexpr auto sizeOfT = sizeof(T);
  typedef uint64_t IndexType;
  IndexType N;
  _ORAM::Cache<T,CACHE_SIZE> cache;

  EnclaveFileServer(IndexType _N) : N(_N)  {
    TRACE_FUNCTION(_N);
    ocall_InitServer(sizeOfT, _N);
  }

  T& Access(const IndexType i) {
    Assert(i < N, i, N);
    if (!cache.CheckContains(i)) {
      bool wasFull = false;
      IndexType evictedIndex;
      typename T::Encrypted_t evictedEnc;

      if (cache.IsFull()) {
        {
          T& evicted = cache.GetNextToEvict(evictedIndex);          
          evictedEnc.Encrypt(evicted);
        }
        cache.EvictLRU(evictedIndex);
      }
      typename T::Encrypted_t retEnc;
      T ret;

      if (wasFull) {
        Swap(evictedIndex, i, evictedEnc, retEnc);
      } else {
        Read(i, retEnc);
      }
      retEnc.Decrypt(ret);
      std::ignore = cache.Insert(i, ret);
    }

    return cache.Access(i);
  }

  void Swap(const IndexType index_in, const IndexType index_out, const typename T::Encrypted_t & in, typename T::Encrypted_t & out) {
    Assert(index_in < N);
    Assert(index_out < N);
    
    uint8_t aux[4096];
    memcpy(aux, &index_in, sizeOfT);
    ocall_WritePage(sizeOfT, index_in, aux);
    ocall_ReadPage(sizeOfT, index_out, aux);
    memcpy(&out, aux, sizeOfT);
  }

  void Write(const IndexType i, const typename T::Encrypted_t & in) {
    Assert(i < N);
    uint8_t aux[4096];
    memcpy(aux, &in, sizeOfT);
    ocall_WritePage(sizeOfT, i, aux);
  }

  void Read(const IndexType i, typename T::Encrypted_t & out) {
    Assert(i < N);
    uint8_t aux[4096];
    ocall_ReadPage(sizeOfT, i, aux);
    memcpy(&out, aux, sizeOfT);
  }
}; // struct EnclaveFileServer
} // namespace EnclaveFileServer
} // namespace _ORAM