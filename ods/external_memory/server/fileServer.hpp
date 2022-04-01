#pragma once
#include <cinttypes>

#include "common/lrucache.hpp"

#define DEFAULT_FILENAME "./storage.bin"
// This has the same interface as an array that transfers data from a file
//
namespace EM {
namespace FileServer {
// T is typically a LargeBucket
//
template<typename T, uint64_t CACHE_SIZE=ORAM_SERVER__CACHE_SIZE>
struct FileServer {
  static inline constexpr auto sizeOfT = sizeof(T);
  typedef uint64_t IndexType;
  IndexType N;
  _ORAM::Cache<T,CACHE_SIZE> cache;
  std::fstream lbios;

  FileServer(IndexType _N) : N(_N)  {
    PROFILE_F();
    
    lbios.open(DEFAULT_FILENAME, std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::trunc);
    Assert(lbios.is_open());
    for (IndexType i=0; i<N; i++) {
      uint8_t zeros[sizeOfT] = {0};
      memset(zeros, 0, sizeOfT);
      lbios.write((char*)zeros, sizeOfT);
    }
    Assert(lbios.is_open());
    lbios.flush();
    X_LOG("ORAMClientInterface: Done allocating file (", lbios.tellp(), " bytes written)");
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
    std::streampos filePos = i * sizeOfT;
    lbios.seekp(filePos);
    lbios.write((char*)&in, sizeOfT);
    // lbios.flush();
  }

  void Read(const IndexType i, typename T::Encrypted_t & out) {
    Assert(i < N);
    std::streampos filePos = i * sizeOfT;
    lbios.seekg(filePos);
    lbios.read((char*)&out, sizeOfT);
    // lbios.flush();
  }
}; // struct FileServer
} // namespace FileServer
} // namespace _ORAM