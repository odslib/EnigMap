#pragma once
#include <cinttypes>
#include <vector>

#include "common/lrucache.hpp"
#include "common/tracing/tracer.hpp"
#include "external_memory/server/cached.hpp"

namespace EM {
namespace MemServer {
// T is typically a LargeBucket
//
template <typename T, uint64_t MEM_SIZE = (3UL << 25)>
struct MemServerBackend {
  static inline constexpr auto sizeOfT = sizeof(T);
  static inline constexpr auto cacheSize = MEM_SIZE / sizeOfT;
  typedef uint64_t IndexType;
  IndexType N;
  std::vector<typename T::Encrypted_t> data;

  explicit MemServerBackend(IndexType _N) : N(_N), data(_N) {
    // PROFILE_F();

    uint64_t byteCount = ((uint64_t)sizeof(typename T::Encrypted_t)) * _N;
    // std::cerr << "Lbsize: " << sizeof(typename T::Encrypted_t) << std::endl;
    // std::cerr << "Datasize: " << _N << std::endl;
    // std::cerr << "Cachesize: " << sizeof(T) * 512 << " bytes" << std::endl;
    // std::cerr << "ORAMClientInterface: Creating in memory file (" <<
    // byteCount << " bytes)" << std::endl;

    memset(&data[0], 0, byteCount);
    X_LOG("ORAMClientInterface: Done creating in memory file (", byteCount,
          " bytes)");
  }

  MemServerBackend(IndexType _N, const T& defaultVal) : N(_N), data(_N) {
    // PROFILE_F();

    uint64_t byteCount = ((uint64_t)sizeof(typename T::Encrypted_t)) * _N;
    // std::cerr << "Lbsize: " << sizeof(typename T::Encrypted_t) << std::endl;
    // std::cerr << "Datasize: " << _N << std::endl;
    // std::cerr << "Cachesize: " << sizeof(T) * 512 << " bytes" << std::endl;
    // std::cerr << "ORAMClientInterface: Creating in memory file (" <<
    // byteCount << " bytes)" << std::endl;

    for (uint64_t i = 0; i < _N; i++) {
      data[i].Encrypt(defaultVal);
    }
    X_LOG("ORAMClientInterface: Done creating in memory file (", byteCount,
          " bytes)");
  }

  void Write(const IndexType i, const T& in) {
    PERFCTR_INCREMENT(writeCount);
    Assert(i < N);
    data[i].Encrypt(in);
  }

  void Read(const IndexType i, T& out) {
    PERFCTR_INCREMENT(readCount);
    Assert(i < N);
    data[i].Decrypt(out);
  }
};  // struct MemServerBackend

template <typename T, uint64_t MEM_SIZE = (3UL << 20)>
using MemServer =
    _ORAM::Cached<T, MemServerBackend<T, MEM_SIZE>, ORAM_SERVER__CACHE_SIZE, 2>;
}  // namespace MemServer
}  // namespace EM