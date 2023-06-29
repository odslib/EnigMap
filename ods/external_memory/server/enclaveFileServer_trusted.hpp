#pragma once
#include <cinttypes>

#include "common/lrucache.hpp"
#include "external_memory/server/cached.hpp"

// Similar to fileServer.hpp, but separated to be used within an enclave.
//
namespace EM {
namespace EnclaveFileServer {
// T is typically a LargeBucket
//
template <typename T>
struct EnclaveFileServerBackend {
  static inline constexpr auto sizeOfT = sizeof(typename T::Encrypted_t);
  typedef uint64_t IndexType;
  IndexType N;

  EnclaveFileServerBackend(IndexType _N) : N(_N) {
    TRACE_FUNCTION(_N);
    ocall_InitServer(nullptr, 4096, _N);
  }

  void Swap(const IndexType index_in, const IndexType index_out, const T& in,
            T& out) {
    Write(index_in, in);
    Read(index_out, out);
  }

  void Write(const IndexType i, const T& in) {
    Assert(i < N);
    Assert(sizeOfT <= 4096);

    uint8_t aux[4096] = {0};
    typename T::Encrypted_t inEnc;
    inEnc.Encrypt(in);
    memcpy(aux, &inEnc, sizeOfT);
    printf("prepare enclave write page %ld of size %ld\n", i, sizeOfT);
    // for (size_t i = 0; i < 4096; ++i) {
    //   printf("%d ", aux[i]);
    // }
    // printf("\n");
    ocall_Write(i, 4096, aux);
    printf("end enclave write page %ld of size %ld\n", i, sizeOfT);
  }

  void Read(const IndexType i, T& out) {
    Assert(i < N);
    Assert(sizeOfT <= 4096);

    uint8_t aux[4096];
    typename T::Encrypted_t inEnc;
    printf("prepare enclave read page %ld of size %ld\n", i, sizeOfT);
    ocall_Read(i, 4096, aux);
    printf("end enclave read page %ld of size %ld\n", i, sizeOfT);
    memcpy(&inEnc, aux, sizeOfT);
    inEnc.Decrypt(out);
  }
};  // struct EnclaveFileServerBackend

template <typename T, uint64_t CACHE_SIZE = (1UL << 16)>
using EnclaveFileServer =
    _ORAM::Cached<T, EnclaveFileServerBackend<T>, CACHE_SIZE, 2>;

}  // namespace EnclaveFileServer
}  // namespace EM