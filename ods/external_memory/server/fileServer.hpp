#pragma once
#include <cinttypes>

#include "common/lrucache.hpp"
#include "external_memory/server/cached.hpp"

#define DEFAULT_FILENAME "./storage.bin"
// This has the same interface as an array that transfers data from a file
//
namespace EM {
namespace FileServer {
// T is typically a LargeBucket
//
template <typename T>
struct FileServerBackend {
  static inline constexpr auto sizeOfT = sizeof(typename T::Encrypted_t);
  typedef uint64_t IndexType;
  IndexType N;
  std::fstream lbios;

  FileServerBackend(IndexType _N) : N(_N) {
    PROFILE_F();

    lbios.open(DEFAULT_FILENAME, std::fstream::in | std::fstream::out |
                                     std::fstream::binary |
                                     std::fstream::trunc);
    Assert(lbios.is_open());
    for (IndexType i = 0; i < N; i++) {
      uint8_t zeros[sizeOfT] = {0};
      memset(zeros, 0, sizeOfT);
      lbios.write((char*)zeros, sizeOfT);
    }
    Assert(lbios.is_open());
    lbios.flush();
    X_LOG("ORAMClientInterface: Done allocating file (", lbios.tellp(),
          " bytes written)");
  }

  void Write(const IndexType i, const T& in) {
    Assert(i < N);
    std::streampos filePos = i * sizeOfT;
    // UNDONE: encrypt directly to file
    //
    typename T::Encrypted_t inEnc;
    inEnc.Encrypt(in);
    lbios.seekp(filePos);
    lbios.write((char*)&inEnc, sizeOfT);
    // lbios.flush();
  }

  void Read(const IndexType i, T& out) {
    Assert(i < N);
    std::streampos filePos = i * sizeOfT;
    // UNDONE: decrypt directly from file
    //
    typename T::Encrypted_t outEnc;
    lbios.seekg(filePos);
    lbios.read((char*)&outEnc, sizeOfT);
    outEnc.Decrypt(out);
    // lbios.flush();
  }
};  // struct FileServerBackend
template <typename T, uint64_t CACHE_SIZE = (1UL << 26)>
using FileServer = _ORAM::Cached<T, FileServerBackend<T>, CACHE_SIZE, 2>;
}  // namespace FileServer
}  // namespace EM