
#pragma once

#include "external_memory/server/memServer.hpp"
// In this file we define a vector that is partially in internal memory
// and partially in external memory. The interface for this vector is similar
// to the interface for stl's vectors and for ORAM, with the exception that 
// indexes on this vector might have to be read from main memory rather than
// external memory.
// This vector will leak the position of every access.
// This vector can leak it's content or not, based on the parameter ENCRYPTED.
// The encrypted also prevents from knowing if modifications were made 
// to the external memory, as blocks will be encrypted using AES GCM with a random 
// new IV every time they are commited to external memory.
// We provide functions to assert that no external memory accesses were done
//

#include <vector>
namespace EM::Vector {
  template<
    typename T
    // T's per page: the number of T's that should go into a page of external memory.
    // UNDONE: server -> external memory server type.
    // ENCRYPTED: store T encrypted or as plaintext.
  >
  struct Vector {
    // UNDONE(): actually implement the external fetches, rathan than this std::vector<>.
    //
    std::vector<T> mem;

    // UNDONE: fileserver here.
    //    
    uint64_t N;

    // default:
    Vector(uint64_t N_=0) : N(N_), mem(N_) { }
    // fill:
    Vector(uint64_t N_, const T& default_val) : N(N_), mem(N_, default_val) { }
    // UNDONE: range and copy.

    T& At(uint64_t index) {
      Assert(index < N);
      return mem[index];
    }
    const T& At(uint64_t index) const {
      Assert(index < N);
      return mem[index];
    }

    T& operator[](uint64_t index) {
      return At(index);
    }
    const T& operator[](uint64_t index) const {
      return At(index);
    }

    uint64_t size() const { return N; }

    void LRUTouch(uint64_t index) {
      // UNDONE(): touches the given index in the LRU cache.
      //
    }

    #ifndef NDEBUG
    // UNDONE(): don't forget to add these assertions to the cache once implemented.
    //
    bool assertOnCache = false;
    void SetAssertOnCache(bool enable=true) {
      assertOnCache = enable;
    }
    #endif
  };
}