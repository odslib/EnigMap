
#pragma once

// In this file we define a vector that is partially in internal memory
// and partially in external memory. The interface for this vector is similar
// to the interface for stl's vectors and for ORAM, with the exception that
// indexes on this vector might have to be read from main memory rather than
// external memory.
// This vector will leak the position of every access.
// This vector can leak it's content or not, based on the parameter ENCRYPTED.
// The encrypted also prevents from knowing if modifications were made
// to the external memory, as blocks will be encrypted using AES GCM with a
// random new IV every time they are commited to external memory. We provide
// functions to assert that no external memory accesses were done
//

#include <vector>
namespace EM::Vector {
template <typename T
          // T's per page: the number of T's that should go into a page of
          // external memory. UNDONE: server -> external memory server type.
          // ENCRYPTED: store T encrypted or as plaintext.
          >
struct Vector {
  // UNDONE(): actually implement the external fetches, rathan than this
  // std::vector<>.
  //
  std::vector<T> mem;

  // UNDONE: fileserver here.
  //
  uint64_t N;

  // https://www.internalpointers.com/post/writing-custom-iterators-modern-cpp
  struct Iterator {
    // Iterator tags here...
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = typename std::vector<T>::iterator;
    using reference = T&;

    // Iterator constructors here...
    explicit Iterator(pointer ptr) : m_ptr(ptr) {}

    reference operator*() const { return *m_ptr; }

    pointer operator->() { return m_ptr; }

    const pointer operator->() const { return m_ptr; }

    // Prefix increment
    Iterator& operator++() {
      ++m_ptr;
      return *this;
    }

    // Prefix decrement
    Iterator& operator--() {
      --m_ptr;
      return *this;
    }

    // Prefix increment
    Iterator& operator+=(int n) {
      m_ptr += n;
      return *this;
    }

    // Postfix increment
    Iterator operator++(int) {
      Iterator tmp = *this;
      ++(*this);
      return tmp;
    }

    friend bool operator==(const Iterator& a, const Iterator& b) {
      return a.m_ptr == b.m_ptr;
    };
    friend bool operator!=(const Iterator& a, const Iterator& b) {
      return a.m_ptr != b.m_ptr;
    };
    friend bool operator<(const Iterator& a, const Iterator& b) {
      return a.m_ptr < b.m_ptr;
    };
    friend bool operator<=(const Iterator& a, const Iterator& b) {
      return a.m_ptr <= b.m_ptr;
    };
    friend size_t operator-(const Iterator& it1, const Iterator& it2) {
      return it1.m_ptr - it2.m_ptr;
    }
    friend Iterator operator+(const Iterator& it, size_t size) {
      return Iterator(it.m_ptr + size);
    }
    friend Iterator operator-(const Iterator& it, size_t size) {
      return Iterator(it.m_ptr - size);
    }

   private:
    pointer m_ptr;
  };

  // default:
  explicit Vector(uint64_t N_ = 0) : N(N_), mem(N_) {}
  // fill:
  Vector(uint64_t N_, const T& default_val) : N(N_), mem(N_, default_val) {}
  // UNDONE: range and copy.

  T& At(uint64_t index) {
    Assert(index < N);
    return mem[index];
  }
  const T& At(uint64_t index) const {
    Assert(index < N);
    return mem[index];
  }

  T& operator[](uint64_t index) { return At(index); }
  const T& operator[](uint64_t index) const { return At(index); }

  uint64_t size() const { return N; }

  Iterator begin() { return Iterator(mem.begin()); }

  Iterator end() { return Iterator(mem.end()); }

  void LRUTouch(uint64_t index) {
    // UNDONE(): touches the given index in the LRU cache.
    //
  }

#ifndef NDEBUG
  // UNDONE(): don't forget to add these assertions to the cache once
  // implemented.
  //
  bool assertOnCache = false;
  void SetAssertOnCache(bool enable = true) { assertOnCache = enable; }
#endif
};
}  // namespace EM::Vector