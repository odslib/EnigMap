
#pragma once

#include "external_memory/server/batchFrontend.hpp"
#include "external_memory/server/memServer.hpp"
#include "external_memory/server/serverFrontend.hpp"

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

#include "common/encrypted.hpp"
namespace EM::ExtVector {
template <typename T,
          // UNDONE: parametrize here the backend / frontend
          uint64_t page_size = std::max((1UL << 12) - 32, sizeof(T)),  // 4kB
          bool ENCRYPTED = true, bool AUTH = true,
          uint64_t cache_size = ORAM_SERVER__CACHE_SIZE

          // T's per page: the number of T's that should go into a page of
          // external memory. UNDONE: server -> external memory server type.
          // ENCRYPTED: store T encrypted or as plaintext.
          >
struct Vector {
  // UNDONE(): actually implement the external fetches, rathan than this
  // std::vector<>.
  //

  static constexpr uint64_t page_count = page_size / sizeof(T);

  struct Dim2D {
    size_t row;
    size_t col;
    size_t unit;
  };
  std::vector<Dim2D> transposeData;

  struct Page {
    T pages[page_count];
    using Encrypted_t = std::conditional_t<
        AUTH, FreshEncrypted<Page>,
        std::conditional_t<ENCRYPTED, Encrypted<Page>, NonEncrypted<Page>>>;
  };

  static constexpr uint64_t DMCacheSize =
      cache_size != GetNextPowerOfTwo(cache_size)
          ? cache_size
          : std::max(1UL, cache_size - 1);

  using Server = EM::MemoryServer::ServerFrontendInstance<
      Page, ::EM::Backend::MemServerBackend, ENCRYPTED, AUTH, DMCacheSize>;

  Server server;

  // UNDONE: fileserver here.
  //
  uint64_t N;

  // https://www.internalpointers.com/post/writing-custom-iterators-modern-cpp
  struct Iterator {
    // Iterator tags here...
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = uint64_t;
    using page_idx_type = uint64_t;
    using page_offset_type = uint64_t;
    using reference = T&;
    using const_reference = const T&;

    // Iterator constructors here...
    explicit Iterator(pointer ptr, Vector& vec) : m_ptr(ptr), vec_ptr(&vec) {}

    Iterator() : m_ptr(0), vec_ptr(NULL) {}

    reference operator*() {
      Assert(m_ptr < vec_ptr->N);
      // const size_t realIdx = vec_ptr->mapTransposeIdx(m_ptr);
      const size_t realIdx = m_ptr;
      const size_t pageIdx = realIdx / page_count;
      const size_t pageOffset = realIdx % page_count;
      return vec_ptr->server.Access(pageIdx).pages[pageOffset];
    }

    const_reference operator*() const {
      Assert(m_ptr < vec_ptr->N);
      // const size_t realIdx = vec_ptr->mapTransposeIdx(m_ptr);
      const size_t realIdx = m_ptr;
      const size_t pageIdx = realIdx / page_count;
      const size_t pageOffset = realIdx % page_count;
      return vec_ptr->server.AccessReadOnly(pageIdx).pages[pageOffset];
    }

    // don't write back the page
    const_reference derefNoWriteBack() const {
      Assert(m_ptr < vec_ptr->N);
      // const size_t realIdx = vec_ptr->mapTransposeIdx(m_ptr);
      const size_t realIdx = m_ptr;
      const size_t pageIdx = realIdx / page_count;
      const size_t pageOffset = realIdx % page_count;
      return vec_ptr->server.AccessNoWriteBack(pageIdx).pages[pageOffset];
    }

    // always skip read page
    reference derefWriteOnly() {
      Assert(m_ptr < vec_ptr->N);
      // const size_t realIdx = vec_ptr->mapTransposeIdx(m_ptr);
      const size_t realIdx = m_ptr;
      const size_t pageIdx = realIdx / page_count;
      const size_t pageOffset = realIdx % page_count;
      return vec_ptr->server.AccessWriteOnly(pageIdx).pages[pageOffset];
    }

    T* operator->() { return &(**this); }

    const T* operator->() const { return &(**this); }

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

    page_idx_type get_page_idx() const { return m_ptr / page_count; }

    page_offset_type get_page_offset() const { return m_ptr % page_count; }

    auto& getVector() { return *vec_ptr; }

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
      return Iterator(it.m_ptr + size, *it.vec_ptr);
    }
    friend Iterator operator-(const Iterator& it, size_t size) {
      return Iterator(it.m_ptr - size, *it.vec_ptr);
    }

   private:
    Vector* vec_ptr;
    pointer m_ptr;
    // Server* server;
  };

  // default:
  explicit Vector(uint64_t N_ = 0, typename Server::BackendType& _backend =
                                       *Backend::g_DefaultBackend)
      : N(N_), server(_backend, N_ / page_count + 1, makeDefaultPage()) {}
  // explicit Vector(uint64_t N_=0, typename Server::BackendType&
  // _backend=*Backend::g_DefaultBackend) : N(N_), server(_backend,
  // N_/page_count+1) { } fill:
  Page makeDefaultPage() {
    // UNDONE(): Make T::DUMMY() a template function so it supports
    // DUMMY<int>().
    //
    T tdummy;
    return makeDefaultPage(tdummy);
  }
  Page makeDefaultPage(const T& defaultVal) {
    Page defaultPage;
    std::fill_n(defaultPage.pages, page_count, defaultVal);
    return defaultPage;
  }
  Vector(uint64_t N_, const T& defaultVal,
         typename Server::BackendType& _backend = *Backend::g_DefaultBackend)
      : N(N_),
        server(_backend, N_ / page_count + 1, makeDefaultPage(defaultVal)) {}
  // UNDONE: range and copy.

  Vector(Iterator begin, Iterator end,
         typename Server::BackendType& _backend = *Backend::g_DefaultBackend)
      : N(end - begin),
        server(_backend, (end - begin) / page_count + 1, makeDefaultPage()) {
    auto outIt = Iterator(0, *this);
    std::copy(begin, end, outIt);
  }

  T& AtForLateInit(uint64_t index) {
    return Iterator(index, *this).derefWriteOnly();
  }

  T& At(uint64_t index) { return *Iterator(index, *this); }
  const T& At(uint64_t index) const { return *Iterator(index, *this); }

  T& operator[](uint64_t index) { return At(index); }
  const T& operator[](uint64_t index) const { return At(index); }

  uint64_t size() const { return N; }

  Iterator begin() { return Iterator(0, *this); }

  Iterator end() { return Iterator(N, *this); }

  void LRUTouch(uint64_t index) {
    // UNDONE(): touches the given index in the LRU cache.
    //
  }

  void logicalTranspose(size_t row, size_t col) {
#ifndef NDEBUG
    Assert(N % (row * col) == 0);
#endif
    std::reverse(transposeData.begin(), transposeData.end());
    transposeData.push_back({row, col, N / (row * col)});
    std::reverse(transposeData.begin(), transposeData.end());
  }

  inline size_t mapTransposeIdx(size_t idx) const {
    PROFILE_F();
    for (const Dim2D& dim : transposeData) {
      const size_t bucketIdx = idx / dim.unit;
      const size_t offset = idx % dim.unit;
      const size_t rowIdx = bucketIdx / dim.col;
      const size_t colIdx = bucketIdx % dim.col;
      const size_t newBucketIdx = colIdx * dim.row + rowIdx;
      idx = newBucketIdx * dim.unit + offset;
    }

    return idx;
  }

  struct Reader {
    Iterator it;
    Iterator end;
    Reader(Iterator _begin, Iterator _end) : it(_begin), end(_end) {}

    const T& get() {
      const auto& it_const = it;
      return *it_const;
    }

    const T& read() {
      const T& val = get();
      ++it;
      return val;
    }
    bool eof() { return end <= it; }
  };

  struct Writer {
    Iterator it;
    Iterator end;
    Writer() {}
    Writer(Iterator _begin, Iterator _end) { init(_begin, _end); }
    void init(Iterator _begin, Iterator _end) {
      this->it = _begin;
      this->end = _end;
    }
    void write(const T& element) {
      *it = element;
      ++it;
    }
    bool eof() { return end <= it; }
    void flush() {}
  };

#ifndef NDEBUG
  // UNDONE(): don't forget to add these assertions to the cache once
  // implemented.
  //
  bool assertOnCache = false;
  void SetAssertOnCache(bool enable = true) { assertOnCache = enable; }
#endif
};

template <class InputIterator, class OutputIterator>
static OutputIterator Copy(InputIterator begin, InputIterator end,
                           OutputIterator to) {
  for (auto it = begin; it != end; ++it, ++to) {
    const auto& it_const = it;
    *to = *it_const;
  }
  return to;
}

template <class InputIterator, class OutputIterator>
static OutputIterator CopyWithoutWriteBackInput(InputIterator begin,
                                                InputIterator end,
                                                OutputIterator to) {
  for (auto it = begin; it != end; ++it, ++to) {
    *to = it.derefNoWriteBack();
  }
  return to;
}

template <class InputIterator, class OutputIterator>
static OutputIterator CopyForLateInit(InputIterator begin, InputIterator end,
                                      OutputIterator to) {
  auto it = begin;
  auto toEnd = to + (end - begin);
  auto fullPageEnd = toEnd - toEnd.get_page_offset();
  if (to < fullPageEnd) {
    for (; to.get_page_offset() != 0; ++it, ++to) {
      const auto& it_const = it;
      *to = *it_const;
    }
    for (; to != fullPageEnd; ++it, ++to) {
      const auto& it_const = it;
      to.derefWriteOnly() = *it_const;
    }
  }

  for (; to != toEnd; ++it, ++to) {
    const auto& it_const = it;
    *to = *it_const;
  }
  return to;
}

template <class Iterator, typename T>
static void FillForLateInit(Iterator begin, Iterator end, const T& val) {
  auto to = begin;
  auto fullPageEnd = end - end.get_page_offset();
  if (to < fullPageEnd) {
    for (; to.get_page_offset() != 0; ++to) {
      *to = val;
    }
    for (; to != fullPageEnd; ++to) {
      to.derefWriteOnly() = val;
    }
  }
  for (; to != end; ++to) {
    *to = val;
  }
}
}  // namespace EM::ExtVector