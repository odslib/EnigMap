
#pragma once

#include "external_memory/server/memServer.hpp"
#ifdef DISK_IO
#include "external_memory/server/batchFrontend.hpp"
#else
#include "external_memory/server/serverFrontend.hpp"
#endif

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
namespace EM::NonCachedVector {
template <typename T,
          // UNDONE: parametrize here the backend / frontend
          uint64_t page_size =
              std::max((1UL << 12) - 32, 4 * sizeof(T)),  // 4kB
          bool ENCRYPTED = true, bool AUTH = true

          // T's per page: the number of T's that should go into a page of
          // external memory. UNDONE: server -> external memory server type.
          // ENCRYPTED: store T encrypted or as plaintext.
          >
struct Vector {
  // UNDONE(): actually implement the external fetches, rathan than this
  // std::vector<>.
  //
  static constexpr bool encrypted = ENCRYPTED;
  static constexpr bool auth = AUTH;
  static constexpr uint64_t item_per_page =
      GetNextPowerOfTwo(page_size / sizeof(T) / 2 + 1);

  uint64_t referenceCount;

  struct Page {
    T pages[item_per_page];
    //   uint8_t padding[page_size - item_per_page * sizeof(T)];
    using Encrypted_t = std::conditional_t<
        AUTH, FreshEncrypted<Page>,
        std::conditional_t<ENCRYPTED, Encrypted<Page>, NonEncrypted<Page>>>;
  };

  Page cachePage;

  Page& getPageInstance() { return cachePage; }

// using Server =
// EM::MemoryServer::NonCachedServerFrontendInstance<Page,::EM::Backend::MemServerBackend,ENCRYPTED>;
#ifdef DISK_IO
  using Server = EM::MemoryServer::NonCachedBatchServerFrontend<
      Page, ::EM::Backend::MemServerBackend, ENCRYPTED, AUTH>;
#else
  using Server = EM::MemoryServer::NonCachedServerFrontendInstance<
      Page, ::EM::Backend::MemServerBackend, ENCRYPTED, AUTH>;
#endif

  // fileserver here.
  //
  uint64_t N;
  Server server;

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

    // should not be called unless absolutely necessary
    T operator*() {

      const size_t pageIdx = m_ptr / item_per_page;
      const size_t pageOffset = m_ptr % item_per_page;
      Page page;
      if constexpr (AUTH) {
        vec_ptr->server.Read(pageIdx, page, 0);
      } else {
        vec_ptr->server.Read(pageIdx, page);
      }
      return page.pages[pageOffset];
    }

    const T* operator->() { return &(**this); }

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

    page_idx_type get_page_idx() const { return m_ptr / item_per_page; }

    page_offset_type get_page_offset() const { return m_ptr % item_per_page; }

    pointer get_m_ptr() const { return m_ptr; }

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
    pointer m_ptr;
    Vector* vec_ptr;
  };

  struct Reader {
    Page cache;
    Iterator it;
    Iterator end;
    T* curr = (T*)UINT64_MAX;
    uint32_t counter;
    Reader(Iterator _begin, Iterator _end, uint32_t counter = 0)
        : it(_begin), end(_end), counter(counter) {}

    T& get() {
      Assert(!eof());
      if (curr >= cache.pages + item_per_page) {
        auto& vec = it.getVector();
        const size_t pageIdx = it.get_page_idx();
        if constexpr (AUTH) {
          vec.server.Read(pageIdx, cache, counter);
        } else {
          vec.server.Read(pageIdx, cache);
        }

        if (curr != cache.pages + item_per_page) {
          curr = cache.pages + it.get_page_offset();
        } else {
          curr = cache.pages;
        }
      }
      return *curr;
    }

    const T& read() {
      const T& val = get();
      ++it;
      ++curr;
      return val;
    }
    bool eof() { return end <= it; }
  };

#ifdef DISK_IO
  struct PrefetchReader {
    static const uint64_t cacheCapacity = Server::bufferCapacity;
    Page cache[cacheCapacity];
    Iterator it;
    Iterator end;
    T* curr = (T*)UINT64_MAX;
    const size_t endPageIdx;
    uint32_t counter;
    PrefetchReader(Iterator _begin, Iterator _end, uint32_t counter = 0)
        : it(_begin),
          end(_end),
          endPageIdx((_end - 1).get_page_idx() + 1),
          counter(counter) {}

    T& get() {
      Assert(!eof());
      if (curr >= (T*)(cache + cacheCapacity)) {
        auto& vec = it.getVector();
        const size_t pageIdx = it.get_page_idx();
        size_t cachePageBeginIdx = pageIdx;
        size_t cachePageEndIdx =
            std::min(cachePageBeginIdx + cacheCapacity, endPageIdx);
        for (size_t i = cachePageBeginIdx, cacheIdx = 0; i < cachePageEndIdx;
             ++i, ++cacheIdx) {
          if constexpr (AUTH) {
            vec.server.ReadLazy(i, cache[cacheIdx], counter);
          } else {
            vec.server.ReadLazy(i, cache[cacheIdx]);
          }
        }
        vec.server.flushRead();
        if (curr != (T*)(cache + cacheCapacity)) {
          curr = cache[0].pages + it.get_page_offset();
        } else {
          curr = (T*)cache;
        }
      }
      return *curr;
    }

    const T& read() {
      Assert(!eof());
      const T& val = get();
      ++it;
      ++curr;
      return val;
    }
    bool eof() { return end <= it; }
  };
#else
  struct PrefetchReader : Reader {
    PrefetchReader(Iterator _begin, Iterator _end, uint32_t counter = 0)
        : Reader(_begin, _end, counter) {}
  };
#endif

#ifdef DISK_IO
  // when there are multiple readers for the same frontend, e.g., external merge
  // sort. the reader always send read requests for all its cache pages in
  // advanced, and check if the job is done by the time it read those pages the
  // advantage is that the requests for multipler readers can be interleaved and
  // the most on-demand pages can be handled first
  struct LazyPrefetchReader {
    std::vector<Page> cache;  // ring cache

    std::vector<uint64_t> jobCounters;
    const uint64_t cacheCapacity;
    Iterator it;
    Iterator end;
    T* curr = (T*)UINT64_MAX;
    T* currPageEnd = (T*)UINT64_MAX;
    const size_t endPageIdx;
    uint32_t counter;
#ifndef NDEBUG
    bool hasInit = false;
#endif
    LazyPrefetchReader(Iterator _begin, Iterator _end, uint32_t counter = 0,
                       size_t cacheCapacity = 2)
        : it(_begin),
          end(_end),
          cacheCapacity(cacheCapacity),
          endPageIdx((_end - 1).get_page_idx() + 1),
          counter(counter) {
      Assert(cacheCapacity >= 2);
    }

    void init() {
#ifndef NDEBUG
      Assert(!hasInit);
      hasInit = true;
#endif
      cache.resize(cacheCapacity);
      jobCounters.resize(cacheCapacity);
      curr = cache[0].pages + it.get_page_offset();
      currPageEnd = cache[1].pages;
      auto& vec = it.getVector();
      for (size_t cacheIdx = 0, pageIdx = it.get_page_idx();
           cacheIdx < cacheCapacity && pageIdx < endPageIdx;
           ++cacheIdx, ++pageIdx) {
        // printf("read %ld\n", pageIdx);
        if constexpr (AUTH) {
          jobCounters[cacheIdx] =
              vec.server.ReadLazy(pageIdx, cache[cacheIdx], counter);
        } else {
          jobCounters[cacheIdx] = vec.server.ReadLazy(pageIdx, cache[cacheIdx]);
        }
      }
      vec.server.flushRead();
    }

    T& get() {
      Assert(hasInit);
      Assert(!eof());
      if (curr == currPageEnd) {
        size_t prevCacheIdx = (Page*)currPageEnd - &cache[0] - 1;
        size_t currCacheIdx = (prevCacheIdx + 1) % cacheCapacity;
        auto& vec = it.getVector();
        if (!vec.server.isJobDone(jobCounters[currCacheIdx])) {
          vec.server.flushRead();
          Assert(vec.server.isJobDone(jobCounters[currCacheIdx]));
        }
        size_t nextPageIdx = it.get_page_idx() + cacheCapacity - 1;
        if (nextPageIdx < endPageIdx) {
          // overwrite the previously cached page
          // printf("read %ld\n", nextPageIdx);
          if constexpr (AUTH) {
            jobCounters[prevCacheIdx] =
                vec.server.ReadLazy(nextPageIdx, cache[prevCacheIdx], counter);
          } else {
            jobCounters[prevCacheIdx] =
                vec.server.ReadLazy(nextPageIdx, cache[prevCacheIdx]);
          }
        }
        curr = (T*)(&cache[0] + currCacheIdx);  // at the begining of the page
        currPageEnd = curr + item_per_page;
      }
      return *curr;
    }

    const T& read() {
      Assert(!eof());
      const T& val = get();
      ++it;
      ++curr;
      return val;
    }
    bool eof() { return end <= it; }
  };
#else
  struct LazyPrefetchReader : Reader {
    LazyPrefetchReader(Iterator _begin, Iterator _end, uint32_t counter = 0)
        : Reader(_begin, _end, counter) {}
    void init() {}
  };
#endif

  // writer always overwrites data between begin and end
  struct Writer {
    Page cache;
    Iterator it;
    Iterator end;
    T* curr;
    uint32_t counter;
    Writer() {}
    Writer(Iterator _begin, Iterator _end, uint32_t counter = 0)
        : counter(counter) {
      init(_begin, _end, counter);
    }
    void init(Iterator _begin, Iterator _end, uint32_t counter = 0) {
      this->counter = counter;
      Assert(_begin.getVector().begin() <= _begin);
      Assert(_end <= _end.getVector().end());
      this->it = _begin;
      this->end = _end;
      auto& vec = it.getVector();
      size_t pageIdx = it.get_page_idx();
      size_t pageOffset = it.get_page_offset();
      if (pageOffset != 0) {
        if constexpr (AUTH) {
          vec.server.Read(pageIdx, cache, counter);
        } else {
          vec.server.Read(pageIdx, cache);
        }
      }
      curr = cache.pages + pageOffset;
    }
    void write(const T& element) {
      *curr = element;
      ++curr;
      ++it;
      if (curr == cache.pages + item_per_page) {
        auto& vec = it.getVector();
        Assert(it.get_page_offset() == 0);
        size_t pageIdx = it.get_page_idx() - 1;
        if constexpr (AUTH) {
          vec.server.WriteLazy(pageIdx, cache, counter);
        } else {
          vec.server.WriteLazy(pageIdx, cache);
        }
        curr = cache.pages;
      }
    }
    bool eof() { return end <= it; }
    void flush() {
      size_t pageOffset = it.get_page_offset();
      size_t pageIdx = it.get_page_idx();
      auto& vec = it.getVector();
      if (pageOffset != 0) {
        Page originalPage;

        if constexpr (AUTH) {
          if (counter == 0) {
            vec.server.Read(pageIdx, originalPage, counter);
          }
          // for counter greater than 0, we have to overwrite the last page.
          // this issue should be addressed when calling the writer
        } else {
          vec.server.Read(pageIdx, originalPage);
        }

        std::memcpy(originalPage.pages, cache.pages, pageOffset * sizeof(T));
        if constexpr (AUTH) {
          vec.server.WriteLazy(pageIdx, originalPage, counter);
        } else {
          vec.server.WriteLazy(pageIdx, originalPage);
        }
      }
      vec.server.flushWrite();
    }
  };

  // default:
  explicit Vector(uint64_t N_ = 0, typename Server::BackendType& _backend =
                                       *Backend::g_DefaultBackend)
      : N(N_), server(_backend, N_ / item_per_page + 1) {}
  // fill:

  Page makeDefaultPage(const T& defaultVal) {
    Page defaultPage;
    std::fill_n(defaultPage.pages, item_per_page, defaultVal);
    return defaultPage;
  }
  Vector(uint64_t N_, const T& defaultVal,
         typename Server::BackendType& _backend = *Backend::g_DefaultBackend)
      : N(N_),
        server(_backend, N_ / item_per_page + 1, makeDefaultPage(defaultVal)) {}
  // UNDONE: range and copy.

  Vector(Vector&& vec)
      : server(vec.server), referenceCount(vec.referenceCount), N(vec.N) {
    vec.server.preventFree();
  }

  Vector(Vector& vec) = delete;

  template <class InputIterator>
  Vector(InputIterator inputBegin, InputIterator inputEnd,
         typename Server::BackendType& _backend = *Backend::g_DefaultBackend)
      : N(inputEnd - inputBegin),
        server(_backend, (inputEnd - inputBegin) / item_per_page + 1) {
    CopyOut(inputBegin, inputEnd, begin());
  }

  // these functions are slow, as they need to read the whole page
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
};

template <class InputIterator, class OutputIterator>
static OutputIterator CopyOut(InputIterator begin, InputIterator end,
                              OutputIterator to, uint32_t counter = 0) {
  using T = typename std::iterator_traits<OutputIterator>::value_type;
  auto& vec = to.getVector();
  using Vec = std::remove_reference_t<decltype(vec)>;
  auto& page = vec.getPageInstance();
  constexpr bool auth = Vec::auth;
  using Page = typename Vec::Page;
  size_t inputSize = end - begin;
  static size_t item_per_page = Vec::item_per_page;
  size_t toBeginOffset = to.get_page_offset();
  size_t toBeginPageIdx = to.get_page_idx();
  auto toEnd = to + inputSize;
  Assert(to.getVector().begin() <= to);
  Assert(toEnd <= toEnd.getVector().end());
  size_t toEndOffset = toEnd.get_page_offset();
  size_t toEndPageIdx = toEnd.get_page_idx();

  auto from = begin;

  // if constexpr (auth) {
  //     if (toBeginOffset || toEndOffset) {
  //         dbg_printf("Copy out not aligned\n");
  //     }
  // }

  if (toBeginPageIdx == toEndPageIdx) {  // within a single page
    if (inputSize) {
      if constexpr (auth) {
        Assert(counter == 0);  // fine if counter is always 0
        vec.server.Read(toBeginPageIdx, page, counter);
      } else {
        vec.server.Read(toBeginPageIdx, page);
      }
      for (size_t pageOffset = toBeginOffset; pageOffset != toEndOffset;
           ++pageOffset) {
        page.pages[pageOffset] = *from;
        ++from;
      }
      if constexpr (auth) {
        vec.server.Write(toBeginPageIdx, page, counter);
      } else {
        vec.server.Write(toBeginPageIdx, page);
      }
    }
    return toEnd;
  }

  if (toBeginOffset != 0) {
    if constexpr (auth) {
      Assert(counter == 0);  // fine if counter is always 0
      vec.server.Read(toBeginPageIdx, page, counter);
    } else {
      vec.server.Read(toBeginPageIdx, page);
    }
    for (size_t pageOffset = toBeginOffset; pageOffset != item_per_page;
         ++pageOffset) {
      page.pages[pageOffset] = *from;
      ++from;
    }
    if constexpr (auth) {
      vec.server.Write(toBeginPageIdx++, page, counter);
    } else {
      vec.server.Write(toBeginPageIdx++, page);
    }
  }
  for (size_t pageIdx = toBeginPageIdx; pageIdx < toEndPageIdx; ++pageIdx) {
    if constexpr (auth) {
      vec.server.WriteLazy(pageIdx, *(Page*)(&(*from)), counter);
    } else {
      vec.server.WriteLazy(pageIdx, *(Page*)(&(*from)));
    }

    from += item_per_page;
  }
  vec.server.flushWrite();
  if (toEndOffset != 0) {
    if constexpr (auth) {
      Assert(counter == 0);
      vec.server.Read(toEndPageIdx, page, counter);
    } else {
      vec.server.Read(toEndPageIdx, page);
    }
    for (size_t pageOffset = 0; pageOffset != toEndOffset; ++pageOffset) {
      page.pages[pageOffset] = *from;
      ++from;
    }
    if constexpr (auth) {
      vec.server.Write(toEndPageIdx, page, counter);
    } else {
      vec.server.Write(toEndPageIdx, page);
    }
  }
  Assert(from == end);
  return toEnd;
}

template <typename T, class OutputIterator>
static void Fill(OutputIterator begin, OutputIterator end, const T& val,
                 uint32_t counter = 0) {
  auto& vec = begin.getVector();
  using Vec = std::remove_reference_t<decltype(vec)>;
  auto& page = vec.getPageInstance();
  constexpr bool auth = Vec::auth;
  using Page = typename Vec::Page;

  size_t inputSize = end - begin;
  static size_t item_per_page = Vec::item_per_page;
  size_t beginOffset = begin.get_page_offset();
  size_t beginPageIdx = begin.get_page_idx();
  size_t endOffset = end.get_page_offset();
  size_t endPageIdx = end.get_page_idx();

  // if constexpr (auth) {
  //     if (beginOffset || endOffset) {
  //         dbg_printf("Fill not aligned\n");
  //     }
  // }

  if (beginPageIdx == endPageIdx) {  // within a single page
    if (inputSize) {
      if constexpr (auth) {
        Assert(counter == 0);  // fine if counter is always 0
        vec.server.Read(beginPageIdx, page, counter);
      } else {
        vec.server.Read(beginPageIdx, page);
      }
      for (size_t pageOffset = beginOffset; pageOffset != endOffset;
           ++pageOffset) {
        page.pages[pageOffset] = val;
      }
      if constexpr (auth) {
        vec.server.Write(beginPageIdx, page, counter);
      } else {
        vec.server.Write(beginPageIdx, page);
      }
    }
    return;
  }

  if (beginOffset != 0) {
    if constexpr (auth) {
      Assert(counter == 0);
      vec.server.Read(beginPageIdx, page, counter);
    } else {
      vec.server.Read(beginPageIdx, page);
    }
    for (size_t pageOffset = beginOffset; pageOffset != item_per_page;
         ++pageOffset) {
      page.pages[pageOffset] = val;
    }
    if constexpr (auth) {
      vec.server.Write(beginPageIdx++, page, counter);
    } else {
      vec.server.Write(beginPageIdx++, page);
    }
  }
  for (size_t pageOffset = 0; pageOffset != item_per_page; ++pageOffset) {
    page.pages[pageOffset] = val;
  }
  for (size_t pageIdx = beginPageIdx; pageIdx < endPageIdx; ++pageIdx) {
    if constexpr (auth) {
      vec.server.WriteLazy(pageIdx, page, counter);
    } else {
      vec.server.WriteLazy(pageIdx, page);
    }
  }
  vec.server.flushWrite();
  if (endOffset != 0) {
    if constexpr (auth) {
      Assert(counter == 0);
      vec.server.Read(endPageIdx, page, counter);
    } else {
      vec.server.Read(endPageIdx, page);
    }
    for (size_t pageOffset = 0; pageOffset != endOffset; ++pageOffset) {
      page.pages[pageOffset] = val;
    }
    if constexpr (auth) {
      vec.server.Write(endPageIdx, page, counter);
    } else {
      vec.server.Write(endPageIdx, page);
    }
  }
}

template <class InputIterator, class OutputIterator>
static OutputIterator CopyIn(InputIterator begin, InputIterator end,
                             OutputIterator to, uint32_t counter = 0) {
  Assert(end <= end.getVector().end());
  using T = typename std::iterator_traits<OutputIterator>::value_type;
  auto& vec = begin.getVector();
  using Vec = std::remove_reference_t<decltype(vec)>;
  auto& page = vec.getPageInstance();
  constexpr bool auth = Vec::auth;
  using Page = typename Vec::Page;

  static size_t item_per_page = Vec::item_per_page;
  size_t beginOffset = begin.get_page_offset();
  size_t beginPageIdx = begin.get_page_idx();
  size_t endOffset = end.get_page_offset();
  size_t endPageIdx = end.get_page_idx();

  if (beginOffset != 0) {
    size_t endPageOffset =
        beginPageIdx == endPageIdx ? endOffset : item_per_page;
    if constexpr (auth) {
      vec.server.Read(beginPageIdx, page, counter);
    } else {
      vec.server.Read(beginPageIdx, page);
    }
    for (size_t pageOffset = beginOffset; pageOffset != endPageOffset;
         ++pageOffset) {
      *to = page.pages[pageOffset];
      ++to;
    }
    if (beginPageIdx == endPageIdx) {
      // vec.server.flushRead();
      return to;
    }
    ++beginPageIdx;
  }
  for (size_t pageIdx = beginPageIdx; pageIdx < endPageIdx; ++pageIdx) {
    if constexpr (auth) {
      vec.server.ReadLazy(pageIdx, *(Page*)(&(*to)), counter);
    } else {
      vec.server.ReadLazy(pageIdx, *(Page*)(&(*to)));
    }

    to += item_per_page;
  }

  if (endOffset != 0) {
    if constexpr (auth) {
      vec.server.ReadLazy(endPageIdx, page, counter);
    } else {
      vec.server.ReadLazy(endPageIdx, page);
    }
    vec.server.flushRead();
    for (size_t pageOffset = 0; pageOffset != endOffset; ++pageOffset) {
      *to = page.pages[pageOffset];
      ++to;
    }
  } else {
    vec.server.flushRead();
  }

  return to;
}

}  // namespace EM::NonCachedVector