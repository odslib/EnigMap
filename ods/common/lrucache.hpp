#pragma once
#include <list>
#include <unordered_map>

namespace _ORAM {
template<typename T, uint64_t CACHE_SIZE=ORAM_SERVER__CACHE_SIZE>
struct Cache {
  typedef uint64_t IndexType;
  uint64_t size;
  std::unordered_map<IndexType, T> data;
  std::unordered_map<IndexType, std::list<IndexType>::iterator> dataRefs;
  std::list<IndexType> LRU;

  Cache() {
    size = 0;
    data.reserve(CACHE_SIZE);
    dataRefs.reserve(CACHE_SIZE);
  }

  bool CheckContains(const IndexType& rootIdx) {
    return data.count(rootIdx) > 0;
  }

  bool IsFull() {
    return size == CACHE_SIZE;
  }

  T& GetNextToEvict(IndexType& indexToEvict) {
    Assert(size > 0);
    Assert(size == CACHE_SIZE);
    indexToEvict = LRU.back();
    Assert(data.count(indexToEvict) > 0);
    return data[indexToEvict];
  }

  void EvictLRU(const IndexType& evictedIndex) {
    Assert(size > 0);
    Assert(size == CACHE_SIZE);
    size--;
    Assert(evictedIndex == LRU.back());

    LRU.pop_back();
    Assert(data.count(evictedIndex) > 0);
    data.erase(evictedIndex);
    dataRefs.erase(evictedIndex);
  }

  T& Insert(const IndexType& newIndex, const T& val) {
    Assert(size < CACHE_SIZE);
    Assert(data.count(newIndex) == 0);
    size++;
    LRU.push_front(newIndex);

    data[newIndex] = val;
    dataRefs[newIndex] = (LRU.begin());

    return data[newIndex];
  }

  T& Access(const IndexType& accessedIndex) {
    Assert(size > 0);
    Assert(data.count(accessedIndex) > 0);
    
    LRU.erase(dataRefs[accessedIndex]);
    LRU.push_front(accessedIndex);
    dataRefs[accessedIndex] = LRU.begin();

    return data[accessedIndex];
  }
};
}