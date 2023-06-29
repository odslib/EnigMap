#pragma once
#include "sort.hpp"
// NOTE dummy flag is set at the least significant bit in this example
namespace EM::Algorithm {
using EM::NonCachedVector::Vector;
template <typename T, class Iterator>
void binPacking(Iterator begin, Iterator end, Iterator outputBegin, size_t Z,
                uint64_t bitMask) {
  size_t size = end - begin;
  size_t numBucket = size / Z;
  size_t logNumBucket = GetLogBaseTwo(numBucket);
  uint64_t baseMask = (bitMask >> logNumBucket) << 1;
  uint64_t rangeMask = baseMask * (numBucket - 1);
  uint64_t rangeMaskWithDummyFlag = rangeMask + 1UL;
  std::vector<TaggedT<T>> temp(size * 2);
  std::copy(begin, end, temp.begin());
  auto iter = temp.begin() + size;
  for (size_t i = 0; i < numBucket; ++i) {
    for (size_t j = 0; j < Z; ++j) {
      // dummyFlag at the end
      iter->tag = i * baseMask + 1UL;
      ++iter;
    }
  }
  const auto cmpTag = [rangeMaskWithDummyFlag = rangeMaskWithDummyFlag](
                          const auto& element1, const auto& element2) {
    return (element1.tag & rangeMaskWithDummyFlag) <
           (element2.tag & rangeMaskWithDummyFlag);
  };
  BitonicSort(temp, cmpTag);
  size_t count = 0;
  size_t curr = 0;
  for (auto& element : temp) {
    size_t newCurr = (element.tag & rangeMask);
    CMOV(newCurr != curr, count, 0UL);
    curr = newCurr;
    // excess real elements
    Assert(!(count >= Z && !(element.tag & 1UL)));
    // mark excessive element with -1
    CMOV(count >= Z, element.tag, -1UL);
    ++count;
  }
  BitonicSort(temp, cmpTag);  // could be replaced by OrCompact
  std::copy(temp.begin(), temp.begin() + size, outputBegin);
}

template <typename T, class Iterator>
void butterflyRouting(Iterator begin, Iterator end, Iterator outputBegin,
                      size_t Z, uint64_t bitMask, size_t gamma) {
  size_t size = end - begin;
  size_t numBucket = size / Z;
  if (numBucket <= gamma) {
    binPacking<T>(begin, end, outputBegin, Z, bitMask);
    return;
  }
  size_t num_layer = GetLogBaseTwo(numBucket);
  size_t half_num_layer = num_layer / 2;
  size_t numBucketPerBatch = 1UL << half_num_layer;
  size_t batchSize = Z * numBucketPerBatch;
  size_t batchCount = size / batchSize;
  std::vector<TaggedT<T>> nextLayer(size);
  for (size_t batchIdx = 0; batchIdx < batchCount; ++batchIdx) {
    auto batchBegin = begin + batchIdx * batchSize;
    auto batchEnd = begin + (batchIdx + 1) * batchSize;
    butterflyRouting<T>(batchBegin, batchEnd, batchBegin, Z, bitMask, gamma);
    for (size_t bucketIdx = 0; bucketIdx < numBucketPerBatch; ++bucketIdx) {
      std::copy(batchBegin + bucketIdx * Z, batchBegin + (bucketIdx + 1) * Z,
                nextLayer.begin() + bucketIdx * size / numBucketPerBatch +
                    batchIdx * Z);
    }
  }
  size_t nextLayerBucketPerBatch = 1UL << (num_layer - half_num_layer);
  size_t nextLayerBatchSize = nextLayerBucketPerBatch * Z;
  size_t nextLayerBatchCount = size / nextLayerBatchSize;
  for (size_t batchIdx = 0; batchIdx < nextLayerBatchCount; ++batchIdx) {
    butterflyRouting<T>(nextLayer.begin() + batchIdx * nextLayerBatchSize,
                        nextLayer.begin() + (batchIdx + 1) * nextLayerBatchSize,
                        outputBegin + batchIdx * nextLayerBatchSize, Z,
                        bitMask >> half_num_layer, gamma);
  }
}

template <typename T>
Vector<TaggedT<T>> tagAndPad(typename Vector<T>::Iterator begin,
                             typename Vector<T>::Iterator end, uint64_t Z) {
  size_t inputSize = end - begin;
  uint64_t intermidiateSize = GetNextPowerOfTwo(2 * (end - begin));
  Vector<TaggedT<T>> tv(intermidiateSize);
  uint64_t numBucket = intermidiateSize / Z;
  uint64_t numRealPerBucket = divRoundUp(inputSize, numBucket);
  typename Vector<T>::PrefetchReader inputReader(begin, end);
  typename Vector<TaggedT<T>>::Writer taggedTWriter(tv.begin(), tv.end());

  for (uint64_t bucketIdx = 0; bucketIdx < numBucket; ++bucketIdx) {
    for (uint64_t offset = 0; offset < Z; ++offset) {
      TaggedT<T> tt;
      tt.tag =
          UniformRandom() & -2UL;  // UNDONE: change to a secure rand function
      if (offset < numRealPerBucket && !inputReader.eof()) {
        tt.v = inputReader.read();
      } else {
        tt.tag |= 1UL;
      }
      taggedTWriter.write(tt);
    }
  }
  taggedTWriter.flush();
  return tv;
}

template <typename T>
void externalButterflyRouting(typename Vector<TaggedT<T>>::Iterator begin,
                              typename Vector<TaggedT<T>>::Iterator end,
                              typename Vector<TaggedT<T>>::Iterator outputBegin,
                              size_t Z, uint64_t bitMask,
                              uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  size_t size = end - begin;
  const size_t gamma = GetLogBaseTwo(size);
  if ((size * 2 + gamma * Z) * sizeof(TaggedT<T>) <= heapSize) {
    std::vector<TaggedT<T>> mem(size);
    CopyIn(begin, end, mem.begin());
    butterflyRouting<T>(mem.begin(), mem.end(), mem.begin(), Z, bitMask, gamma);
    CopyOut(mem.begin(), mem.end(), outputBegin);
    return;
  }
  size_t num_layer = GetLogBaseTwo(size / Z);
  size_t half_num_layer = num_layer / 2;
  size_t numBucketPerBatch = 1UL << half_num_layer;
  size_t batchSize = Z * numBucketPerBatch;
  size_t batchCount = size / batchSize;
  Vector<TaggedT<T>> nextLayer(size);
  std::vector<TaggedT<T>> temp(Z);
  for (size_t batchIdx = 0; batchIdx < batchCount; ++batchIdx) {
    auto batchBegin = begin + batchIdx * batchSize;
    auto batchEnd = begin + (batchIdx + 1) * batchSize;
    externalButterflyRouting<T>(batchBegin, batchEnd, batchBegin, Z, bitMask,
                                heapSize);
    for (size_t bucketIdx = 0; bucketIdx < numBucketPerBatch; ++bucketIdx) {
      CopyIn(batchBegin + bucketIdx * Z, batchBegin + (bucketIdx + 1) * Z,
             temp.begin());
      CopyOut(temp.begin(), temp.end(),
              nextLayer.begin() + bucketIdx * size / numBucketPerBatch +
                  batchIdx * Z);
    }
  }
  size_t nextLayerBucketPerBatch = 1UL << (num_layer - half_num_layer);
  size_t nextLayerBatchSize = nextLayerBucketPerBatch * Z;
  size_t nextLayerBatchCount = size / nextLayerBatchSize;
  for (size_t batchIdx = 0; batchIdx < nextLayerBatchCount; ++batchIdx) {
    externalButterflyRouting<T>(
        nextLayer.begin() + batchIdx * nextLayerBatchSize,
        nextLayer.begin() + (batchIdx + 1) * nextLayerBatchSize,
        outputBegin + batchIdx * nextLayerBatchSize, Z,
        bitMask >> half_num_layer, heapSize);
  }
}

template <class Iterator>
void CABucketShuffle(Iterator begin, Iterator end,
                     uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  const uint64_t Z = 512;
  Vector<TaggedT<T>> tv = tagAndPad<T>(begin, end, Z);
  externalButterflyRouting<T>(tv.begin(), tv.end(), tv.begin(), Z, 1UL << 63,
                              heapSize);
  typename Vector<T>::Writer outputWriter(begin, end);
  std::vector<TaggedT<T>> temp(Z);
  uint64_t prev = 0;
  for (size_t bucketIdx = 0; bucketIdx < tv.size() / Z; ++bucketIdx) {
    CopyIn(tv.begin() + bucketIdx * Z, tv.begin() + (bucketIdx + 1) * Z,
           temp.begin());
    BitonicSort(
        temp.begin(), temp.end(),
        [](const TaggedT<T>& a, const TaggedT<T>& b) { return a.tag < b.tag; });
    for (const TaggedT<T>& element : temp) {
      if (!(element.tag & 1UL)) {
        outputWriter.write(element.v);
        Assert(element.tag >= prev && (prev = element.tag || true));
      }
    }
  }
  outputWriter.flush();
}

template <class Iterator>
void CABucketSort(Iterator begin, Iterator end,
                  uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  CABucketShuffle(begin, end, heapSize);
  using T = typename std::iterator_traits<Iterator>::value_type;
  using Reader = typename Vector<T>::LazyPrefetchReader;
  size_t size = end - begin;
  Vector<T> batchSorted(size);
  std::vector<Reader> mergeReaders;
  size_t batchSize = heapSize / sizeof(T);
  size_t batchCount = divRoundUp(size, batchSize);
  mergeReaders.reserve(batchCount);
  std::vector<T> mem(batchSize);
  for (size_t batchIdx = 0; batchIdx != batchCount; ++batchIdx) {
    size_t len = std::min(batchSize, size - batchIdx * batchSize);
    CopyIn(begin + batchIdx * batchSize, begin + batchIdx * batchSize + len,
           mem.begin());
    std::sort(mem.begin(), mem.begin() + len);
    CopyOut(mem.begin(), mem.begin() + len,
            batchSorted.begin() + batchIdx * batchSize);
    mergeReaders.emplace_back(
        batchSorted.begin() + batchIdx * batchSize,
        batchSorted.begin() + (batchIdx * batchSize + len));
  }
  typename Vector<T>::Writer outputWriter(begin, end, 1);
  auto cmpmerge = [&](const auto& a, const auto& b) {
    return *b.second < *a.second;
  };
  using Reader = typename Vector<T>::LazyPrefetchReader;
  std::vector<std::pair<Reader*, T*>> heap;
  heap.reserve(mergeReaders.size() + 1);
  for (auto& reader : mergeReaders) {
    reader.init();
    heap.push_back({&reader, &reader.get()});
  }
  std::make_heap(heap.begin(), heap.end(), cmpmerge);
  while (!heap.empty()) {
    Reader* top = heap[0].first;
    outputWriter.write(top->read());
    if (!top->eof()) {
      // add a top at the end, which will be swapped to the top by pop_heap
      heap.emplace_back(top, &top->get());
    }
    std::pop_heap(heap.begin(), heap.end(), cmpmerge);
    heap.resize(heap.size() - 1);
  }
  outputWriter.flush();
}

template <typename Vec>
void CABucketSort(Vec& vec, uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  CABucketSort(vec.begin(), vec.end(), heapSize);
}

template <typename Vec>
void CABucketShuffle(Vec& vec, uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  CABucketShuffle(vec.begin(), vec.end(), heapSize);
}
}  // namespace EM::Algorithm