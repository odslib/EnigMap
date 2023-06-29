#pragma once
#include <iostream>
#include <vector>

#include "bitonic.hpp"
#include "block_for_sort.hpp"
#include "edge_rec.hpp"
#include "external_memory/extemvector.hpp"
#include "external_memory/noncachedvector.hpp"
#include "or_compact_shuffle.hpp"
#include "param_select.hpp"
#include "static_sort.hpp"

namespace EM::Algorithm {
using EM::NonCachedVector::Vector;

enum SortMethod {
  CABUCKETSORT,
  BITONICSORT,
  ORSHUFFLE,
  CABUCKETSHUFFLE,
  BITONICSHUFFLE,
  DISTRIBUTIONOSORT,
  KWAYBUTTERFLYOSORT,
  KWAYBUTTERFLYOSHUFFLE
};

enum PartitionMethod {
  INTERLEAVE_PARTITION,
  OR_COMPACT,
  GOODRICH_COMPACT,
  BITONIC
};

template <class Iterator>
void BucketObliviousSort(Iterator begin, Iterator end,
                         uint64_t heapSize = DEFAULT_HEAP_SIZE);

template <class Iterator>
void KWayButterflySort(Iterator begin, Iterator end, uint64_t heapSize);

template <class Iterator, class Check>
void InterleaveTwoWay(Iterator begin, Iterator end, const Check& isMarked) {
  size_t size = end - begin;
  Assert(size % 2 == 0);
  if (size == 2) {
    bool swapFlag = !isMarked(*begin);
    condSwap(swapFlag, *begin, *(begin + 1));
    return;
  }
  size_t leftHalfSize = (size + 2) / 4 * 2;
  bool targetBit;
  Iterator mid = begin + leftHalfSize;
  Iterator left = begin, right = mid;
  targetBit = !isMarked(*left);
  switch (size % 4) {
    case 0:
      ++left;
      ++right;
      break;
    case 2:
      left += 2;
      break;
  }
  // every time left tag bit != right tag bit, we swap so that the left tag bit
  // becomes target bit, then we negate the target bit. Eventually the number of
  // 0 and 1 will be the same on both half.
  for (; right != end; ++left, ++right) {
    bool leftTag = isMarked(*left);
    bool rightTag = isMarked(*right);
    bool swapFlag = leftTag != targetBit;
    targetBit = rightTag != swapFlag;
    condSwap(swapFlag, *left, *right);
  }
  InterleaveTwoWay(begin, mid, isMarked);
  InterleaveTwoWay(mid, end, isMarked);
}

template <class Iterator, class Check>
void InterleaveTwoWay(Iterator beginLeft, Iterator beginRight, size_t Z,
                      const Check& isMarked) {
  Assert(Z % 2 == 0);
  bool targetBit;
  Iterator endLeft = beginLeft + Z;
  Iterator endRight = beginRight + Z;
  Iterator left = beginLeft;
  Iterator right = beginRight;

  targetBit = isMarked(*right);
  ++left;
  ++right;

  // every time left tag bit != right tag bit, we swap so that the left tag bit
  // becomes target bit, then we negate the target bit. Eventually the number of
  // 0 and 1 will be the same on both half.
  for (; right != endRight; ++left, ++right) {
    bool leftTag = isMarked(*left);
    bool rightTag = isMarked(*right);
    bool swapFlag = leftTag != targetBit;
    targetBit = rightTag != swapFlag;
    condSwap(swapFlag, *left, *right);
  }

  InterleaveTwoWay(beginLeft, endLeft, isMarked);
  InterleaveTwoWay(beginRight, endRight, isMarked);
}

template <class Iterator, class MarkIterator>
void Permute(Iterator begin, Iterator end, MarkIterator marksBegin,
             MarkIterator marksEnd, int level = 0) {
  uint64_t size = marksEnd - marksBegin;
  Assert(size <= 16);
  Assert(marksEnd - marksBegin == end - begin);
  if (size == 1) {
    return;
  } else if (size == 2) {
    bool swapFlag = *marksBegin & 0xFU;
    condSwap(swapFlag, *begin, *(begin + 1));
    condSwap(swapFlag, *marksBegin, *(marksBegin + 1));
    return;
  }
  size_t halfSize = (size + 1) / 2;
  EdgeRec<uint64_t> rec(halfSize);
  Iterator mid = begin + halfSize;
  const MarkIterator marksMid = marksBegin + halfSize;
  MarkIterator marksLeft = marksBegin, marksRight = marksMid;
  for (; marksLeft != marksMid - 1; ++marksLeft, ++marksRight) {
    uint8_t leftLargeMask = -(uint8_t)((*marksLeft & 0xF) >= halfSize);
    *marksLeft -= leftLargeMask & halfSize;

    uint8_t rightLargeMask = -(uint8_t)((*marksRight & 0xF) >= halfSize);
    *marksRight -= rightLargeMask & halfSize;
    rec.flipEdge(*marksLeft & 0xF, *marksRight & 0XF);
    *marksLeft |= (0X10U << level) & leftLargeMask;
    *marksRight |= (0X10U << level) & rightLargeMask;
  }

  uint8_t leftLargeMask = -(uint8_t)((*marksLeft & 0xF) >= halfSize);
  *marksLeft -= leftLargeMask & halfSize;
  *marksLeft |= (0X10U << level) & leftLargeMask;

  uint8_t lastLeft = *marksLeft & 0xFU;
  uint8_t lastRight;
  if (size & 1) {
    lastRight = size - halfSize;
  } else {
    uint8_t rightLargeMask = -(uint8_t)((*marksRight & 0xF) >= halfSize);
    *marksRight -= rightLargeMask & halfSize;
    lastRight = *marksRight & 0xF;
    *marksRight |= (0X10U << level) & rightLargeMask;
  }
  rec.flipEdge(lastLeft, lastRight);
  EdgeRec<uint64_t> path = rec.EulerPath(halfSize);
  // if the last edge is in the wrong direction, reverse the path orientation
  // rather than do the swap
  uint64_t reversePathMask = path.retrieveAndFlipEdge(lastLeft, lastRight);
  path.flipEdge(-reversePathMask);
  Iterator left = begin, right = mid;
  for (MarkIterator marksLeft = marksBegin, marksRight = marksMid;
       marksLeft != marksMid - 1; ++marksLeft, ++marksRight, ++left, ++right) {
    uint8_t leftMark = *marksLeft & 0xFU;
    uint8_t rightMark = *marksRight & 0xFU;
    bool swapFlag = path.retrieveAndFlipEdge(leftMark, rightMark);
    condSwap<false>(swapFlag, *marksLeft, *marksRight);
    condSwap(swapFlag, *left, *right);
  }
  Permute(begin, mid, marksBegin, marksMid, level + 1);
  Permute(mid, end, marksMid, marksEnd, level + 1);
  left = begin, right = mid;
  for (MarkIterator marksLeft = marksBegin, marksRight = marksMid;
       marksRight != marksEnd; ++marksLeft, ++marksRight, ++left, ++right) {
    bool swapFlag = !(*marksRight & (0X10U << level));

    condSwap(swapFlag, *left, *right);
    condSwap(swapFlag, *marksLeft, *marksRight);
  }
}

template <class Iterator, class MarkIterator>
INLINE void InterleaveBaseCase(Iterator begin, Iterator end,
                               MarkIterator marksBegin, MarkIterator marksEnd,
                               const uint16_t k) {
  // TODO find a cleaner way to write this
  static constexpr StaticSort<3> staticSort3;
  static constexpr StaticSort<4> staticSort4;
  static constexpr StaticSort<5> staticSort5;
  static constexpr StaticSort<6> staticSort6;
  static constexpr StaticSort<7> staticSort7;
  static constexpr StaticSort<8> staticSort8;
  auto iter_pair = std::make_pair(begin, marksBegin);
  // should be branchless. k doesn't need to be oblivious
  switch (k) {
    case 3:
      staticSort3(iter_pair);
      break;
    case 4:
      staticSort4(iter_pair);
      break;
    case 5:
      staticSort5(iter_pair);
      break;
    case 6:
      staticSort6(iter_pair);
      break;
    case 7:
      staticSort7(iter_pair);
      break;
    case 8:
      staticSort8(iter_pair);
      break;
    default:
      X_LOG("wrong way partition\n");
      abort();
  }
}

template <class Iterator, class MarkIterator>
void Interleave(Iterator begin, Iterator end, MarkIterator marksBegin,
                MarkIterator marksEnd, const uint16_t k) {
  uint64_t size = marksEnd - marksBegin;
  Assert(size % k == 0);
  Assert(size / k == GetNextPowerOfTwo(size / k));
  Assert(marksEnd - marksBegin == end - begin);
  if (size == k) {
    // 23% of runtime
    InterleaveBaseCase(begin, end, marksBegin, marksEnd, k);
    return;
  }
  EdgeRec<uint64_t> rec(k);
  Iterator mid = begin + size / 2;
  MarkIterator marksMid = marksBegin + size / 2;

  // 3% of runtime
  for (MarkIterator marksLeft = marksBegin, marksRight = marksMid;
       marksRight != marksEnd; ++marksLeft, ++marksRight) {
    rec.flipEdge(*marksLeft, *marksRight);
  }

  // 13% of runtime
  EdgeRec<uint64_t> path = rec.EulerPath(size / 2);

  // if the first edge is in the wrong direction, reverse the path orientation
  // rather than do the swap
  uint64_t reversePathMask = path.retrieveAndFlipEdge(*marksBegin, *marksMid);
  path.flipEdge(-reversePathMask);

  Iterator left = begin + 1, right = mid + 1;

  for (MarkIterator marksLeft = marksBegin + 1, marksRight = marksMid + 1;
       marksRight != marksEnd; ++marksLeft, ++marksRight, ++left, ++right) {
    bool swapFlag = path.retrieveAndFlipEdge(*marksLeft, *marksRight);

    // 5%
    condSwap<false>(swapFlag, *marksLeft, *marksRight);

    // 56%
    condSwap(swapFlag, *left, *right);
  }
  Interleave(begin, mid, marksBegin, marksMid, k);
  Interleave(mid, end, marksMid, marksEnd, k);
}

template <class Iterator, typename Indicator>
void MergeSplitTwoWay(Iterator beginLeft, Iterator beginRight, size_t Z,
                      Indicator indicator,
                      const PartitionMethod method = OR_COMPACT);

template <typename Iterator>
void MergeSplitKWay(std::vector<Iterator>& begins, const size_t Z) {
  const size_t k = begins.size();
  Assert(k >= 2);
  Assert(k <= 8);
  if (k == 2) {
    MergeSplitTwoWay(begins[0], begins[1], Z, 1UL, INTERLEAVE_PARTITION);
    for (auto it = begins[0]; it != begins[0] + Z; ++it) {
      it->tag = (int64_t)it->tag >> 1;
      // arithmetic shift to keep the dummy flag
    }
    for (auto it = begins[1]; it != begins[1] + Z; ++it) {
      it->tag = (int64_t)it->tag >> 1;
    }
    return;
  }
  __m256i remainCounts = _mm256_set1_epi32(Z);  // could change to a larger reg
  using T = typename std::iterator_traits<Iterator>::value_type;
  T* temp = (T*)malloc(k * Z * sizeof(T));
  uint8_t* marks = (uint8_t*)malloc(k * Z);
  if (!temp) {
    printf("run out of memory when malloc temporary array for mergesplit");
    abort();
  }
  if (!marks) {
    printf("run out of memory when malloc temporary array for marks");
    abort();
  }
  auto* tempIt = temp;
  auto* marksIt = marks;
  for (auto begin : begins) {
    for (auto it = begin; it != begin + Z; ++it, ++marksIt, ++tempIt) {
      std::memcpy(tempIt, &(*it), sizeof(T));
      uint8_t isDummy = tempIt->isDummy();
      uint8_t mark =
          tempIt->getMarkAndUpdate(k) | (-isDummy);  // set to -1 if it's dummy

      Assert(isDummy <= 1);
      *marksIt = mark;
      remainCounts = mm256_decrement_epi32_var_indx(remainCounts, mark);
    }
  }

  if constexpr (IO_ROUND > 0) {
    if (_mm256_movemask_epi8(
            _mm256_cmpgt_epi32(_mm256_set1_epi32(1), remainCounts))) {
      printf("Not enough ones\n");
      std::abort();  // failed due to the lack of dummy elements
    }
  }
  int currMark = 0;
  uint32_t currRemain = mm256_extract_epi32_var_indx(remainCounts, 0);
  // assume that there's at least one dummy for each mark
  for (marksIt = marks; marksIt != marks + Z * k; ++marksIt) {
    uint8_t isDummy = (*marksIt == (uint8_t)-1);
    CMOV(isDummy, *marksIt, (uint8_t)currMark);
    currRemain -= isDummy;
    currMark += !currRemain;
    uint32_t remainCount = mm256_extract_epi32_var_indx(remainCounts, currMark);
    CMOV(!currRemain, currRemain, remainCount);
    Assert(currRemain > 0);
  }
  Interleave(temp, temp + k * Z, marks, marks + k * Z, k);
  free(marks);
  for (size_t i = 0; i < Z; ++i) {
    for (size_t j = 0; j < k; ++j) {
      memcpy(&*(begins[j] + i), temp + (i * k + j), sizeof(T));
    }
  }
  free(temp);
}

template <class Iterator, class Check>
void InterleaveTwoWayPartition(Iterator begin, Iterator end,
                               const Check& isMarked) {
  uint64_t size = (end - begin);
  Assert(size % 2 == 0);
  uint64_t Z = size / 2;
  InterleaveTwoWay(begin, end, isMarked);
  for (auto it1 = begin + 1, it2 = begin + Z + Z % 2; it2 < end;
       it1 += 2, it2 += 2) {
    swap(*it1, *it2);
  }
}

template <class Iterator, class Check>
void InterleaveTwoWayPartition(Iterator beginLeft, Iterator beginRight,
                               size_t Z, const Check& isMarked) {
  InterleaveTwoWay(beginLeft, beginRight, Z, isMarked);
  auto endRight = beginRight + Z;
  for (auto it1 = beginLeft + 1, it2 = beginRight + Z % 2; it2 < endRight;
       it1 += 2, it2 += 2) {
    swap(*it1, *it2);
  }
}

/**
 * Partition array such that all elements with tag & bitMask == 0 are at the
 * left half, and all elements with tag & bitMask == 1 are at the right half
 * (i.e., 0000...1111...). Note that the left and right half of the subarrays
 * doesn't need to be consecutive. This allows us to merge split in place.
 * @param[in] bitMask the mask of a single bit to partition
 * @param[in] beginLeft the beginning iterator (inclusive) of the left half to
 * partition
 * @param[in] beginRight the beginning iterator (inclusive) of the right half to
 * partition
 * @tparam Z the bucketsize
 * @pre size of the (sub)array must be a power of two
 *
 * */
template <class Iterator, typename Indicator>
void MergeSplitInPlace(Iterator begin, Iterator end, Indicator indicator,
                       const PartitionMethod method = INTERLEAVE_PARTITION) {
  // Balance the number of 0 and 1 at the indicator bit
  uint64_t markCount = 0;
  uint64_t n = end - begin;
  uint64_t Z = n / 2;
  for (auto it = begin; it != end; ++it) {
    markCount += it->setAndGetMarked(indicator);
    // for struct with separate flag, this
    // should also update the flag
  }
  const bool dir = markCount > Z;
  uint64_t diff = Z - markCount;
  CMOV(dir, diff, -diff);

  for (auto it = begin; it != end; ++it) {
    bool isMarked = it->isMarked(indicator);
    bool isDummy = it->isDummy();
    bool changeMark = isDummy & (!!diff) & (isMarked == dir);
    it->condChangeMark(changeMark, indicator);
    diff -= (uint64_t)changeMark;
  }
  if constexpr (IO_ROUND > 0) {
    if (diff) {
      printf("Not enough ones\n");
      std::abort();  // failed due to the lack of dummy elements
    }
  }

  const auto isMarked = [indicator = indicator](const auto& element) {
    return element.isMarked(indicator);
  };
  // Partition 0 and 1 to even and odd bits
  switch (method) {
    case INTERLEAVE_PARTITION:
      InterleaveTwoWayPartition(begin, end, isMarked);
      break;
    case OR_COMPACT:
      OrCompact(begin, end, isMarked);
      break;
    case GOODRICH_COMPACT:
      GoodrichCompact(begin, end, isMarked);
      break;
    case BITONIC:
      auto cmp = [=](const auto& element1, const auto& element2) {
        return (element1.isMarked(indicator) > (element2.isMarked(indicator))) |
               ((element2.isDummy() & !(element2.isMarked(indicator))));
      };
      BitonicSort(begin, end, cmp);
  }
}

template <class Iterator, typename Indicator>
void MergeSplitTwoWay(Iterator beginLeft, Iterator beginRight, size_t Z,
                      Indicator indicator, const PartitionMethod method) {
  // Balance the number of 0 and 1 at the indicator bit
  if (method != INTERLEAVE_PARTITION) {
    using T = typename std::iterator_traits<Iterator>::value_type;
    std::vector<T> temp(2 * Z);
    std::copy(beginLeft, beginLeft + Z, temp.begin());
    std::copy(beginRight, beginRight + Z, temp.begin() + Z);
    MergeSplitInPlace(temp.begin(), temp.end(), indicator, method);
    std::copy(temp.begin(), temp.begin() + Z, beginLeft);
    std::copy(temp.begin() + Z, temp.end(), beginRight);
    return;
  }
  uint64_t markCount = 0;
  Iterator endLeft = beginLeft + Z, endRight = beginRight + Z;
  Iterator end = endLeft;
  for (auto it = beginLeft;; ++it) {
    if (it == end) {
      if (it == endRight) {
        break;
      }
      it = beginRight;
      end = endRight;
    }
    markCount +=
        it->setAndGetMarked(indicator);  // for struct with separate flag, this
                                         // should also update the flag
  }
  const bool dir = markCount > Z;
  uint64_t diff = Z - markCount;
  CMOV(dir, diff, -diff);

  end = endLeft;
  for (auto it = beginLeft;; ++it) {
    if (it == end) {
      if (it == endRight) {
        break;
      }
      it = beginRight;
      end = endRight;
    }
    bool isMarked = it->isMarked(indicator);
    bool isDummy = it->isDummy();
    bool changeMark = isDummy & (!!diff) & (isMarked == dir);
    it->condChangeMark(changeMark, indicator);
    diff -= (uint64_t)changeMark;
  }
  if constexpr (IO_ROUND > 0) {
    if (diff) {
      printf("First level not enough ones\n");
      std::abort();  // failed due to the lack of dummy elements
    }
  }

  const auto isMarked = [indicator = indicator](const auto& element) {
    return element.isMarked(indicator);
  };
  // Partition 0 and 1 to even and odd bits

  InterleaveTwoWayPartition(beginLeft, beginRight, Z, isMarked);
}

/**
 * Partition dummy elements to the back
 * @param[in] begin the begin iterator of the range
 * @param[in] end the end iterator of the range
 * @return the end iterator of elements that are not dummy
 */
template <typename Iterator>
Iterator partitionDummy(Iterator begin, Iterator end) {
  for (auto left = begin, right = end - 1;;) {
    while (!left->isDummy()) {
      ++left;
      if (right <= left) {
        return left;
      }
    }

    while (right->isDummy()) {
      --right;
      if (right <= left) {
        return left;
      }
    }
    swap(*left, *right);
  }
}

/**
  begin and end defines the input array to sample from,
  alpha is the sampling ratio,
  M is the number of elements to sample into enclave each time
  returns the array of samples
  */
template <typename T, const bool reservoir = false>
Vector<Block<T>> sampleForTight(typename Vector<T>::Iterator begin,
                                typename Vector<T>::Iterator end, size_t M,
                                double alpha, double slack_sampling) {
  Assert(alpha > 0);
  Assert(slack_sampling * alpha < 1);
  Assert(slack_sampling > 1);
  Assert(begin < end);
  size_t N = end - begin;
  size_t expectedSampleSize = alpha * (double)N;
  size_t sampleSize = 0;
  if constexpr (reservoir) {
    sampleSize = expectedSampleSize;
  }
  size_t numBatch = N / M + 1;
  Vector<Block<T>> Gamma((size_t)(slack_sampling * M * alpha) * numBatch);
  auto outIt = Gamma.begin();
  {
    std::vector<Block<T>> Mem(M);  // add 64 bytes to make sure malloc work
    size_t Np = N;
    size_t np = expectedSampleSize;

    typename Vector<T>::PrefetchReader inputReader(begin, end);
    auto isMarked = [](const Block<T>& element) { return !element.isDummy(); };

    for (size_t i = 0; i < numBatch; ++i) {
      size_t batchSize = std::min(M, N - i * M);
      for (size_t j = 0; j < batchSize; ++j) {
        if constexpr (reservoir) {
          uint64_t z = UniformRandom(Np - 1);
          bool chosen = z < np;
          CMOV(!chosen, Mem[j], DUMMY<Block<T>>());
          auto realData = Block<T>();
          realData.setData(inputReader.read());
          CMOV(chosen, Mem[j], realData);
          np -= chosen;
          --Np;
        } else {
          uint64_t z = UniformRandom(1, N);
          bool chosen = z <= expectedSampleSize;
          auto realData = Block<T>();
          realData.setData(inputReader.read());
          CMOV(chosen, Mem[j], realData);
          CMOV(!chosen, Mem[j], DUMMY<Block<T>>());
          CMOV(chosen, sampleSize, sampleSize + 1);
        }
      }
      OrCompact(Mem.begin(), Mem.begin() + batchSize, isMarked);
      size_t len = std::min(batchSize, (size_t)(slack_sampling * alpha * M));
      CopyOut(Mem.begin(), Mem.begin() + len, outIt);
      outIt += len;
    }
    // next, sort the samples
    size_t prefixLen = outIt - Gamma.begin();
    if (prefixLen < M) {
      // can be solved in memory
      CopyIn(Gamma.begin(), outIt, Mem.begin());
      BitonicSort(Mem.begin(), Mem.begin() + prefixLen);
      return Vector<Block<T>>(Mem.begin(), Mem.begin() + sampleSize);
    }
  }
  KWayButterflySort(Gamma.begin(), outIt, M * sizeof(Block<T>));
  /** hold some heap place to avoid fragmentation */
  void* heapPlaceholder = malloc(32 + (M + 1) * sizeof(Block<T>));
  Assert(heapPlaceholder);

  Vector<Block<T>> samples(sampleSize);
  typename Vector<Block<T>>::PrefetchReader reader(
      Gamma.begin(), Gamma.begin() + sampleSize, 1);
  typename Vector<Block<T>>::Writer writer(samples.begin(), samples.end());
  uint64_t prevTag = 0;
  while (!reader.eof()) {
    const auto& data = reader.read();
    Assert(!data.isDummy());
    Assert(data.data.key >= prevTag);
    prevTag = data.data.key;
    writer.write(data);
  }
  writer.flush();
  free(heapPlaceholder);

  return samples;
}

template <typename T>
std::vector<T> getPQuantile(typename Vector<T>::Iterator begin,
                            typename Vector<T>::Iterator end, size_t p) {
  size_t n = end - begin;
  std::vector<T> quantiles(p - 1, DUMMY<T>());
  for (size_t k = 1; k < p; ++k) {
    quantiles[k - 1] = *(begin + (n * k) / p);
  }
  return quantiles;
}

/**
 * Used to sort tags
 */
template <typename T, SortMethod task, typename WrappedT>
class ButterflySorter {
 private:
  using Iterator = typename Vector<T>::Iterator;

  uint64_t Z;
  uint64_t numTotalBucket;
  uint64_t numRealPerBucket;

  // the number of buckets that fit in heapSize memory
  uint64_t heapSize;
  uint64_t numBucketFit;
  uint64_t numElementFit;

  DistriOSortParams distributionOParams = {};
  KWayButterflyParams KWayParams = {};

  // all pivots for distribution osort
  std::vector<WrappedT> allPivots;

  Vector<T> mergeSortFirstLayer;
  std::vector<std::pair<Iterator, Iterator>> mergeSortRanges;

  typename Vector<T>::Writer mergeSortFirstLayerWriter;
  typename Vector<T>::PrefetchReader inputReader;
  typename Vector<T>::Writer outputWriter;

 public:
  ButterflySorter(typename Vector<T>::Iterator inputBeginIt,
                  typename Vector<T>::Iterator inputEndIt,
                  uint64_t _heapSize = DEFAULT_HEAP_SIZE)
      : inputReader(inputBeginIt, inputEndIt),
        mergeSortFirstLayer(
            task == KWAYBUTTERFLYOSORT ? inputEndIt - inputBeginIt : 0),
        heapSize(_heapSize),
        numElementFit(_heapSize / sizeof(WrappedT)) {
    size_t size = inputEndIt - inputBeginIt;
    if constexpr (task == KWAYBUTTERFLYOSORT || task == KWAYBUTTERFLYOSHUFFLE) {
      KWayParams = bestKWayButterflyParams(size, numElementFit);
      Z = KWayParams.Z;
      Assert(numElementFit > 8 * Z);
      numElementFit -= Z * 8;
      numBucketFit = numElementFit / Z;
      numTotalBucket = KWayParams.totalBucket;
      numRealPerBucket = 1 + (size - 1) / numTotalBucket;
    } else {
      static_assert(task == DISTRIBUTIONOSORT);
      distributionOParams = bestDistriOSortParams(size, numElementFit);
      Z = distributionOParams.Z;
      Assert(Z % 2 == 0);
      numBucketFit = 1UL << GetLogBaseTwo(numElementFit / Z);
      numTotalBucket = 2UL << GetLogBaseTwo(size / Z);
      size_t p = numTotalBucket / numBucketFit;
      Vector<WrappedT> samples = sampleForTight<T>(
          inputBeginIt, inputEndIt,
          numElementFit / (sizeof(WrappedT) + 8) * sizeof(WrappedT),
          distributionOParams.alpha,
          distributionOParams
              .slack_sampling);  // leave some space for prefix sums
      allPivots = getPQuantile<WrappedT>(samples.begin(), samples.end(), p);

      numRealPerBucket = 1 + (size - 1) / numTotalBucket;
    }

    if constexpr (task == KWAYBUTTERFLYOSORT) {
      mergeSortFirstLayerWriter.init(mergeSortFirstLayer.begin(),
                                     mergeSortFirstLayer.end());
    } else {
      outputWriter.init(inputBeginIt, inputEndIt, 1);
    }
  }

  const std::vector<std::pair<Iterator, Iterator>>& getMergeSortBatchReaders() {
    return mergeSortRanges;
  }

  template <class Iterator>
  void KWayButterflySortBasic(Iterator begin, Iterator end, size_t ioLayer,
                              size_t innerLayer) {
    uint64_t numElement = end - begin;
    uint64_t numBucket = numElement / Z;
    uint64_t way = KWayParams.ways[ioLayer][innerLayer];
    Assert(numElement % Z == 0);
    Assert(numBucket % way == 0);
    uint64_t waySize = numElement / way;
    uint64_t wayBucket = numBucket / way;
    if (innerLayer > 0) {
      for (uint64_t i = 0; i < way; ++i) {
        KWayButterflySortBasic(begin + i * waySize, begin + (i + 1) * waySize,
                               ioLayer, innerLayer - 1);
      }
    } else {
      Assert(numBucket == way);
      if (ioLayer == 0) {
        Assert(waySize == Z);
        // tag and pad input
        for (uint64_t i = 0; i < way; ++i) {
          auto it = begin + i * waySize;
          for (uint64_t offset = 0; offset < Z; ++offset, ++it) {
            if (offset < numRealPerBucket && !inputReader.eof()) {
              it->v = inputReader.read();
              it->setTag(UniformRandom());
            } else {
              it->setDummy();
            }
          }
        }
      }
    }
    std::vector<Iterator> KWayIts(way);
    for (uint64_t j = 0; j < wayBucket; ++j) {
      for (uint64_t i = 0; i < way; ++i) {
        KWayIts[i] = begin + (i * wayBucket + j) * Z;
      }
      MergeSplitKWay(KWayIts, Z);
    }
  }

  template <class Iterator, typename PivotType = typename std::iterator_traits<
                                Iterator>::value_type>
  void DistriOSortBasic(Iterator begin, Iterator end, uint64_t layer,
                        const std::vector<PivotType>& pivots) {
    uint64_t numElement = end - begin;
    uint64_t numBucket = numElement / Z;
    uint64_t pivotNum = pivots.size();
    auto mid = begin + numElement / 2;

    if (numBucket > 1) {
      DistriOSortBasic(begin, mid, layer - 1, pivots);
      DistriOSortBasic(mid, end, layer - 1, pivots);
    } else {
      if (layer == 0) {
        // tag and pad input
        auto it = begin;
        for (uint64_t offset = 0; offset < Z; ++offset, ++it) {
          if (offset < numRealPerBucket && !inputReader.eof()) {
            it->setData(inputReader.read());
          } else {
            it->setDummy();
          }
        }
      }
      return;
    }
    if (pivotNum + 1 < numBucket) {
      return;
    }
    auto outputIt = begin;  // used to compact non-dummies to the front

    uint64_t pivotIdxOffset = (pivotNum + 1) / numBucket - 1;
    uint64_t i = 0;
    for (auto beginLeft = begin, beginRight = mid; beginLeft != mid;
         beginLeft += Z, beginRight += Z, ++i) {
      uint64_t pivotIdx =
          pivotIdxOffset +
          (reverseLowest32Bits(i & (numBucket / 2 - 1)) * (pivotNum + 1) >> 32);

      MergeSplitTwoWay(beginLeft, beginRight, Z, pivots[pivotIdx]);
    }
  }

  struct Dim2D {
    size_t row;
    size_t col;
  };
  std::vector<Dim2D> transposeData;

  void logicalTranspose(size_t row, size_t col) {
    std::reverse(transposeData.begin(), transposeData.end());
    transposeData.push_back({row, col});
    std::reverse(transposeData.begin(), transposeData.end());
  }

  template <class Iterator>
  inline Iterator mapIterator(Iterator it) const {
    size_t idx = it.get_m_ptr();
    size_t bucketIdx = idx / Z;

    const size_t offset = idx % Z;
    for (const Dim2D& dim : transposeData) {
      size_t transposeIdx = bucketIdx % (dim.row * dim.col);
      bucketIdx -= transposeIdx;
      size_t rowIdx = transposeIdx / dim.col;
      if constexpr (task == DISTRIBUTIONOSORT) {
        rowIdx = reverseLowest32Bits(rowIdx) * dim.row >> 32;
        // reverse the bits
      }
      const size_t colIdx = transposeIdx % dim.col;
      transposeIdx = colIdx * dim.row + rowIdx;
      bucketIdx += transposeIdx;
    }
    idx = bucketIdx * Z + offset;
    return Iterator(idx, it.getVector());
  }

  /**
   * Sort the tags of elments the vector between [begin, end)
   * When the problem can fit into heapSize, we use call sortBasic, otherwise we
   * either: Divide the problem into two subproblems, solve both, and run
   * merge-splits for the last layer; or Divide the problem to ~sqrt(bucketNum)
   * subproblems, run transpose after solving all of them, and solving
   * ~sqrt(bucketNum) subproblems again. Note that when bucketNum is not a
   * square, we solve sqrt(bucketNum * 2) subproblems, and group two consecutive
   * buckets together in transposition. The best strategy is calculated during
   * class initilization
   *
   * @param[in] layer the layer to start with in the butterfly network, which
   * decides the bit mask in MergeSplit
   * @param[in] begin the begin iterator of v to sort
   * @param[in] end the end iterator of v to sort
   * @param[in] stride the distance between the two buckets that mergesplit
   * receives. E.g., in layer 0, we always mergesplit adjacent buckets, so
   * stride = 1.
   * @post the tags within each bucket is not sorted, a bitonic sort should be
   * applied before we filter dummy elements.
   *
   */
  void DistriOSort(typename Vector<WrappedT>::Iterator begin,
                   typename Vector<WrappedT>::Iterator end, uint64_t layer) {
    uint64_t numElement = end - begin;
    uint64_t numBucket = numElement / Z;
    uint64_t totalLayer = GetLogBaseTwo(allPivots.size() + 1);
    X_LOG("total layer", totalLayer);
    const uint64_t numBucketPerBatch = std::min(numBucket, numBucketFit);
    const uint64_t logNumBucketPerBatch = GetLogBaseTwo(numBucketPerBatch);
    const uint64_t batchCount = numBucket / numBucketPerBatch;
    const uint64_t numElementPerBatch = numBucketPerBatch * Z;
    std::vector<WrappedT> batch(numElementPerBatch);
    for (;;) {
      size_t logWay = std::min(totalLayer - layer, logNumBucketPerBatch);
      for (uint64_t batchIdx = 0; batchIdx < batchCount; ++batchIdx) {
        auto batchBegin = begin + numElementPerBatch * batchIdx;
        auto batchEnd = batchBegin + numElementPerBatch;
        if (layer != 0) {
          auto inMemBatchIt = batch.begin();
          for (auto originalIt = batchBegin; originalIt != batchEnd;
               originalIt += Z, inMemBatchIt += Z) {
            auto actualIt = mapIterator(originalIt);
            Assert(actualIt.get_m_ptr() % Z == 0);
            CopyIn(actualIt, actualIt + Z, inMemBatchIt);
          }
        }
        {
          if (logWay == 0) {
            static const auto notDummyMark = [](const auto& element) {
              return !element.isDummy();
            };
            static const auto notDummyComp = [](const auto& element,
                                                const auto& unused) {
              return !element.isDummy();
            };
            OrCompact(batch.begin(), batch.end(), notDummyMark);
            auto realEnd =
                std::lower_bound(batch.begin(), batch.end(), 0, notDummyComp);
            if constexpr (IO_ROUND == 0) {  // mock
              realEnd = batch.begin() + uint64_t((batch.end() - batch.begin()) *
                                                 (numRealPerBucket - 1) / Z);
            }
            BitonicSort(batch.begin(), realEnd);
            for (auto it = batch.begin(); it != realEnd; ++it) {
              outputWriter.write(it->getData());
            }
            continue;
          }
          size_t pivotRangeSize = (allPivots.size() + 1) >> layer;
          uint64_t pivotRangeIdx = (batchIdx << layer) / batchCount;
          size_t pivotInterval = pivotRangeSize >> logWay;
          std::vector<WrappedT> pivots;
          pivots.reserve((1UL << logWay) - 1);
          for (size_t i = 1; i < (1UL << logWay); ++i) {
            uint64_t pivotIdx =
                pivotRangeIdx * pivotRangeSize + i * pivotInterval - 1;
            pivots.push_back(allPivots[pivotIdx]);
          }
          DistriOSortBasic(batch.begin(), batch.end(),
                           layer + logNumBucketPerBatch, pivots);
        }

        auto inMemBatchIt = batch.begin();
        for (auto originalIt = batchBegin; originalIt != batchEnd;
             originalIt += Z, inMemBatchIt += Z) {
          auto actualIt = mapIterator(originalIt);
          CopyOut(inMemBatchIt, inMemBatchIt + Z, actualIt);
        }
      }

      layer += logWay;

      if (logWay == 0) {
        dbg_printf("pivot osort done\n");
        break;
      }
      dbg_printf("logical transpose\n");
      logicalTranspose(1UL << logWay, numTotalBucket >> layer);
    }
    outputWriter.flush();
  }

  template <typename NumType>
  static NumType getVecProduct(const std::vector<NumType>& vec) {
    NumType product = 1;
    for (const NumType& x : vec) {
      product *= x;
    }
    return product;
  }

  void KWayButterflySort(typename Vector<WrappedT>::Iterator begin,
                         typename Vector<WrappedT>::Iterator end) {
    size_t size = end - begin;
    size_t numIoLayer = KWayParams.ways.size();
    auto mergeFirstLayerIt = mergeSortFirstLayer.begin();

    std::vector<WrappedT> batch(numElementFit);
    for (size_t ioLayer = 0; ioLayer < numIoLayer; ++ioLayer) {
      size_t numInternalWay = getVecProduct(KWayParams.ways[ioLayer]);
      size_t batchSize = numInternalWay * Z;
      if (task == KWAYBUTTERFLYOSORT && ioLayer == numIoLayer - 1) {
        batchSize = numElementFit / batchSize * batchSize;
        // maximize the chunksize at the last layer
      }
      size_t batchCount = divRoundUp(size, batchSize);

      for (uint64_t batchIdx = 0; batchIdx < batchCount; ++batchIdx) {
        auto extBatchBegin = begin + batchSize * batchIdx;
        batchSize = std::min(batchSize, end - extBatchBegin);
        auto extBatchEnd = extBatchBegin + batchSize;
        if (ioLayer != 0) {  // fetch from intermediate ext vector
          auto inMemBatchIt = batch.begin();
          for (auto originalIt = extBatchBegin; originalIt != extBatchEnd;
               originalIt += Z, inMemBatchIt += Z) {
            auto actualIt = mapIterator(originalIt);
            Assert(actualIt.get_m_ptr() % Z == 0);

            CopyIn(actualIt, actualIt + Z, inMemBatchIt, ioLayer - 1);
          }
        }

        for (auto groupBegin = batch.begin();
             groupBegin < batch.begin() + batchSize;
             groupBegin += numInternalWay * Z) {
          KWayButterflySortBasic(groupBegin, groupBegin + numInternalWay * Z,
                                 ioLayer, KWayParams.ways[ioLayer].size() - 1);
          if (ioLayer == numIoLayer - 1) {
            // last layer, combine with bitonic sort and output
            auto cmpTag = [](const auto& a, const auto& b) {
              return a.tag < b.tag;
            };
            for (size_t i = 0; i < numInternalWay; ++i) {
              auto it = groupBegin + i * Z;
              Assert(it + Z <= batch.end());
              BitonicSort(it, it + Z, cmpTag);
              // for shuffling, output directly
              if constexpr (task == KWAYBUTTERFLYOSHUFFLE) {
                for (auto fromIt = it; fromIt != it + Z; ++fromIt) {
                  if (!fromIt->isDummy()) {
                    outputWriter.write(fromIt->getData());
                  }
                }
              }
            }
          }
        }
        if (ioLayer == numIoLayer - 1) {
          if constexpr (task == KWAYBUTTERFLYOSORT) {
            // sort the batch and write to first layer of merge sort
            auto cmpVal = [](const auto& a, const auto& b) {
              return a.v < b.v;
            };
            Assert(batch.begin() + batchSize <= batch.end());
            auto realEnd =
                partitionDummy(batch.begin(), batch.begin() + batchSize);
            // partition dummies to the end
            Assert(realEnd <= batch.end());
            std::sort(batch.begin(), realEnd, cmpVal);
            auto mergeSortReaderBeginIt = mergeSortFirstLayerWriter.it;
            for (auto it = batch.begin(); it != realEnd; ++it) {
              mergeSortFirstLayerWriter.write(it->getData());
            }

            mergeSortRanges.emplace_back(mergeSortReaderBeginIt,
                                         mergeSortFirstLayerWriter.it);
          }

        } else {
          auto inMemBatchIt = batch.begin();
          for (auto originalIt = extBatchBegin; originalIt != extBatchEnd;
               originalIt += Z, inMemBatchIt += Z) {
            auto actualIt = mapIterator(originalIt);

            CopyOut(inMemBatchIt, inMemBatchIt + Z, actualIt, ioLayer);
          }
        }
      }
      if (ioLayer < numIoLayer - 1) {
        logicalTranspose(numInternalWay, numTotalBucket / numInternalWay);
      }
    }
    if constexpr (task == KWAYBUTTERFLYOSORT) {
      mergeSortFirstLayerWriter.flush();
    } else {
      outputWriter.flush();
    }
  }

  void sort(typename Vector<WrappedT>::Iterator begin,
            typename Vector<WrappedT>::Iterator end) {
    if constexpr (task == DISTRIBUTIONOSORT) {
      DistriOSort(begin, end, 0);
    } else {
      KWayButterflySort(begin, end);
    }
  }

  size_t getOutputSize() { return numTotalBucket * Z; }
};

/**
 * Obliviously permute elements of a vector uniform randomly
 * @param[inout] b the vector to permute
 */
template <class Iterator>
void KWayBucketObliviousShuffle(Iterator begin, Iterator end,
                                uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  using T = typename std::iterator_traits<Iterator>::value_type;

  ButterflySorter<T, KWAYBUTTERFLYOSHUFFLE, TaggedT<T>> sorter(begin, end,
                                                               heapSize);
  Vector<TaggedT<T>> v(sorter.getOutputSize());
  sorter.sort(v.begin(), v.end());
}

template <typename Iterator>
void ExtMergeSort(Iterator begin, Iterator end,
                  const std::vector<std::pair<Iterator, Iterator>>& mergeRanges,
                  uint32_t outputCounter = 0) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  using Reader = typename Vector<T>::LazyPrefetchReader;
  typename Vector<T>::Writer outputWriter(begin, end, outputCounter);
  // for merge sort
  auto cmpmerge = [&](const auto& a, const auto& b) {
    return *(b.second) < *(a.second);
  };
  std::vector<Reader> mergeReaders;
  mergeReaders.reserve(mergeRanges.size());
  for (const std::pair<Iterator, Iterator>& range : mergeRanges) {
    mergeReaders.emplace_back(range.first, range.second);
  }
  std::vector<std::pair<Reader*, T*>> heap;
  heap.reserve(mergeRanges.size() + 1);
  for (auto& reader : mergeReaders) {
    reader.init();
    heap.push_back({&reader, &reader.get()});
  }
  std::make_heap(heap.begin(), heap.end(), cmpmerge);
  while (!heap.empty()) {
    Reader* top = heap[0].first;
    outputWriter.write(top->read());
    if (!top->eof()) {
      heap.emplace_back(top, &top->get());
      // add a top at the end, which will be swapped to the top by pop_heap
    }
    std::pop_heap(heap.begin(), heap.end(), cmpmerge);
    heap.resize(heap.size() - 1);
  }
  outputWriter.flush();
}

template <class Iterator>
void KWayButterflySort(Iterator begin, Iterator end, uint64_t heapSize) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  const uint64_t N = end - begin;
  ButterflySorter<T, KWAYBUTTERFLYOSORT, TaggedT<T>> sorter(begin, end,
                                                            heapSize);
  Vector<TaggedT<T>> v(sorter.getOutputSize());
  sorter.sort(v.begin(), v.end());
  using Reader = typename Vector<T>::LazyPrefetchReader;
  const auto& mergeRanges = sorter.getMergeSortBatchReaders();
  ExtMergeSort(begin, end, mergeRanges, 1);
}

template <typename Vec>
void KWayButterflySort(Vec& vec, uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  KWayButterflySort(vec.begin(), vec.end(), heapSize);
}

template <typename Vec>
void KWayBucketObliviousShuffle(Vec& vec,
                                uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  KWayBucketObliviousShuffle(vec.begin(), vec.end(), heapSize);
}

template <class Iterator>
void DistriOSort(Iterator begin, Iterator end,
                 uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  size_t N = end - begin;
  if (N <= heapSize / sizeof(T)) {
    std::vector<T> Mem(N);
    CopyIn(begin, end, Mem.begin());
    BitonicSort(Mem);
    CopyOut(Mem.begin(), Mem.end(), begin);
    return;
  }
  ButterflySorter<T, DISTRIBUTIONOSORT, Block<T>> sorter(begin, end, heapSize);
  Vector<Block<T>> v(sorter.getOutputSize());
  sorter.sort(v.begin(), v.end());
}

template <typename Vec>
void DistriOSort(Vec& vec, uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  DistriOSort(vec.begin(), vec.end(), heapSize);
}

}  // namespace EM::Algorithm