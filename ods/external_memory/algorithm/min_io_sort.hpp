#pragma once
#include <math.h>

#include <iostream>
#include <numeric>
#include <queue>
#include <vector>

#include "block_for_sort.hpp"
#include "sort.hpp"

namespace EM::Algorithm {
using EM::NonCachedVector::Vector;

template <typename Iterator>
void minIOSort(Iterator begin, Iterator end,
               uint64_t heapSize = DEFAULT_HEAP_SIZE);

// get p+1 iterator for pivots evenly between begin and end
template <typename Iterator>
std::vector<Iterator> getPQuantileIterators(Iterator begin, Iterator end,
                                            size_t p) {
  size_t n = end - begin;
  std::vector<Iterator> quantileIts(p + 1, begin);
  for (size_t k = 0; k < p; ++k) {
    quantileIts[k] = begin + (n * k) / p;
  }
  quantileIts[p] = end;
  return quantileIts;
}

// partition elements in [begin, end) into p segments of equal size in place
// using pivots
template <class Iterator, class PivotIterator>
void pWayPartitionHelper(Iterator begin, Iterator end, PivotIterator pivotBegin,
                         PivotIterator pivotEnd) {
  size_t p = pivotEnd - pivotBegin + 1;
  size_t size = end - begin;
  if (p == 1) {
    return;
  }
  size_t pLeft = p / 2;
  size_t pRight = p - pLeft;
  size_t leftSize = size * pLeft / p;
  size_t rightSize = size - leftSize;
  auto pivotIt = pivotBegin + (pLeft - 1);
  auto pivotVal = *pivotIt;
  size_t countLess = 0, countDummy = 0;
  for (auto it = begin; it != end; ++it) {
    const auto& it_const = it;
    bool dummyFlag = it_const->isDummy();
    bool lessFlag = (!dummyFlag) & (*it_const < pivotVal);
    Assert(!(dummyFlag && lessFlag));
    CMOV(lessFlag, countLess, countLess + 1);  // UNDONE
    CMOV(dummyFlag, countDummy, countDummy + 1);
  }

  for (auto it = begin; it != end; ++it) {
    const auto& it_const = it;
    bool isDummy = it_const->isDummy();
    bool addLessFlag = (countLess < leftSize) & isDummy;  // UNDONE
    it->setLessFlag(addLessFlag | ((!isDummy) & (*it_const < pivotVal)));
    CMOV(addLessFlag, countLess, countLess + 1);
  }
  if constexpr (IO_ROUND > 0) {
    if (countLess != leftSize) {
      dbg_printf("UNABLE to pad dummy elements with pivot\n");
      abort();
    }
  }
  auto isMarked = [](const auto& val) { return val.isLess(); };
  // uint64_t Z = (end - begin) / p;
  // for (auto it = begin; it <= end - 2 * Z; it += 2 * Z) {
  //     if ((p & 1) == 0) {
  //         InterleaveTwoWayPartition(it, it + 2 * Z, isMarked);
  //     } else {
  //         OrCompact(it, it + 2 * Z, isMarked);
  //     }
  // }
  if ((p & 1) == 0) {
    InterleaveTwoWayPartition(begin, end, isMarked);
  } else {
    OrCompact(begin, end, isMarked);
  }

  auto midIt = begin + leftSize;

  pWayPartitionHelper(begin, midIt, pivotBegin, pivotIt);
  pWayPartitionHelper(midIt, end, pivotIt + 1, pivotEnd);
}

/* partition elements in [begin, end) into p segments of equal size in place
  using pivots in [pivotsBegin, pivotsEnd) and output to [outputBegin,
  outputEnd) Note that output size can be larger than input size. In that case,
  fillers will be padded
*/
template <typename T>
void pWayPartition(typename Vector<T>::Iterator begin,
                   typename Vector<T>::Iterator end,
                   typename std::vector<T>::iterator pivotBegin,
                   typename std::vector<T>::iterator pivotEnd,
                   typename std::vector<T>::iterator outputBegin,
                   typename std::vector<T>::iterator outputEnd) {
  Assert(outputEnd - outputBegin == end - begin);
  CopyIn(begin, end, outputBegin);
  pWayPartitionHelper(outputBegin, outputEnd, pivotBegin, pivotEnd);
}

/* partition elements in [begin, end) into p segments of equal size in place
   using pivots in [pivotsBegin, pivotsEnd) and output to [outputBegin,
   outputEnd) Note that output size can be larger than input size. In that case,
   fillers will be padded
*/
template <typename T, typename InputIterator>
void pWayPartitionFirstLayer(InputIterator begin, InputIterator end,
                             typename std::vector<T>::iterator pivotBegin,
                             typename std::vector<T>::iterator pivotEnd,
                             typename std::vector<T>::iterator outputBegin,
                             typename std::vector<T>::iterator outputEnd) {
  Assert(begin <= end);
  Assert(outputBegin < outputEnd);
  Assert(outputEnd - outputBegin >= end - begin);
  using OriginalT = typename std::iterator_traits<InputIterator>::value_type;
  size_t p = pivotEnd - pivotBegin + 1;
  size_t Z = (outputEnd - outputBegin) / p;
  Assert(p * Z == (outputEnd - outputBegin));
  size_t Zr = divRoundUp(end - begin, p);
  // static_assert(std::is_same<Block<OriginalT, false>, T>::value); // input
  typename Vector<OriginalT>::Reader inputReader(begin, end);
  for (size_t j = 0; j < p; ++j) {
    size_t numReal = 0;
    auto bucketBegin = outputBegin + j * Z;
    auto bucketEnd = outputBegin + (j + 1) * Z;
    for (; numReal < Zr && !inputReader.eof(); ++numReal) {
      auto block = Block<OriginalT>();
      block.setData(inputReader.read());
      *(bucketBegin + numReal) = block;
    }
    std::fill(bucketBegin + numReal, bucketEnd, DUMMY<T>());
  }

  pWayPartitionHelper(outputBegin, outputEnd, pivotBegin, pivotEnd);
}

/* partition elements in [begin, end) into p segments of equal size in place
   using pivots and output to [outBegin, outEnd). Note that output size can be
   larger than input size. In that case, fillers will be padded. The partition
   is done in batchCount number of batches. M defines the output size of each
   batch.
   */

template <typename T, typename InputIterator>
void partitionInBatch(InputIterator begin, InputIterator end,
                      typename Vector<T>::Iterator outBegin,
                      typename Vector<T>::Iterator outEnd,
                      std::vector<T>& pivots, size_t batchCount,
                      std::vector<T>& batchOutput) {
  size_t N = end - begin;
  size_t outN = outEnd - outBegin;

  Assert(N <= outN);
  size_t Mp = (N - 1) / batchCount + 1;
  size_t p = pivots.size() + 1;
  Assert(outN % batchCount == 0);
  size_t M = outN / batchCount;
  Assert(M % p == 0);
  Assert(Mp <= M);
  size_t perBatchPartitionSize = M / p;
  size_t perPartitionSize = outN / p;
  size_t perBatchPartitionInterval =
      perPartitionSize /
      batchCount;  // we want dummies to be splited evenly within each partition
  constexpr bool isFirstLayer =
      !std::is_same<typename Vector<T>::Iterator, InputIterator>::value;
  for (size_t i = 0; i < batchCount; ++i) {
    auto batchBeginIt = begin + i * Mp;
    auto batchEndIt = begin + (i + 1) * Mp;
    if (end < batchBeginIt) {
      batchBeginIt = end;
    }
    if (end < batchEndIt) {
      batchEndIt = end;
    }
    if constexpr (isFirstLayer) {
      pWayPartitionFirstLayer<T>(batchBeginIt, batchEndIt, pivots.begin(),
                                 pivots.end(), batchOutput.begin(),
                                 batchOutput.begin() + M);
    } else {
      pWayPartition<T>(batchBeginIt, batchEndIt, pivots.begin(), pivots.end(),
                       batchOutput.begin(), batchOutput.begin() + M);
    }
    for (size_t j = 0; j < p; ++j) {
      auto toBegin =
          outBegin + (j * perPartitionSize + i * perBatchPartitionInterval);
      auto fromBegin = batchOutput.begin() + j * perBatchPartitionSize;
      auto fromEnd = fromBegin + perBatchPartitionSize;
      Assert(fromEnd <= batchOutput.begin() + M);
      CopyOut(fromBegin, fromEnd, toBegin);
    }
  }
}

/**
recursive function to partition [begin, end) given samples and output to
[outputBegin, outputEnd) using p way partition and batch size M
*/
template <typename T>
void multiLevelPartitionHelper(typename Vector<T>::Iterator begin,
                               typename Vector<T>::Iterator end,
                               typename Vector<T>::Iterator samplesBegin,
                               typename Vector<T>::Iterator samplesEnd,
                               typename Vector<T>::Iterator outputBegin,
                               typename Vector<T>::Iterator outputEnd, size_t p,
                               size_t M, std::vector<T>& batchOutput) {
  size_t N = end - begin;
  if (N <= M) {
    Assert(outputEnd - outputBegin >= end - begin);
    std::vector<T> temp(end - begin);
    CopyIn(begin, end, temp.begin());
    CopyOut(temp.begin(), temp.end(), outputBegin);
    // CopyForLateInit(begin, end, outputBegin); // UNDONE copy directly between
    // noncachedvectors
    Fill(outputBegin + N, outputEnd, DUMMY<T>());
    return;
  }
  std::vector<T> pivots = getPQuantile<T>(samplesBegin, samplesEnd, p);
  size_t batchCount = (N - 1) / M + 1;
  size_t perPartitionSize = N / p;
  // if (perPartitionSize * p != N) {
  //     std::cerr << "partition not divisible !!!!!" << std::endl;
  //     std::cerr << "perPartitionSize = " << perPartitionSize << std::endl;
  //     std::cerr << "N = " << N << std::endl;
  // }

  auto levelOutBegin = outputBegin;
  if (perPartitionSize > M) {
    Vector<T> X(N);
    levelOutBegin = X.begin();
    partitionInBatch(begin, end, X.begin(), X.end(), pivots, batchCount,
                     batchOutput);
    auto quantileIts = getPQuantileIterators(samplesBegin, samplesEnd, p);
    for (size_t j = 0; j < p; ++j) {
      multiLevelPartitionHelper<T>(
          X.begin() + j * perPartitionSize,
          X.begin() + (j + 1) * perPartitionSize,  // input
          quantileIts[j], quantileIts[j + 1],      // pivots
          outputBegin + j * perPartitionSize,
          outputBegin + (j + 1) * perPartitionSize,  // output
          p, M, batchOutput);
    }
  } else {
    partitionInBatch(begin, end, outputBegin, outputEnd, pivots, batchCount,
                     batchOutput);
  }

  // #ifndef NDEBUG
  // for (size_t j = 0; j < p; ++j) {
  //     size_t dummyCount = 0;
  //     for (size_t k = 0; k < perPartitionSize; ++k) {
  //         size_t idx = j * perPartitionSize + k;
  //         auto val = *(levelOutBegin + idx);
  //         if (!val.isDummy()) {
  //             if (j != 0 && val < pivots[j-1]) {
  //                 std::cerr << "Incorrect partition; partition "
  //                 << j << " offset " << k << " val " << val
  //                 << " smaller than prev pivot " << pivots[j-1] << std::endl;
  //                 // return;
  //             }
  //             if (j != p - 1 && pivots[j] < val) {
  //                 std::cerr << "Incorrect partition; partition "
  //                 << j << " offset " << k << " val " << val
  //                 << " larger than next pivot " << pivots[j] << std::endl;
  //                 // return;
  //             }

  //         } else {
  //             ++dummyCount;
  //         }
  //     }
  //     // std::cout << "partition " << j << " contains " << dummyCount << "
  //     dummies out of " << perPartitionSize << " elements" << std::endl;
  // }
  // #endif
}

/**
partition [begin, end) given samples and output to [outputBegin, outputEnd)
using p way partition except on the first layer and batch size M
NOTE THE FUNCTION WILL CHANGE MEMSIZE M!!!
*/
template <typename T>
Vector<Block<T>> multiLevelPartition(typename Vector<T>::Iterator begin,
                                     typename Vector<T>::Iterator end,
                                     Vector<Block<T>>& Gamma, size_t& M,
                                     std::vector<Block<T>>& Mem,
                                     const OQSortParams oqsort_params) {
  size_t N = end - begin;
  if (M >= N) {
    Vector<Block<T>> res(N);
    typename Vector<T>::Reader inputReader(begin, end);
    typename Vector<Block<T>>::Writer outputWriter(res.begin(), res.end());
    for (size_t i = 0; i < N; ++i) {
      T data = inputReader.read();
      Block<T> blockData = Block<T>();
      blockData.setData(data);
      outputWriter.write(blockData);
    }
    outputWriter.flush();
    return res;
  }
  size_t N_ = N * (1 + oqsort_params.eps);
  auto params = oqsort_params;
  dbg_printf("N_ = %ld\n", N_);
  size_t batchCount, p0, p, Mmax = M;
  int r = 1;
  for (size_t N_M = divRoundUp(N_, Mmax);;
       ++N_M) {  // in case M exceed limit after rounding up

    r = std::max(1, (int)ceil(log(N_M) / log(oqsort_params.p_max)));
    p = (size_t)ceil(pow(N_M, 1.0 / r));
    p0 = ceil(N_M / pow(p, r - 1));
    dbg_printf("pivots generated p0 = %ld, p = %ld, r = %d\n", p0, p, r);
    batchCount = p0 * (size_t)pow(p, r - 1);
    M = divRoundUp(N_, batchCount);
    size_t lcmpp0 = r == 1 ? p0 : std::lcm(p, p0);
    M = divRoundUp(M, lcmpp0) * lcmpp0;  // want M to be multiple of p and p0
    if (M <= Mmax) {
      break;
    }
  }

  size_t finalOutputSize = M * batchCount;
  dbg_printf("final OutputSize is %ld\n", finalOutputSize);
  // fisherYatesShuffle(begin, end); // UNDONE disable enclave mode here

  // std::cout << "finalOutputSize: " << finalOutputSize << std::endl;
  Vector<Block<T>> X1(finalOutputSize, DUMMY<Block<T>>());
  std::vector<Block<T>> pivots =
      getPQuantile<Block<T>>(Gamma.begin(), Gamma.end(), p0);
  // for (const auto& pivot: pivots) {
  //     std::cout << pivot << " ";
  // }
  // std::cout << std::endl;
  partitionInBatch(begin, end, X1.begin(), X1.end(), pivots, batchCount, Mem);
  dbg_printf("first level partition done\n");

  if (r == 1) {
    return X1;
  }
  size_t perPartitionSize = finalOutputSize / p0;

  Vector<Block<T>> finalOutput(finalOutputSize);
  auto quantileIts = getPQuantileIterators(Gamma.begin(), Gamma.end(), p0);
  for (size_t j = 0; j < p0; ++j) {
    multiLevelPartitionHelper<Block<T>>(
        X1.begin() + j * perPartitionSize,
        X1.begin() + (j + 1) * perPartitionSize,  // input
        quantileIts[j], quantileIts[j + 1],       // pivots
        finalOutput.begin() + j * perPartitionSize,
        finalOutput.begin() + (j + 1) * perPartitionSize,  // output
        p, M, Mem);
  }

  return finalOutput;
}

template <typename Iterator>
void minIOSort(Iterator begin, Iterator end, uint64_t heapSize) {
  using T = typename std::iterator_traits<Iterator>::value_type;
  size_t N = end - begin;
  if (N <= heapSize / sizeof(T)) {
    std::vector<T> Mem(N);
    CopyIn(begin, end, Mem.begin());
    BitonicSort(Mem);
    CopyOut(Mem.begin(), Mem.end(), begin);
    return;
  }

  size_t M = heapSize / sizeof(Block<T>);
  dbg_printf("M = %ld\n", M);

  OQSortParams oqsort_params = bestOQSortParams(N, M);

  dbg_printf(
      "alpha=%f, eps=%f, M=%zu, slack_sampling=%f, layer = %ld, p_max=%ld\n",
      oqsort_params.alpha, oqsort_params.eps, oqsort_params.M,
      oqsort_params.slack_sampling, oqsort_params.layer, oqsort_params.p_max);
  M = oqsort_params.M;

  Vector<Block<T>> Gamma = sampleForTight<T>(begin, end, M, oqsort_params.alpha,
                                             oqsort_params.slack_sampling);
  // if (heapPlaceholder) {
  //     free(heapPlaceholder);
  // }
  /** don't allocate anything here on heap*/
  std::vector<Block<T>> Mem(M);
  dbg_printf("sampling done\n");

  auto res = multiLevelPartition<T>(begin, end, Gamma, M, Mem, oqsort_params);
  // M here can be smaller than the original M, so don't use Mem.end()
  // freeing Mem can potentially cause fragmentation
  dbg_printf("partition done\n");

  typename Vector<T>::Writer outputWriter(begin, end);
  static const auto notDummyMark = [](const auto& element) {
    return !element.isDummy();
  };
  static const auto notDummyComp = [](const auto& element, const auto& unused) {
    return !element.isDummy();
  };
  for (auto it = res.begin(); it < res.end();) {
    auto batchEnd = it + M;
    if (res.end() < batchEnd) {
      batchEnd = res.end();
      std::fill(Mem.begin(), Mem.begin() + M, DUMMY<Block<T>>());
    }

    CopyIn(it, batchEnd, Mem.begin());
    OrCompact(Mem.begin(), Mem.begin() + M, notDummyMark);
    auto realEnd =
        std::lower_bound(Mem.begin(), Mem.begin() + M, 0, notDummyComp);
    if constexpr (IO_ROUND == 0) {
      realEnd = Mem.begin() + uint64_t(M / (1 + oqsort_params.eps));
    }
    BitonicSort(Mem.begin(), realEnd);
    for (auto it = Mem.begin(); it != realEnd; ++it) {
      outputWriter.write(it->data);
    }
    it = batchEnd;
  }
  outputWriter.flush();
}

template <typename Vec>
void minIOSort(Vec& vec, uint64_t heapSize = DEFAULT_HEAP_SIZE) {
  minIOSort(vec.begin(), vec.end(), heapSize);
}
}  // namespace EM::Algorithm