#pragma once

namespace EM::Algorithm {
/**
 * Wrapper to parse a Vector as a square matrix and perform transposition
 */
template <typename Iterator>
class MatTransposer {
 private:
  uint64_t unitSize;
  uint64_t sideLen;
  uint64_t cacheEfficientSideLen;
  Iterator begin;
  using T = typename std::iterator_traits<Iterator>::value_type;
  /**
   * get the element at mat[row][col][offset]
   */
  T& get(uint64_t row, uint64_t col, uint64_t offset) {
    return *(begin + ((row * sideLen + col) * unitSize + offset));
  };

  /**
   * swap mat[row][col] and mat[col][row]
   * unitSize elements of v are swapped
   */
  void swapUnit(uint64_t row1, uint64_t col1, uint64_t row2, uint64_t col2) {
    for (uint64_t offset = 0; offset != unitSize; ++offset) {
      swap(get(row1, col1, offset), get(row2, col2, offset));
    }
  }

  /**
   * transpose two regions of size N * N and swap the two regions
   */
  void transposeAndSwap(uint64_t rowBegin1, uint64_t colBegin1,
                        uint64_t rowBegin2, uint64_t colBegin2, uint64_t N) {
    if (N <= cacheEfficientSideLen) {
      for (uint64_t i = 0; i != N; ++i) {
        for (uint64_t j = 0; j != N; ++j) {
          swapUnit(rowBegin1 + i, colBegin1 + j, rowBegin2 + j, colBegin2 + i);
        }
      }
      return;
    }
    uint64_t halfN = N / 2;
    transposeAndSwap(rowBegin1, colBegin1, rowBegin2, colBegin2, halfN);
    transposeAndSwap(rowBegin1, colBegin1 + halfN, rowBegin2 + halfN, colBegin2,
                     halfN);
    transposeAndSwap(rowBegin1 + halfN, colBegin1, rowBegin2, colBegin2 + halfN,
                     halfN);
    transposeAndSwap(rowBegin1 + halfN, colBegin1 + halfN, rowBegin2 + halfN,
                     colBegin2 + halfN, halfN);
  }

  /**
   * transpose an N * N region
   */
  void transpose(uint64_t rowBegin, uint64_t colBegin, uint64_t N) {
    if (N <= cacheEfficientSideLen) {
      for (uint64_t i = 0; i != N; ++i) {
        for (uint64_t j = i + 1; j != N; ++j) {
          swapUnit(rowBegin + i, colBegin + j, rowBegin + j, colBegin + i);
        }
      }
      return;
    }
    uint64_t halfN = N / 2;
    transpose(rowBegin, colBegin, halfN);
    transposeAndSwap(rowBegin + halfN, colBegin, rowBegin, colBegin + halfN,
                     halfN);
    transpose(rowBegin + halfN, colBegin + halfN, halfN);
  }

 public:
  /** Specify the parameters for transposition
   * @param[in] begin is the beginning iterator of v to transpose
   * @param[in] end is the ending iterator of v to transpose
   * @param[in] unitSize is the number of element considered as a unit in the
   * transposition. In our application, it's either the bucket size or twice the
   * bucket size
   * @param[in] cacheSize is the size of cache in byte. Naive approach will be
   * applied if the (sub)problem fits in cacheSize.
   *
   */
  MatTransposer(Iterator _begin, Iterator _end, uint64_t _unitSize,
                uint64_t cacheSize)
      : unitSize(_unitSize), begin(_begin) {
    sideLen = 1UL << (GetLogBaseTwo((_end - _begin) / _unitSize) / 2);
    cacheEfficientSideLen =
        1UL << (GetLogBaseTwo(cacheSize / (sizeof(T) * _unitSize * 2)) / 2);
    Assert(sideLen * sideLen * _unitSize == _end - _begin);
  }

  /**
   * Get the side length of the square matrix
   */
  uint64_t getSideLen() const { return sideLen; }

  /**
   * transpose the vector according to the specifications obtained from the
   * constructor
   */
  void transpose() { return transpose(0, 0, sideLen); }
};
}  // namespace EM::Algorithm