#pragma once
#include "common/dummy.hpp"
#include "common/utils.hpp"
// Obliviously sort
// cmp(a,b) -> true iff a < b (accepting const refs)
template <typename T, typename Compare>
void ObliSortP2(std::vector<T>& guys, Compare cmp) {
  Assert(1 << CeilLog2(guys.size()) == guys.size());
  uint64_t n = guys.size();
  for (uint64_t k = 2; k <= n; k *= 2) {
    for (uint64_t j = k / 2; j > 0; j /= 2) {
      for (uint64_t i = 0; i < n; ++i) {
        uint64_t l = i ^ j;
        if (l > i) {
          bool llti = cmp(guys[l], guys[i]);
          bool swap = (((i & k) == 0) * llti) + (((i & k) != 0) * (!llti));
          CXCHG(swap, guys[i], guys[l]);
        }
      }
    }
  }
}

template <bool dir, typename T, typename Compare>
void _BitonicMerge(const Compare& cmp, std::vector<T>& Arr, uint64_t lo,
                   uint64_t n) {
  if (n <= 1) return;

  uint64_t m = GetNextPowerOfTwo(n) >> 1;
  for (uint64_t i = lo; i < lo + n - m; i++) {
    bool llti = cmp(Arr[i], Arr[i + m]);
    bool swap = dir ^ llti;
    CXCHG(swap, Arr[i], Arr[i + m]);
  }

  _BitonicMerge<dir, T, Compare>(cmp, Arr, lo, m);
  _BitonicMerge<dir, T, Compare>(cmp, Arr, lo + m, n - m);
}

// UNDONE(): Rewrite this as iterative version
template <bool dir, typename T, typename Compare>
void _BitonicSort(const Compare& cmp, std::vector<T>& Arr, uint64_t lo,
                  uint64_t n) {
  if (n <= 1) return;
  const uint64_t m = n >> 1;
  _BitonicSort<!dir>(cmp, Arr, lo, m);
  _BitonicSort<dir>(cmp, Arr, lo + m, n - m);
  _BitonicMerge<dir>(cmp, Arr, lo, n);
}

template <typename T, typename Compare>
void ObliSort(std::vector<T>& Arr, Compare cmp) {
  _BitonicSort<true>(cmp, Arr, 0, Arr.size());
}

template <typename StashedBlock, typename Compare>
inline void ObliSortBlocks(std::vector<StashedBlock>& blocks, Compare cmp) {
  ObliSort(blocks, cmp);
}
