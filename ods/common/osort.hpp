#pragma once

// Obliviously sort
// cmp(a,b) -> true iff a < b (accepting const refs)
template<typename T, typename Compare>
void ObliSort(std::vector<T>& guys, Compare cmp) {
  Assert(1 << CeilLog2(guys.size()) == guys.size());
  uint64_t n = guys.size();
  for (uint64_t k = 2; k <= n; k *= 2) {
    for (uint64_t j = k/2; j > 0; j /= 2) {
      for (uint64_t i = 0; i < n; ++i) {
        uint64_t l = i ^ j;
        if (l > i) {
          bool llti = cmp(guys[l], guys[i]);
          bool swap = 
            (((i & k) == 0) * llti) +
            (((i & k) != 0) * (!llti));
          CXCHG(swap, guys[i], guys[l]);
        }
      }
    }
  }
}

template <typename StashedBlock, typename Compare>
inline void ObliSortBlocks(std::vector<StashedBlock>& blocks, Compare cmp) {
  // Append dummies to blocks
  uint64_t currSize = blocks.size();
  uint64_t desiredSize = GetNextPowerOfTwo(currSize);
  for (uint64_t i = currSize; i < desiredSize; ++i) {
    blocks.push_back(MakeDummy<StashedBlock>());
  }
  // Sort
  ObliSort(blocks, cmp);
  // Resize
  blocks.resize(currSize);
}
