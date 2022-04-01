#include "common/utils.hpp"
#include "oram/common/indexers.hpp"
#include <gtest/gtest.h>
#include <map>
using namespace _ORAM;
using namespace std;

TEST(Indexers, GetNextPowerOfTwo) {
  ASSERT_EQ(0, GetNextPowerOfTwo(0));
  ASSERT_EQ(1, GetNextPowerOfTwo(1));
  ASSERT_EQ(2, GetNextPowerOfTwo(2));
  ASSERT_EQ(4, GetNextPowerOfTwo(3));
  ASSERT_EQ(4, GetNextPowerOfTwo(4));
  ASSERT_EQ(0x8000'0000'0000'0000, GetNextPowerOfTwo(0x8000'0000'0000'0000));
  for (uint64_t i=2; i<63; i++) {
    for (uint64_t _=0; _<1000; _++) {
      uint64_t test = 1ULL + (1ULL<<i) + (random()%(1<<i));
      // cout << i << " " << test << endl;
      ASSERT_EQ((1ULL<<(i+1ULL)), GetNextPowerOfTwo(test));
    }
  }
}


TEST(Indexers, GetLogBaseTwo) {
  ASSERT_EQ(0, GetLogBaseTwo(0));
  ASSERT_EQ(0, GetLogBaseTwo(1));
  ASSERT_EQ(1, GetLogBaseTwo(2));
  ASSERT_EQ(1, GetLogBaseTwo(3));
  ASSERT_EQ(2, GetLogBaseTwo(4));
  ASSERT_EQ(63, GetLogBaseTwo(0x8000'0000'0000'0000));
  for (uint64_t i=0; i<63; i++) {
    for (uint64_t _=0; _<1000; _++) {
      uint64_t test = (1ULL<<i) + (random()%(1ULL<<i));
      ASSERT_EQ(i, GetLogBaseTwo(test));
    }
  }
}

TEST(Indexers, CeilLog2) {
  ASSERT_EQ(0, CeilLog2(0));
  ASSERT_EQ(0, CeilLog2(1));
  ASSERT_EQ(1, CeilLog2(2));
  ASSERT_EQ(2, CeilLog2(3));
  ASSERT_EQ(2, CeilLog2(4));
  ASSERT_EQ(63, CeilLog2(0x8000'0000'0000'0000));
  for (uint64_t i=0; i<63; i++) {
    for (uint64_t _=0; _<1000; _++) {
      uint64_t test = (1ULL<<i) + 1ULL + (random()%(1ULL<<i));
      ASSERT_EQ(i+1ULL, CeilLog2(test));
    }
  }
}

// UNDONE(test): ToReverseLexicographicalOrder

TEST(Indexers, GetArrIndex) {
  for (Index l=2; l<=10; l++) {
    uint64_t currInd = static_cast<uint64_t>(-1);
    uint64_t step = (1ULL<<l);
    for (Index depth=0; depth<=l; depth++) {
      step = ((1ULL<<(l-depth)));
      for (Position pos=0; pos < (1ULL<<l); pos++) {
        if (pos % step == 0) currInd++;
        // cout << depth << " " << step << " " << pos << " " << currInd << endl;
        ASSERT_EQ(currInd, _ORAM::Indexers::GetArrIndex(l, pos, depth));
      }
    }
    ASSERT_EQ(step, 1);
    ASSERT_EQ(currInd, (2ULL<<l)-2ULL);
  }  
}


// Asserts that if we do GetArrIndex(GetPosDepthFromIndex(GettArrIndex)) == GetArrIndex()
// It is NOT true that GetPosDepthFromIndex(GettArrIndex) = I
//
TEST(Indexers, GetPosDepthFromIndex) {
  constexpr Index LEVELS_PER_PACK = 2;
  constexpr Index BUCKETS_PER_PACK = (1<<LEVELS_PER_PACK)-1;
  for (Index l=2; l<=10; l++) {
    uint64_t currInd = static_cast<uint64_t>(-1);
    uint64_t step = (1ULL<<l);
    for (Index depth=0; depth<=l; depth++) {
      step = ((1ULL<<(l-depth)));
      for (Position pos=0; pos < (1ULL<<l); pos++) {
        if (pos % step == 0) currInd++;
        ASSERT_EQ(currInd, _ORAM::Indexers::GetArrIndex(l, pos, depth));
        Position compPos;
        Index compDepth;
        _ORAM::Indexers::GetPosDepthFromIndex(l, currInd, compPos, compDepth);
        // cout << currInd << " " << depth << " " << pos << " " << compDepth << " " << compPos << endl;
        ASSERT_EQ(depth, compDepth);
        ASSERT_EQ(currInd, _ORAM::Indexers::GetArrIndex(l, compPos, compDepth));
      }
    }
    ASSERT_EQ(step, 1);
    ASSERT_EQ(currInd, (2ULL<<l)-2ULL);
  }
}

// UNDONE(test): GetHBIndex
// UNDONE(test): GetLBIndex
// UNDONE(test): PathsIntercept