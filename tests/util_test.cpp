#include "common/utils.hpp"
#include "oram/common/block.hpp"
#include "common/osort.hpp"
#include <gtest/gtest.h>
#include <map>
// using namespace _OBST;
// using namespace _ORAM;
using namespace std;

TEST(TestUtil, Log2) {
  EXPECT_EQ(CeilLog2(0), 0);
  EXPECT_EQ(CeilLog2(1), 0);
  EXPECT_EQ(CeilLog2(2), 1);
  EXPECT_EQ(CeilLog2(3), 2);
  EXPECT_EQ(CeilLog2(4), 2);
  EXPECT_EQ(CeilLog2(31), 5);
}

TEST(TestCXCHG, PrimitiveFunctionality) {
  uint64_t max = std::numeric_limits<uint64_t>::max();
  uint64_t a = 15213;
  uint64_t b = max;
  CXCHG(false, a, b);
  EXPECT_EQ(15213, a);
  EXPECT_EQ(max, b);
  CXCHG(true, a, b);
  EXPECT_EQ(max, a);
  EXPECT_EQ(15213, b);
}

TEST(TestCXCHG, Uint64Test) {
  uint64_t max = std::numeric_limits<uint64_t>::max();
  uint64_t a = 15213;
  uint64_t b = max;
  CXCHG(false, a, b);
  EXPECT_EQ(15213, a);
  EXPECT_EQ(max, b);
  CXCHG(true, a, b);
  EXPECT_EQ(max, a);
  EXPECT_EQ(15213, b);
}

TEST(TestCXCHG, BlockTest) {
  using Block_t = typename _ORAM::Block::Block<>;
  Block_t b1;
  Block_t b2;
  Block_t b1_copy;
  Block_t b2_copy;
  memset(&b1, 0, sizeof(Block_t));
  memset(&b1_copy, 0, sizeof(Block_t));
  memset(&b2, 0x42, sizeof(Block_t));
  memset(&b2_copy, 0x42, sizeof(Block_t));
  CXCHG(false, b1, b2);  
  EXPECT_EQ(b1,b1_copy);
  EXPECT_EQ(b2,b2_copy);
  CXCHG(true, b2, b1);
  EXPECT_EQ(b1,b2_copy);
  EXPECT_EQ(b2,b1_copy);
}

TEST(TestCMOV, PrimitiveFunctionality) {
  uint64_t max = std::numeric_limits<uint64_t>::max();
  uint64_t a = 15213;
  uint64_t b = max;
  CMOV(false, a, b);
  EXPECT_EQ(15213, a);
  EXPECT_EQ(max, b);
  CMOV(true, a, b);
  EXPECT_EQ(max, a);
  EXPECT_EQ(max, b);
}

TEST(TestCMOV, Uint64Test) {
  uint64_t max = std::numeric_limits<uint64_t>::max();
  uint64_t a = 15213;
  uint64_t b = max;
  CMOV(false, a, b);
  EXPECT_EQ(15213, a);
  EXPECT_EQ(max, b);
  CMOV(true, a, b);
  EXPECT_EQ(max, a);
  EXPECT_EQ(max, b);
}

TEST(TestCMOV, ServerBlockTest) {
  using Block_t = typename _ORAM::Block::Block<>;
  Block_t b1;
  Block_t b2;
  Block_t b1_copy;
  Block_t b2_copy;
  memset(&b1, 0, sizeof(Block_t));
  memset(&b1_copy, 0, sizeof(Block_t));
  memset(&b2, 0x42, sizeof(Block_t));
  memset(&b2_copy, 0x42, sizeof(Block_t));
  CMOV(false, b1, b2);
  EXPECT_EQ(b1,b1_copy);
  EXPECT_EQ(b2,b2_copy);
  CMOV(true, b1, b2);
  EXPECT_EQ(b2,b2_copy);
  EXPECT_EQ(b1,b2_copy);
}

TEST(TestObliSort, Uint64Test) {
  std::vector<uint64_t> v = {3, 2, 5, 6, 7, 8, 1, 4};
  ObliSort(v, std::less<int>());
  for (uint64_t i = 0; i < 8; ++i) {
    EXPECT_EQ(i + 1, v[i]);
  }
}

TEST(TestObliSort, ServerBlockTest) {
  using StashedBlock_t = _ORAM::StashedBlock::StashedBlock<>;
  std::vector<StashedBlock_t> v;
  v.reserve(128);
  for (uint64_t i = 0; i < v.capacity(); ++i) {
    v.push_back(MakeDummy<StashedBlock_t>());
  }
  srand(time(0));
  for (uint64_t i = 0; i < v.size(); ++i) {
    v[i].oaddress.address = rand() % 256;
  }
  auto cmp = [](const StashedBlock_t &A, const StashedBlock_t &B) {
    return A.oaddress.address < B.oaddress.address;
  };
  ObliSort(v, cmp);
  for (uint64_t i = 0; i < v.size() - 1; ++i) {
    EXPECT_LE(v[i].oaddress.address, v[i + 1].oaddress.address);
  }
}
