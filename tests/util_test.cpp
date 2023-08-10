#include <gtest/gtest.h>

#include <iostream>
#include <map>

#include "common/osort.hpp"
#include "common/tracing/tracer.hpp"
#include "common/utils.hpp"
#include "external_memory/algorithm/param_select.hpp"
#include "oram/common/block.hpp"
#include "testutils.hpp"

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
  EXPECT_EQ(b1, b1_copy);
  EXPECT_EQ(b2, b2_copy);
  CXCHG(true, b2, b1);
  EXPECT_EQ(b1, b2_copy);
  EXPECT_EQ(b2, b1_copy);
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
  EXPECT_EQ(b1, b1_copy);
  EXPECT_EQ(b2, b2_copy);
  CMOV(true, b1, b2);
  EXPECT_EQ(b2, b2_copy);
  EXPECT_EQ(b1, b2_copy);
}

TEST(TestCMOV, ServerStashedBlockTest) {
  using StashedBlock_t = _ORAM::StashedBlock::StashedBlock<>;
  StashedBlock_t b1;
  StashedBlock_t b2;
  StashedBlock_t b1_copy;
  StashedBlock_t b2_copy;
  memset(&b1, 0, sizeof(StashedBlock_t));
  memset(&b1_copy, 0, sizeof(StashedBlock_t));
  memset(&b2, 0x42, sizeof(StashedBlock_t));
  memset(&b2_copy, 0x42, sizeof(StashedBlock_t));
  b2.cached = b2_copy.cached = 1;
  CMOV(false, b1, b2);
  EXPECT_EQ(b1, b1_copy);
  EXPECT_EQ(b2, b2_copy);
  CMOV(true, b1, b2);
  EXPECT_EQ(b2, b2_copy);
  EXPECT_EQ(b1, b2_copy);
}

TEST(TestObliSort, Uint64Test) {
  std::vector<uint64_t> v = {3, 2, 5, 6, 7, 8, 1, 4};
  ObliSort(v, std::less<int>());
  for (uint64_t i = 0; i < 8; ++i) {
    EXPECT_EQ(i + 1, v[i]);
  }
}

TEST(TestObliSortP2, ServerBlockTest) {
  using StashedBlock_t = _ORAM::StashedBlock::StashedBlock<>;
  std::vector<StashedBlock_t> v;
  v.reserve(128);
  for (uint64_t i = 0; i < v.capacity(); ++i) {
    v.push_back(StashedBlock_t::DUMMY());
  }
  srand(time(0));
  for (uint64_t i = 0; i < v.size(); ++i) {
    v[i].oaddress.address = rand() % 256;
  }
  auto cmp = [](const StashedBlock_t &A, const StashedBlock_t &B) {
    return A.oaddress.address < B.oaddress.address;
  };
  ObliSortP2(v, cmp);
  for (uint64_t i = 0; i < v.size() - 1; ++i) {
    EXPECT_LE(v[i].oaddress.address, v[i + 1].oaddress.address);
  }
}

TEST(TestObliSort, ServerBlockTest) {
  using StashedBlock_t = _ORAM::StashedBlock::StashedBlock<>;
  srand(time(0));

  for (int iter = 0; iter < 10000; iter++) {
    uint64_t sz = 1 + rand() % 1234;
    std::vector<StashedBlock_t> v;
    v.reserve(sz);
    for (uint64_t i = 0; i < v.capacity(); ++i) {
      v.push_back(StashedBlock_t::DUMMY());
    }
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
}

void ASSERT_ROUGH_EQUAL(double a, double b) {
  ASSERT_LE(abs(a - b) / max(abs(a), abs(b)), 1e-8);
}
using namespace EM::Algorithm;

TEST(TestProbCalc, TestBinomial) {
  ASSERT_ROUGH_EQUAL(pow(2, logCombin(2, 5)), 10);
  ASSERT_ROUGH_EQUAL(pow(2, logCombin(4, 5)), 5);
  ASSERT_ROUGH_EQUAL(pow(2, logCombin(0, 5)), 1);
  ASSERT_ROUGH_EQUAL(binomLogPmf(3, 10, 0.2), -2.3123903532651036);
  ASSERT_ROUGH_EQUAL(binomLogPmf(5, 1e6, 1e-7), -23.660814287355144);
  ASSERT_ROUGH_EQUAL(binomLogSf(3, 10, 0.2), -3.048425553831353);
  ASSERT_ROUGH_EQUAL(binomLogSf(5, 1e6, 1e-7), -29.546991149020247);
  ASSERT_ROUGH_EQUAL(binomLogCdf(3, 10, 0.2), -0.18585794733617342);
  ASSERT_ROUGH_EQUAL(binomLogCdf(5, 1e8, 1e-6), -117.88396022944082);
}

TEST(TestProbCalc, TestHypergeom) {
  ASSERT_ROUGH_EQUAL(hypergeomLogPmf(3, 100, 10, 20), -2.2569894866783438);
  ASSERT_ROUGH_EQUAL(hypergeomLogPmf(5, 1e8, 1e6, 10), -25.31451279737036);
  ASSERT_ROUGH_EQUAL(hypergeomLogSf(3, 100, 10, 20), -3.190050877769886);
  ASSERT_ROUGH_EQUAL(hypergeomLogSf(5, 1e8, 1e6, 10), -32.198576291056035);
}

TEST(TestProbCalc, TestBounds) {
  ASSERT_ROUGH_EQUAL(logCombin(12345, 123456), logCombin(12345, 123456));
  ASSERT_LE(logCombin(12345, 123456), logCombin(12345, 123456));
  ASSERT_LE(logCombin(12345, 123456), logCombin(12345, 123456));
}