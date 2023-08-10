// #define NO_INMEM_SERVER
#include "oram/pathoram/oram.hpp"
#include "otree/otree.hpp"
#include <gtest/gtest.h>
#include <map>
using namespace _ORAM;
using namespace std;

using K = _OBST::K;
using V = _OBST::V;

#define TTHEADER() \
  using ORAMClient_t = typename TestFixture::ORAMClient_t; \
  using OramClient_t = typename TestFixture::OramClient_t; \
  using OBST = typename TestFixture::OBST;

template<
    typename _ORAMClient>
struct TestParameter {
  using ORAMClient = _ORAMClient; 
};
template <typename T>
class PerfOTree : public testing::Test {
  public:
  using ORAMClient_t = typename T::ORAMClient;
  using OramClient_t = typename _OBST::OramClient::OramClient<ORAMClient_t>;
  using OBST = typename _OBST::OBST::OBST<OramClient_t>;
};

typedef ::testing::Types<
    TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<_OBST::Node,ORAM__Z,false,1,10> >
> TestedTypes;

TYPED_TEST_SUITE(PerfOTree, TestedTypes);

TYPED_TEST(PerfOTree, PointSearchPartial) {
  TTHEADER();
  return;
  int sizes[27];
  const int testcases = sizeof(sizes) / sizeof(int); 
  for (int i=0; i<testcases; i++) {
    sizes[i] = (1<<(i+2))-3;
  }
  const int mult = 10;
  for (int iter = 0; iter < testcases; iter++) {
    const int size = 1 + sizes[iter];

    const K maxK = size * 100;
    std::map<K, V> m;
    OBST t(size);
    K k;
    V v1, v2;

    ASSERT_GE(t.maxNodes, size);

    TRACER_SET(false);

    // Make a few inserts to initially fill ORAM.
    //
    for (int i = 0; i < std::min(1000, size); i++) {
      k = (random() % maxK);
      v1 = random();
      t.Insert(k, v1);
      m[k] = v1;
    }
    
    TRACER_SET(true);

    const int cycles = 1000 * mult;
    const clock_t t0 = clock();
    for (int i = 0; i < cycles; i++) {
      k = (random() % maxK);
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    const clock_t t1 = clock();

    double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
    cout << endl;
    cout << "[Report]" << endl;
    cout << "max depth: " << t.maxDepth << endl;
    cout << "max nodes: " << t.maxNodes << endl;
    cout << "oram.N: " << t.oram.N_ << endl;
    cout << "oram.L: " << t.oram.L_ << endl;
    cout << "oram.Z: " << ORAM__Z << endl;
    cout << "queries: " << cycles << endl;
    cout << "total time: " << s << endl;
    cout << "us/query: " << s / cycles * 1'000'000 << endl;
    cout << endl;
    TRACER_RESET();
    cout << "[/Report]" << endl;
  }
}


TYPED_TEST(PerfOTree, InsertionsPartial) {
  TTHEADER();
  int sizes[27];
  const int testcases = sizeof(sizes) / sizeof(int); 
  for (int i=0; i<testcases; i++) {
    sizes[i] = (1<<(i+2))-5;
  }
  const int mult = 10;
  for (int iter = 5; iter < testcases; iter++) {
    const int size = 1 + sizes[iter];

    const K maxK = 1<<30;
    OBST t(size);
    K k;
    V v1, v2;

    ASSERT_GE(t.maxNodes, size);

    TRACER_SET(true);

    const int cycles = std::min(1000, size-4);
    const clock_t t0 = clock();
    for (int i = 0; i < cycles; i++) {
      k = (random() % maxK);
      v1 = random();
      t.Insert(k, v1);
    }
    const clock_t t1 = clock();

    double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
    cout << endl;
    cout << "[Report]" << endl;
    cout << "max depth: " << t.maxDepth << endl;
    cout << "max nodes: " << t.maxNodes << endl;
    cout << "oram.N: " << t.oram.N_ << endl;
    cout << "oram.L: " << t.oram.L_ << endl;
    cout << "oram.Z: " << ORAM__Z << endl;
    cout << "queries: " << cycles << endl;
    cout << "total time: " << s << endl;
    cout << "us/query: " << s / cycles * 1'000'000 << endl;
    cout << endl;
    TRACER_RESET();
    cout << "[/Report]" << endl;
  }
}
