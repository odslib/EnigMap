#include "oram/ringoram/oram.hpp"
#include "oram/pathoram/oram.hpp"
#include "oram/notoram/oram.hpp"
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


template <typename T>
class TestOTree : public testing::Test {
  public:
  using ORAMClient_t = typename T::ORAMClient;
  using OramClient_t = typename _OBST::OramClient::OramClient<ORAMClient_t>;
  using OBST = typename _OBST::OBST::OBST<OramClient_t>;
};
TYPED_TEST_SUITE_P(TestOTree);

TYPED_TEST_P(TestOTree, BasicAssertions) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 1 + (random() % 100);
  const K maxK = size * 100;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);

  for (int i = 0; i < size; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      cerr << k << endl;
      k = (random() % maxK);
      v1 = 0;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }

    k = (random() % maxK);
    v1 = random();
    cout << i << " " << k << endl;
    t.Insert(k, v1);
    m[k] = v1;
  }

  // This asserts that insertions worked.
  //
  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimpleBalancingRot2) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  const K maxK = size;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);

  for (int i = 0; i < size; i++) {
    t.RecursivePrint();
    cout << i << " " << k << endl;
    for (int _ = 0; _ < 10; _++) {
      k = (random() % maxK);
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }

    k = i;
    v1 = random();
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimpleBalancingRot3) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  const K maxK = size;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);

  for (int i = 0; i < size; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      k = (random() % maxK);
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }

    k = (i % 2) * 2 * size + i;
    cout << i << " " << k << endl;
    v1 = random();
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimplestTest) {
  TTHEADER();
  const int size = 11;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);
  cout << endl;
  const int vals[] = {1, 3, 2, 10};
  for (int i = 0; i < 4; i++) {
    k = vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    cout << endl;

    t.RecursivePrint();
    cout << endl;
  }
}


TYPED_TEST_P(TestOTree, SimplestRot2) {
  TTHEADER();
  // This is actually doing a rot3 that is equivalent to a rot2
  //
  std::map<K, V> m;
  const int size = 11;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);

  const int vals[] = {1, 3, 2, 10};
  for (int i = 0; i < 4; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      k = _;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    k = vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimplestRot2_2) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);

  const int vals[] = {3, 1, 2, 100};
  for (int i = 0; i < 4; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      k = _;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    k = vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimplestBalanced) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);

  const int vals[] = {3, 2, 5, 1, 4, 6, 8, 7, 100};
  for (int i = 0; i < 9; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      k = _;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    k = vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimplestRot3) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);

  const int vals[] = {2, 1, 6, 4, 7, 5, 10};
  for (int i = 0; i < 7; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      k = _;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    k = vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimplestRot3_2) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);

  const int vals[] = {2, 1, 6, 4, 7, 5, 10};
  for (int i = 0; i < 7; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      k = _;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    k = 100 - vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimplestRot3_3) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);
  const int vals[] = {5, 2, 6, 1, 3, 4, 100};
  for (int i = 0; i < 7; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      k = _;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    k = vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimplestRot3_4) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);
  const int vals[] = {2, 1, 5, 4, 6, 3, 100};
  for (int i = 0; i < 7; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 10; _++) {
      k = _;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    k = vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}

TYPED_TEST_P(TestOTree, SimplestRot3_5) {
  TTHEADER();
  std::map<K, V> m;
  const int size = 100;
  OBST t(size);
  K k;
  V v1, v2;

  ASSERT_GE(t.maxNodes, size);
  const int vals[] = {8, 3, 10, 2, 5, 9, 11, 1, 4, 7, 6, 100};
  for (int i = 0; i < 12; i++) {
    t.RecursivePrint();
    for (int _ = 0; _ < 15; _++) {
      k = _;
      EXPECT_EQ(t.Get(k, v1), m.count(k) > 0);
      if (m.count(k) > 0) {
        v2 = m[k];
        EXPECT_EQ(v1, v2);
      }
    }
    k = vals[i];
    v1 = vals[i] * 10;
    t.Insert(k, v1);
    m[k] = v1;
  }

  t.RecursivePrint();
}



template<
    typename _ORAMClient>
struct TestParameter {
  using ORAMClient = _ORAMClient; 
};


typedef ::testing::Types<
    TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<_OBST::Node,ORAM__Z,false,4> >
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<_OBST::Node,ORAM__Z,true,4> >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<_OBST::Node,false,false> >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<_OBST::Node,false,true> >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<_OBST::Node,true,false> >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<_OBST::Node,true,true> >
> TestedTypes;

REGISTER_TYPED_TEST_SUITE_P(TestOTree,
  BasicAssertions,
  SimpleBalancingRot2,
  SimpleBalancingRot3,
  SimplestTest,
  SimplestRot2,
  SimplestRot2_2,
  SimplestBalanced,
  SimplestRot3,
  SimplestRot3_2,
  SimplestRot3_3,
  SimplestRot3_4,
  SimplestRot3_5
);
INSTANTIATE_TYPED_TEST_SUITE_P(OTREE, TestOTree, TestedTypes);
