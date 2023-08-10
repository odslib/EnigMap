#include "oram/notoram/oram.hpp"
#include "oram/pathoram/oram.hpp"
#include "oram/ringoram/oram.hpp"
#include <gtest/gtest.h>
#include <map>
#include "testutils.hpp"

using namespace std;

#define TTHEADER() \
  using ORAMClient_t = typename TestFixture::ORAMClient_t; \
  using Block_t = typename TestFixture::Block_t; \
  constexpr _ORAM::Index sz = TestFixture::sz;

template <typename T>
class TestORAM : public testing::Test {
  public:
  using ORAMClient_t = typename T::ORAMClient;
  using Block_t = typename ORAMClient_t::Block_t;
  static constexpr uint64_t sz = T::sz;
};
TYPED_TEST_SUITE_P(TestORAM);

TYPED_TEST_P(TestORAM, BasicAssertions) {
  TTHEADER();
  ORAMClient_t client{sz};
  _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{0,0};
  Block_t r;
  r.data = 1;
  client.ForceIntoStash(addr, false);
  client.FinishAccess();
  client.BeginAccess(addr, 0, r);
  r.data = 2;
  client.WriteToSomeStashedBlock(addr.address, r);
  ASSERT_EQ(r.data, 2);
  r.data = 3;
  client.FinishAccess();
  ASSERT_EQ(r.data, 3);

  r.data = 4;
  client.BeginAccess(addr, 0, r);
  ASSERT_EQ(r.data, 2);

  r.data = 5;
  client.WriteToSomeStashedBlock(addr.address, r);
  ASSERT_EQ(r.data, 5);
  r.data = 6;
  client.FinishAccess();
  ASSERT_EQ(r.data, 6);
}

TYPED_TEST_P(TestORAM, BasicAssertion2) {
  TTHEADER();
  ORAMClient_t client{sz};
  std::map<_ORAM::Index, uint64_t> vals;
  std::map<_ORAM::Index, _ORAM::Position> positions;

  // UNDONE(): Review this test (I think it has no randomness)
  //
  uint64_t addresses[] = {0,1,2,3,4,5};
  for (uint64_t _=0; _< 10000; _++) {
    for (int i=0; i<sizeof(addresses)/sizeof(uint64_t); i++) {
      bool isNew = false;
      _ORAM::ORAMAddress addr;
      Block_t r;

      addr.address = addresses[i]%sz;
      _ORAM::Position newPos = addresses[i]%sz;
      if (positions.count(addr.address)==0) {
        addr.position = addresses[i]%sz;
        client.ForceIntoStash(addr, false);
        client.FinishAccess();
        positions[addr.address] = addr.position;
        vals[addr.address] = uint64_t{UniformRandom(12345678)};
        isNew = true;
      }

      addr.position = positions[addr.address];
      client.BeginAccess(addr, newPos, r);
      if (!isNew) {
        ASSERT_EQ(r.data, vals[addr.address]);
      }
      r.data = addresses[i];
      vals[addr.address] = r.data;
      client.WriteToSomeStashedBlock(addr.address, r);
      ASSERT_EQ(r.data, vals[addr.address]);
      client.FinishAccess();
      ASSERT_EQ(r.data, vals[addr.address]);
    }
  }
}

TYPED_TEST_P(TestORAM, BasicAssertion3) {
  TTHEADER();
  ORAMClient_t client{sz};
  std::map<_ORAM::Index, uint64_t> vals;
  std::map<_ORAM::Index, _ORAM::Position> positions;

  // UNDONE(): Review this test (I think it has no randomness)
  //
  uint64_t addresses[] = {0,1,0,1,0,1,0,0};
  uint64_t randompos[] = {0,1,sz-3,sz-1, 0, 1, 1, 1};
  for (uint64_t _=0; _< 2; _++) {
    for (int i=0; i<sizeof(addresses)/sizeof(uint64_t); i++) {
      bool isNew = false;
      _ORAM::ORAMAddress addr;
      Block_t r;

      addr.address = addresses[i]%sz;
      cout << "Addr=" << addr.address << " (";

      _ORAM::Position newPos = randompos[i]%sz;
      if (positions.count(addr.address)==0) {
        cout << "new";
        addr.position = randompos[i]%sz;
        client.ForceIntoStash(addr, false);
        client.FinishAccess();
        positions[addr.address] = addr.position;
        vals[addr.address] = uint64_t{UniformRandom(12345678)};
        isNew = true;
      }
      addr.position = positions[addr.address];
      cout << ", pos=" << addr.position << ", newPos=" << newPos << ")" << endl;
      client.BeginAccess(addr, newPos, r);
      if (!isNew) {
        ASSERT_EQ(r.data, vals[addr.address]);
      }
      r.data = addresses[i];
      vals[addr.address] = r.data;
      client.WriteToSomeStashedBlock(addr.address, r);
      ASSERT_EQ(r.data, vals[addr.address]);
      client.FinishAccess();
      ASSERT_EQ(r.data, vals[addr.address]);
    }
  }
}

TYPED_TEST_P(TestORAM, Bug1) {
  TTHEADER();
  ORAMClient_t client{17};
  std::map<_ORAM::Index, uint64_t> vals;
  std::map<_ORAM::Index, _ORAM::Position> positions;
  
  #define FIS(isNew,Addr,pos,newPos,val) { \
    _ORAM::ORAMAddress addr{Addr, pos}; \
    Block_t r; \
    cout << "Addr=" << addr.address << " (" <<  ", pos=" << addr.position << ", newPos=" << newPos << ")" << endl; \
    if constexpr (isNew) { \
      client.ForceIntoStash(addr, /*cached=*/ false); \
      client.FinishAccess(); \
      positions[addr.address] = addr.position; \
      vals[addr.address] = uint64_t{val}; \
    } \
    addr.position = positions[addr.address]; \
    positions[addr.address] = newPos; \
    client.BeginAccess(addr, newPos, r); \
    if (!isNew) { \
      ASSERT_EQ(r.data, vals[addr.address]); \
    } \
    r.data = val; \
    vals[addr.address] = val; \
    client.WriteToSomeStashedBlock(addr.address, r); \
    ASSERT_EQ(r.data, vals[addr.address]); \
    client.FinishAccess(); \
    ASSERT_EQ(r.data, vals[addr.address]); \
  }

  FIS(true, 1, 9, 10, 100);
  FIS(true, 2, 3, 15, 200);
  FIS(true, 3, 5, 11, 300);
  FIS(false, 3, 11, 3, 400);
  FIS(true, 4, 13, 5, 500);
  FIS(false, 1, 10, 1, 600);
  FIS(true, 5, 16, 1, 700);
  FIS(true, 6, 14, 10, 800);
  FIS(false, 5, 1, 2, 900);
}

TYPED_TEST_P(TestORAM, BasicAssertionRepetitive) {
  TTHEADER();
  ORAMClient_t client{sz};
  std::map<_ORAM::Index, uint64_t> vals;
  std::map<_ORAM::Index, _ORAM::Position> positions;

  for (int i=0; i<10000; i++) {
    bool isNew = false;
    _ORAM::ORAMAddress addr;
    Block_t r;

    _ORAM::Position newPos = random()%sz;
    addr.address = random()%sz;

    cout << "Addr=" << addr.address << " (";
    if (positions.count(addr.address)==0) {
      addr.position = random()%sz;
      cout << "new";
      client.ForceIntoStash(addr, /*cached=*/ false);
      client.FinishAccess();
      positions[addr.address] = addr.position;
      vals[addr.address] = uint64_t{UniformRandom(12345678)};
      isNew = true;
    }

    addr.position = positions[addr.address];
    cout << ", pos=" << addr.position << ", newPos=" << newPos << ")" << endl;
    positions[addr.address] = newPos;
    client.BeginAccess(addr, newPos, r);
    if (!isNew) {
      ASSERT_EQ(r.data, vals[addr.address]);
    }
    r.data = random();
    vals[addr.address] = r.data;
    client.WriteToSomeStashedBlock(addr.address, r);
    ASSERT_EQ(r.data, vals[addr.address]);
    client.FinishAccess();
    ASSERT_EQ(r.data, vals[addr.address]);
  }
}

TYPED_TEST_P(TestORAM, DeathTestOOBOldPositionBeginAccess) {
  TTHEADER();
  // UNDONE(35): figure out why death tests are taking too long with perf counters enabled
  GTEST_SKIP();
  DEBUG_ONLY_TEST();
  ORAMClient_t client{sz};
  Block_t tmpBlock = Block_t::DUMMY();
  _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{static_cast<_ORAM::Position>(random())%sz, sz+1}; 
  _ORAM::Position newPos = static_cast<_ORAM::Position>(random())%sz;
  EXPECT_DEATH(client.BeginAccess(addr, newPos, tmpBlock), ".position");
}
TYPED_TEST_P(TestORAM, DeathTestOOBAddressBeginAccess) {
  TTHEADER();
  // UNDONE(35): figure out why death tests are taking too long with perf counters enabled
  GTEST_SKIP();
  DEBUG_ONLY_TEST();
  ORAMClient_t client{sz};
  Block_t tmpBlock = Block_t::DUMMY();
  _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{sz+1,static_cast<_ORAM::Position>(random())%sz}; 
  _ORAM::Position newPos = static_cast<_ORAM::Position>(random())%sz;
  EXPECT_DEATH(client.BeginAccess(addr, newPos, tmpBlock), ".address");
}
TYPED_TEST_P(TestORAM, DeathTestOOBNewPositionBeginAccess) {
  TTHEADER();
  // UNDONE(35): figure out why death tests are taking too long with perf counters enabled
  GTEST_SKIP();
  DEBUG_ONLY_TEST();
  ORAMClient_t client{sz};
  Block_t tmpBlock = Block_t::DUMMY();
  _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{static_cast<_ORAM::Position>(random())%sz, static_cast<_ORAM::Position>(random())%sz}; 
  _ORAM::Position newPos = sz+1;
  EXPECT_DEATH(client.BeginAccess(addr, newPos, tmpBlock), "newPos");
}
TYPED_TEST_P(TestORAM, DeathTestOOBPositionForceIntoStash) {
  TTHEADER();
  // UNDONE(35): figure out why death tests are taking too long with perf counters enabled
  GTEST_SKIP();
  DEBUG_ONLY_TEST();
  ORAMClient_t client{sz};
  _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{static_cast<_ORAM::Position>(random())%sz, sz+1}; 
  EXPECT_DEATH(client.ForceIntoStash(addr, false), ".position");
}
TYPED_TEST_P(TestORAM, DeathTestOOBAddressForceIntoStash) {
  TTHEADER();
  // UNDONE(35): figure out why death tests are taking too long with perf counters enabled
  GTEST_SKIP();
  DEBUG_ONLY_TEST();
  ORAMClient_t client{sz};
  _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{sz+1,static_cast<_ORAM::Position>(random())%sz}; 
  EXPECT_DEATH(client.ForceIntoStash(addr, false), ".address");
}


template<
    typename _ORAMClient
  , _ORAM::Index _sz>
struct TestParameter {
  using ORAMClient = _ORAMClient;
  static constexpr _ORAM::Index sz = _sz;
};



#define BASE (1<<10)
typedef ::testing::Types<
    TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>, 10 >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,true>, 10 >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,true,false>, 10 >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,true,true>, 10 >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>, BASE-2 >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>, BASE-1 >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>, BASE >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>, BASE+1 >
  , TestParameter<_ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>, BASE+2 >

  , TestParameter<_ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>, 6>
  , TestParameter<_ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>, BASE-2>
  , TestParameter<_ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>, BASE-1>
  , TestParameter<_ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>, BASE>
  , TestParameter<_ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>, BASE+1>
  , TestParameter<_ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>, BASE+2>

  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, 17>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4,12,true>, 17>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE-2>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE-1>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE+1>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE+2>
  
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, 17>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4,12,true>, 17>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE-2>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE-1>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE+1>
  , TestParameter<_ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,4>, BASE+2>
> TestedTypes;

REGISTER_TYPED_TEST_SUITE_P(
  TestORAM,
  BasicAssertions,
  BasicAssertion2,
  BasicAssertion3,
  Bug1,
  BasicAssertionRepetitive,
  DeathTestOOBOldPositionBeginAccess,
  DeathTestOOBAddressBeginAccess,
  DeathTestOOBNewPositionBeginAccess,
  DeathTestOOBPositionForceIntoStash,
  DeathTestOOBAddressForceIntoStash
);
INSTANTIATE_TYPED_TEST_SUITE_P(ORAM, TestORAM, TestedTypes);
