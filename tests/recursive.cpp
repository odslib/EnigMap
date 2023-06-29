#include "oram/notoram/oram.hpp"
#include "oram/pathoram/oram.hpp"
#include "oram/ringoram/oram.hpp"
#include "recoram/recursive.hpp"
#include <gtest/gtest.h>
#include <map>
using namespace std;

// This file tests ORAM that doesn't need a position map,
// such as recursive path oram.
//
#ifdef NDEBUG
  #define ISDEBUG false
#else
  #define ISDEBUG true
#endif

#define DEBUG_ONLY_TEST() if (!ISDEBUG) { GTEST_SKIP(); } 


template<
    typename _BaseORAM
  , typename _PositionsORAM
  , _ORAM::Index _sz>
struct TestParameter {
  using BaseORAM = _BaseORAM;
  using PositionsORAM = _PositionsORAM;
  using ORAMClient_t = typename _ORAM::RecursiveORAM::RecursiveORAM<uint64_t, BaseORAM, PositionsORAM>;
  static constexpr _ORAM::Index sz = _sz;
};

#define TTHEADER() \
  using ORAMClient_t = typename TestFixture::ORAMClient_t; \
  using Block_t = typename ORAMClient_t::Block_t; \
  constexpr _ORAM::Index sz = TestFixture::sz;

template <typename T>
class TestRecursiveORAM : public testing::Test {
  public:
  using ORAMClient_t = typename T::ORAMClient_t;
  static constexpr uint64_t sz = T::sz;
};
TYPED_TEST_SUITE_P(TestRecursiveORAM);

TYPED_TEST_P(TestRecursiveORAM, BasicAssertions) {
  TTHEADER();
  ORAMClient_t client{10};
  _ORAM::Address addr = _ORAM::Address{(_ORAM::Address)UniformRandom(9)};
  uint64_t r = 1;
  client.Write(addr, r);
  ASSERT_EQ(r, 1);
  r = 0;
  client.Read(addr, r);
  ASSERT_EQ(r, 1);
  r = 2;
  client.Write(addr, r);
  ASSERT_EQ(r, 2);
  r = 3; 
  client.Write(addr, r);
  ASSERT_EQ(r, 3);
  client.Read(addr, r);
  ASSERT_EQ(r, 3);
}

TYPED_TEST_P(TestRecursiveORAM, BasicAssertion2) {
  TTHEADER();
  ORAMClient_t client{sz};
  std::map<_ORAM::Index, uint64_t> vals;

  uint64_t addresses[] = {0,1,2,3,4,5};
  for (int _=0; _<1111; _++) {
    for (int i=0; i<sizeof(addresses)/sizeof(uint64_t); i++) {
      bool isNew = false;
      _ORAM::Address addr;
      Block_t r;

      addr = addresses[i]%sz;
      if (vals.count(addr)==0) {
        isNew = true;
      }

      client.BeginAccess(addr, r);
      if (!isNew) {
        ASSERT_EQ(r.data, vals[addr]);
      }
      vals[addr] = UniformRandom(12345678);
      r.data = vals[addr];
      
      client.WriteToSomeStashedBlock(addr, r);
      ASSERT_EQ(r.data, vals[addr]);
      client.FinishAccess();
      ASSERT_EQ(r.data, vals[addr]);
    }
  }
}

TYPED_TEST_P(TestRecursiveORAM, BasicAssertionRepetitive) {
  TTHEADER();
  ORAMClient_t client{sz};
  std::map<_ORAM::Index, uint64_t> vals;

  for (int i=0; i<10000; i++) {
    bool isNew = false;
    _ORAM::Address addr;
    uint64_t r;

    addr = UniformRandom(sz-1);
    if (vals.count(addr)==0) {
      isNew = true;
    }

    if (UniformRandom(1) == 0) {
      client.Read(addr, r);
      if (!isNew) {
        ASSERT_EQ(r, vals[addr]);
      }
    }
    
    vals[addr] = UniformRandom(12345678);
    r = vals[addr];
    
    client.Write(addr, r);
    ASSERT_EQ(r, vals[addr]);
  }
}

// TYPED_TEST_P(TestRecursiveORAM, DeathTestOOBOldPositionBeginAccess) {
//   TTHEADER();
//   DEBUG_ONLY_TEST();
//   ORAMClient_t client{sz};
//   Block_t tmpBlock = Block_t::DUMMY();
//   _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{static_cast<_ORAM::Position>(random())%sz, sz+1}; 
//   _ORAM::Position newPos = static_cast<_ORAM::Position>(random())%sz;
//   EXPECT_DEATH(client.BeginAccess(addr, newPos, tmpBlock), ".position");
// }
// TYPED_TEST_P(TestRecursiveORAM, DeathTestOOBAddressBeginAccess) {
//   TTHEADER();
//   DEBUG_ONLY_TEST();
//   ORAMClient_t client{sz};
//   Block_t tmpBlock = Block_t::DUMMY();
//   _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{sz+1,static_cast<_ORAM::Position>(random())%sz}; 
//   _ORAM::Position newPos = static_cast<_ORAM::Position>(random())%sz;
//   EXPECT_DEATH(client.BeginAccess(addr, newPos, tmpBlock), ".address");
// }
// TYPED_TEST_P(TestRecursiveORAM, DeathTestOOBNewPositionBeginAccess) {
//   TTHEADER();
//   DEBUG_ONLY_TEST();
//   ORAMClient_t client{sz};
//   Block_t tmpBlock = Block_t::DUMMY();
//   _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{static_cast<_ORAM::Position>(random())%sz, static_cast<_ORAM::Position>(random())%sz}; 
//   _ORAM::Position newPos = sz+1;
//   EXPECT_DEATH(client.BeginAccess(addr, newPos, tmpBlock), "newPos");
// }
// TYPED_TEST_P(TestRecursiveORAM, DeathTestOOBPositionForceIntoStash) {
//   TTHEADER();
//   DEBUG_ONLY_TEST();
//   ORAMClient_t client{sz};
//   _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{static_cast<_ORAM::Position>(random())%sz, sz+1}; 
//   EXPECT_DEATH(client.ForceIntoStash(addr, false), ".position");
// }
// TYPED_TEST_P(TestRecursiveORAM, DeathTestOOBAddressForceIntoStash) {
//   TTHEADER();
//   DEBUG_ONLY_TEST();
//   ORAMClient_t client{sz};
//   _ORAM::ORAMAddress addr = _ORAM::ORAMAddress{sz+1,static_cast<_ORAM::Position>(random())%sz}; 
//   EXPECT_DEATH(client.ForceIntoStash(addr, false), ".address");
// }


#define BASE (1<<10)
typedef ::testing::Types<
    TestParameter<
          _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>
        , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,false,false>
        , 10 >
  , TestParameter<
        _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,true>
      , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,false,true>
      , 10 >
  , TestParameter<
        _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,true,false>
      , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,true,false>
      , 10 >
  , TestParameter<
        _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,true,true>
      , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,true,true>
      , 10 >
  , TestParameter<
        _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>
      , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,false,false>
      , BASE-2 >
  , TestParameter<
        _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>
      , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,false,false>
      , BASE-1 >
  , TestParameter<
        _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>
      , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,false,false>
      , BASE >
  , TestParameter<
        _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>
      , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,false,false>
      , BASE+1 >
  , TestParameter<
        _ORAM::NotORAM::ORAMClient::ORAMClient<uint64_t,false,false>
      , _ORAM::NotORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,false,false>
      , BASE+2 >

  , TestParameter<
        _ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>
      , _ORAM::RingORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,ORAM__S,false,false,4>
      , 10>
  , TestParameter<
        _ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>
      , _ORAM::RingORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,ORAM__S,false,false,4>
      , BASE-2>
  , TestParameter<
        _ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>
      , _ORAM::RingORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,ORAM__S,false,false,4>
      , BASE-1>
  , TestParameter<
        _ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>
      , _ORAM::RingORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,ORAM__S,false,false,4>
      , BASE>
  , TestParameter<
        _ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>
      , _ORAM::RingORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,ORAM__S,false,false,4>
      , BASE+1>
  , TestParameter<
        _ORAM::RingORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,ORAM__S,false,false,4>
      , _ORAM::RingORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,ORAM__S,false,false,4>
      , BASE+2>

  , TestParameter<
        _ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,false,4>
      , _ORAM::PathORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,false,false,4>
      , 10>
  , TestParameter<
        _ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,false,4>
      , _ORAM::PathORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,false,false,4>
      , BASE-2>
  , TestParameter<
        _ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,false,4>
      , _ORAM::PathORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,false,false,4>
      , BASE-1>
  , TestParameter<
        _ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,false,4>
      , _ORAM::PathORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,false,false,4>
      , BASE>
  , TestParameter<
        _ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,false,4>
      , _ORAM::PathORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,false,false,4>
      , BASE+1>
  , TestParameter<
        _ORAM::PathORAM::ORAMClient::ORAMClient<uint64_t,ORAM__Z,false,false,4>
      , _ORAM::PathORAM::ORAMClient::ORAMClient<_ORAM::RecursiveORAM::Data,ORAM__Z,false,false,4>
      , BASE+2>
> TestedTypes;

REGISTER_TYPED_TEST_SUITE_P(
    TestRecursiveORAM
  , BasicAssertions
  , BasicAssertion2
  , BasicAssertionRepetitive
  //, DeathTestOOBOldPositionBeginAccess
  //, DeathTestOOBAddressBeginAccess
  //, DeathTestOOBNewPositionBeginAccess
  //, DeathTestOOBPositionForceIntoStash
  //,  DeathTestOOBAddressForceIntoStash
);
INSTANTIATE_TYPED_TEST_SUITE_P(RecORAM, TestRecursiveORAM, TestedTypes);
