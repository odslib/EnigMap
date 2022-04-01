#include "oram/common/oram_client_interface.hpp"
#include "otree/node.hpp"
#include <gtest/gtest.h>
#include <map>
using namespace _ORAM;
using namespace std;

namespace NS1 {
template<uint64_t size>
requires (size > 0)
struct VariableSizedStruct {
  uint8_t v[size];
};

template <typename T>
struct Perf_Encryption : public testing::Test{
  using VariableSizedStruct_t = VariableSizedStruct<T::BlockSize>;
  static_assert(sizeof(VariableSizedStruct_t) == T::BlockSize);
};
#define TTHEADER() \
  using VariableSizedStruct_t = typename TestFixture::VariableSizedStruct_t; \
  using Encrypted_t = Encrypted<VariableSizedStruct_t>; \
  using NonEncrypted_t = NonEncrypted<VariableSizedStruct_t>; 

TYPED_TEST_SUITE_P(Perf_Encryption);

TYPED_TEST_P(Perf_Encryption, Encryption) {
  TTHEADER();
  Encrypted_t emd;
  constexpr int count = 1000000;

  VariableSizedStruct_t data1;
  data1.v[0] = 66;
  VariableSizedStruct_t data2;
  data2.v[0] = 88;

  const clock_t t0 = clock();

  for (int i=0; i<count; i++) {
    emd.Encrypt(data1);
    emd.Encrypt(data2);
  }

  const clock_t t1 = clock();

  double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
  cout << endl;
  cout << sizeof(VariableSizedStruct_t) << " encrypted bytes/op" << endl;
  cout << 2*count << " ops" << endl;
  cout << s << " total time" << endl;
  cout << s / (2*count) * 1'000'000 << " us/op" << endl;
  cout << s / (2*count*sizeof(VariableSizedStruct_t)) * 1'000'000 << " us/byte" << endl;
}


TYPED_TEST_P(Perf_Encryption, Decryption) {
  TTHEADER();
 
  Encrypted_t emd1;
  Encrypted_t emd2;
  constexpr int count = 1000000;

  VariableSizedStruct_t data1;
  data1.v[0] = 66;
  VariableSizedStruct_t data2;
  data2.v[0] = 88;
  emd1.Encrypt(data1);
  emd2.Encrypt(data2);

  const clock_t t0 = clock();

  for (int i=0; i<count; i++) {
    emd1.Decrypt(data1);
    emd2.Decrypt(data1);
  }

  const clock_t t1 = clock();

  double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
  cout << endl;
  cout << sizeof(VariableSizedStruct_t) << " decrypted bytes/op" << endl;
  cout << 2*count << " ops" << endl;
  cout << s << " total time" << endl;
  cout << s / (2*count) * 1'000'000 << " us/op" << endl;
  cout << s / (2*count*sizeof(VariableSizedStruct_t)) * 1'000'000 << " us/byte" << endl;
}

TYPED_TEST_P(Perf_Encryption, Mov) {
  TTHEADER();
 
  constexpr int count = 100000000;

  VariableSizedStruct_t data1;
  data1.v[0] = 66;
  volatile VariableSizedStruct_t data2;
  data2.v[0] = 88;

  const clock_t t0 = clock();

  for (int i=0; i<count; i++) {
    if constexpr (sizeof(VariableSizedStruct_t) == 1) {
      data2.v[0] = data1.v[0];
    } else if constexpr (sizeof(VariableSizedStruct_t) == 2) {
      *((volatile uint16_t*)&(data2.v[0])) = *((uint16_t*)&(data2.v[0]));
    } else if constexpr (sizeof(VariableSizedStruct_t) == 4) {
      *((volatile uint32_t*)&(data2.v[0])) = *((uint32_t*)&(data2.v[0]));
    } else {
      for (int w=0; w<sizeof(data1); w+=8) {
        *((volatile uint64_t*)&(data2.v[w])) = *((uint64_t*)&(data2.v[w]));
      }
    } 
  }

  const clock_t t1 = clock();

  double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
  cout << endl;
  cout << sizeof(VariableSizedStruct_t) << " decrypted bytes/op" << endl;
  cout << count << " ops" << endl;
  cout << s << " total time" << endl;
  cout << s / (count) * 1'000'000 << " us/op" << endl;
  cout << s / (count*sizeof(VariableSizedStruct_t)) * 1'000'000 << " us/byte" << endl;
}

TYPED_TEST_P(Perf_Encryption, CMov) {
  TTHEADER();
 
  constexpr int count = 100000000;

  VariableSizedStruct_t data1;
  data1.v[0] = 66;
  volatile VariableSizedStruct_t data2;
  data2.v[0] = 88;

  const clock_t t0 = clock();

  for (int i=0; i<count; i++) {
    if constexpr (sizeof(VariableSizedStruct_t) == 1) {
      CMOV(true, (uint8_t&)data2.v[0], data1.v[0]);
    } else if constexpr (sizeof(VariableSizedStruct_t) == 2) {
      CMOV(true, *((uint16_t*)&(data2.v[0])), *((uint16_t*)&(data2.v[0])));
    } else if constexpr (sizeof(VariableSizedStruct_t) == 4) {
      CMOV(true, *((uint32_t*)&(data2.v[0])), *((uint32_t*)&(data2.v[0])));
    } else {
      for (int w=0; w<sizeof(data1); w+=8) {
        CMOV(true, *((uint64_t*)&(data2.v[w])), *((uint64_t*)&(data2.v[w])));
      }
    } 
  }

  const clock_t t1 = clock();

  double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
  cout << endl;
  cout << sizeof(VariableSizedStruct_t) << " decrypted bytes/op" << endl;
  cout << count << " ops" << endl;
  cout << s << " total time" << endl;
  cout << s / (count) * 1'000'000 << " us/op" << endl;
  cout << s / (count*sizeof(VariableSizedStruct_t)) * 1'000'000 << " us/byte" << endl;
}

template<uint64_t _BlockSize>
struct TestParameter {
  static constexpr uint64_t BlockSize = _BlockSize;
};

TYPED_TEST_P(Perf_Encryption, Swap) {
  TTHEADER();
 
  constexpr int count = 1000000;

  volatile VariableSizedStruct_t data1;
  data1.v[0] = 66;
  std::fstream lbios;
  lbios.open("/tmp/example.txt", std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::trunc);
  ASSERT_TRUE(lbios.is_open());

  const clock_t t0 = clock();
  for (int i=0; i<count; i++) {
    lbios.seekp(0);
    lbios.write((char*)&data1, sizeof(VariableSizedStruct_t));
    lbios.flush();
    lbios.seekg(0);
    lbios.read((char*)&data1, sizeof(VariableSizedStruct_t));
    lbios.sync();
  }

  const clock_t t1 = clock();

  double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
  cout << endl;
  cout << sizeof(VariableSizedStruct_t) << " decrypted bytes/op" << endl;
  cout << count << " ops" << endl;
  cout << s << " total time" << endl;
  cout << s / (count) * 1'000'000 << " us/op" << endl;
  cout << s / (count*sizeof(VariableSizedStruct_t)) * 1'000'000 << " us/byte" << endl;
}


typedef ::testing::Types<
      TestParameter<1>
    , TestParameter<2>
    , TestParameter<4>
    , TestParameter<8>
    , TestParameter<16>
    , TestParameter<32>
    , TestParameter<64>
    , TestParameter<128>
    , TestParameter<256>
    , TestParameter<512>
    , TestParameter<1024>
    , TestParameter<2048>
    , TestParameter<4096>
> TestedTypes;

REGISTER_TYPED_TEST_SUITE_P(Perf_Encryption,
                            Encryption
                            , Decryption
                            , Mov
                            , CMov
                            , Swap
                            );
INSTANTIATE_TYPED_TEST_SUITE_P(TT, Perf_Encryption, TestedTypes);

#undef TTHEADER
}

namespace NS2 {


template <typename T>
class BasicPerf : public testing::Test {
  public:
  using Block_t = typename Block::Block<typename T::Datatype, T::EncryptBlocks>;
  using Bucket_t = typename Bucket::Bucket<Block_t, T::Z, T::S, T::EncryptBlocks>;
  using ORAMServer_t = _ORAM::ORAMClientInterface::ORAMClientInterface<Block_t, Bucket_t,T::EncryptLargeBuckets,T::LargeBucketSize,T::DirectlyCachedLevels>;
  using LargeBucket_t = typename ORAMServer_t::LargeBucket_t;
  using StashedBlock_t = typename StashedBlock::StashedBlock<Block_t>;
  using BucketMetadata_t = typename Bucket_t::BucketMetadata_t;
};
#define TTHEADER() \
  using ORAMServer_t = typename TestFixture::ORAMServer_t; \
  using LargeBucket_t = typename TestFixture::LargeBucket_t; \
  using Block_t = typename TestFixture::Block_t; \
  using StashedBlock_t = typename TestFixture::StashedBlock_t; \
  using Bucket_t = typename TestFixture::Bucket_t; \
  using BucketMetadata_t = typename TestFixture::BucketMetadata_t;
TYPED_TEST_SUITE_P(BasicPerf);

TYPED_TEST_P(BasicPerf, PageEviction) {
  TTHEADER();
 
  constexpr int count = 100000;
  constexpr int depth = 18;
  constexpr int N = (1<<depth);  


  ORAMServer_t server(N);

  const clock_t t0 = clock();

  for (int i=0; i<count; i++) {
    int idx = random() % N;
    Bucket_t& bucket = server.GetBucketRef(idx, depth);
    std::ignore = bucket;
  }

  const clock_t t1 = clock();

  double s = (t1 - t0) / (double)CLOCKS_PER_SEC;
  cout << endl;
  cout << sizeof(typename LargeBucket_t::Encrypted_t) << " evicted bytes/op" << endl;
  cout << count << " ops" << endl;
  cout << s << " total time" << endl;
  cout << s / (count) * 1'000'000 << " us/op" << endl;
  cout << s / (count*sizeof(typename LargeBucket_t::Encrypted_t)) * 1'000'000 << " us/byte" << endl;
}


template<
    typename _Datatype 
  , int _Z
  , int _S
  , bool _EncryptBlocks
  , bool _EncryptLargeBuckets
  , int _LargeBucketSize
  , int _DirectlyCachedLevels>
struct TestParameter {
  using Datatype = _Datatype; 
  static constexpr int Z = _Z;
  static constexpr int S = _S;
  static constexpr bool EncryptBlocks = _EncryptBlocks;
  static constexpr bool EncryptLargeBuckets = _EncryptLargeBuckets;
  static constexpr int LargeBucketSize = _LargeBucketSize;
  static constexpr int DirectlyCachedLevels = _DirectlyCachedLevels;
};


typedef ::testing::Types<
    TestParameter<_OBST::Node,ORAM__Z,0,false,false,1,20>
  , TestParameter<_OBST::Node,ORAM__Z,0,false,false,2,20>
  , TestParameter<_OBST::Node,ORAM__Z,0,false,false,3,21>
  , TestParameter<_OBST::Node,ORAM__Z,0,false,false,4,20>

  , TestParameter<_OBST::Node,ORAM__Z,0,false,true,1,20>
  , TestParameter<_OBST::Node,ORAM__Z,0,false,true,2,20>
  , TestParameter<_OBST::Node,ORAM__Z,0,false,true,3,21>
  , TestParameter<_OBST::Node,ORAM__Z,0,false,true,4,20>

  , TestParameter<_OBST::Node,ORAM__Z,0,true,false,1,20>
  , TestParameter<_OBST::Node,ORAM__Z,0,true,false,2,20>
  , TestParameter<_OBST::Node,ORAM__Z,0,true,false,3,21>
  , TestParameter<_OBST::Node,ORAM__Z,0,true,false,4,20>

> TestedTypes;

REGISTER_TYPED_TEST_SUITE_P(BasicPerf,
                            PageEviction);
INSTANTIATE_TYPED_TEST_SUITE_P(TT, BasicPerf, TestedTypes);

}