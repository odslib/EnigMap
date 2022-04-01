#include "oram/ringoram/oram.hpp"
#include "otree/otree.hpp"
#include <gtest/gtest.h>
#include <map>
using namespace _ORAM;
using namespace std;

using ORAMClient_t = typename _ORAM::RingORAM::ORAMClient::ORAMClient<_OBST::Node,ORAM__Z,ORAM__S,false,false,4>;
using OramClient_t = typename _OBST::OramClient::OramClient<ORAMClient_t>;
using OBST = typename _OBST::OBST::OBST<OramClient_t>;
using ORAMServer_t = typename ORAMClient_t::ORAMClientInterface_t;
using LargeBucket_t = typename ORAMServer_t::LargeBucket_t;
using Block_t = typename ORAMClient_t::Block_t;
using StashedBlock_t = typename ORAMClient_t::StashedBlock_t;
using Bucket_t = typename ORAMClient_t::Bucket_t;
using BucketMetadata_t = typename Bucket_t::BucketMetadata_t;

TEST(TestCache, ServerCacheAssertions) {
  ORAMServer_t server(4);
  BucketMetadata_t bmd1 = BucketMetadata_t::DUMMY();
  bmd1.pub.counter = 66;
  BucketMetadata_t bmd2 = BucketMetadata_t::DUMMY();
  bmd2.pub.counter = 88;

  server.WriteBucketMetadata(2,2,bmd1);
  server.WriteBucketMetadata(1,2,bmd2);
  server.ReadBucketMetadata(2,2,bmd2);
  ASSERT_EQ(bmd1, bmd2);
}

TEST(TestCache, ServerCacheAssertions2) {
  ORAMServer_t server(8);
  BucketMetadata_t bmd1 = BucketMetadata_t::DUMMY();
  bmd1.pub.counter = 66;
  BucketMetadata_t bmd2 = BucketMetadata_t::DUMMY();
  bmd2.pub.counter = 88;

  server.WriteBucketMetadata(2,0,bmd1);
  server.WriteBucketMetadata(7,0,bmd2);
  server.ReadBucketMetadata(2,0,bmd1);
  ASSERT_EQ(bmd1, bmd2);
}
