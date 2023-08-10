#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <list>

#include "common/defs.hpp"
#include "common/utils.hpp"
#include "common/tracing/tracer.hpp"
#include "common/lrucache.hpp"

#include "oram/common/indexers.hpp"
#include "oram/common/block.hpp"
#include "oram/common/bucket.hpp"

#include "external_memory/server/serverFrontend.hpp"

#include "oram/common/oram_client_interface_assertions.hpp"

// This file implements the ORAM server, which turns regular indexes into
// eb-indexes in large buckets, and forward calls the storage class fileServer with the eb index.
//

namespace _ORAM::ORAMClientInterface {
  template<typename Block=_ORAM::Block::Block<>
    , typename Bucket=_ORAM::Bucket::Bucket<>
    , bool ENCRYPT_LARGE_BUCKETS=ORAM_SERVER__ENCRYPT_LARGE_BUCKETS
    , unsigned int LEVELS_PER_PACK=ORAM_SERVER__LEVELS_PER_PACK
    , unsigned int DIRECTLY_CACHED_LEVELS=ORAM_SERVER__DIRECTLY_CACHED_LEVELS>
  requires ::Concepts::Encryptable<LargeBucket::LargeBucket<Bucket,ENCRYPT_LARGE_BUCKETS,LEVELS_PER_PACK>>
  struct ORAMClientInterface {
    using Block_t = Block;
    using BucketMetadata_t = typename Bucket::BucketMetadata_t;
    using LargeBucket_t = LargeBucket::LargeBucket<Bucket,ENCRYPT_LARGE_BUCKETS,LEVELS_PER_PACK>;
    // UNDONE(): Read last level optimizations
    // using LargeBucket_LastLevel_t = LargeBucket::LargeBucket<typename Bucket::Encrypted_t,ENCRYPT_LARGE_BUCKETS,1>;
    // UNDONE(): paramterize ORAM_SERVER__BUCKET_CACHE_SIZE here
    Cache<LargeBucket_t,ORAM_SERVER__BUCKET_CACHE_SIZE> cache;
    
    static_assert( DIRECTLY_CACHED_LEVELS % LEVELS_PER_PACK == 0);
    static constexpr uint64_t DC_LARGE_BUCKET_LEVELS = DIRECTLY_CACHED_LEVELS / LEVELS_PER_PACK;

    uint64_t N_; // Number of blocks
    uint64_t L_; // Height of tree (where single root has height 0)
    uint64_t V_; // Number of buckets (well, technically # nodes + 1)

    // UNDONE(ttt): readd explicit treetop caching
    //
    // std::vector<LargeBucket_t> data;

    EM::MemoryServer::NonCachedServerFrontendInstance<LargeBucket_t, ::EM::Backend::MemServerBackend, true, false, false> server;

    OCI_ONLY(
      ORAMInterfaceVerifier<Block, Bucket> verifier;
    )

    explicit ORAMClientInterface(uint64_t N) :
      N_(N)
      , L_(CeilLog2(N))
      , V_(1ULL << ((L_ + 1 + LEVELS_PER_PACK - 1) / LEVELS_PER_PACK * LEVELS_PER_PACK))
      , server(*EM::Backend::g_DefaultBackend, V_/LargeBucket_t::BUCKETS_PER_PACK + 2, LargeBucket_t::DUMMY())
      OCI_ONLY(, verifier(L_))
    {
      TRACE_FUNCTION(N);
      PROFILE_F();
      // std::cerr << "oci.N_: " << N_ << std::endl;
      // std::cerr << "oci.L_: " << L_ << std::endl;
      // std::cerr << "oci.V_: " << V_ << std::endl;
      // UNDONE(ttt):
      // data.resize(1 << (DIRECTLY_CACHED_LEVELS/LEVELS_PER_PACK), LargeBucket_t());
    }
    

    // Copy not allowed!
    ORAMClientInterface& operator=(const ORAMClientInterface&) = delete;
    ORAMClientInterface(const ORAMClientInterface&) = delete;

    LargeBucket_t& GetLargeBucketRef(Position pos, const Index& depth) {
      Assert(pos <= V_);
      Assert(depth <= L_);

      Index rootDepth = depth - (depth % LEVELS_PER_PACK);
      Index rootIdx = Indexers::GetHBIndex<LEVELS_PER_PACK>(L_, pos, rootDepth);
      Index innerIdx = Indexers::GetLBIndex<LEVELS_PER_PACK>(L_, pos, depth);

      // UNDONE(ttt):
      //
      if (!cache.CheckContains(rootIdx)) {
        if (cache.IsFull()) {
          // UNDONE(): put correct typename here:
          uint64_t evictedIndex;
          {
            LargeBucket_t& evicted = cache.GetNextToEvict(evictedIndex).val;
            server.Write(evictedIndex, evicted);
          }
          cache.EvictLRU(evictedIndex);
        }
        // UNDONE(): How to make perfect forwarding here?
        //
        {
          LargeBucket_t lb;
          server.Read(rootIdx, lb);
          cache.Insert(rootIdx, lb);
        }
      }
      return cache.Access(rootIdx);
    }

    Bucket& GetBucketRef(Position pos, const Index& depth) {
      Assert(pos < V_);
      Assert(depth <= L_);
      // UNDONE(): hack to make noinit not crash:
      //
      if (pos >= V_) {
        // pos = pos % V_;
        pos = V_; // Special case to refer to void bucket.
      }
      Index innerIdx = Indexers::GetLBIndex<LEVELS_PER_PACK>(L_, pos, depth);
      LargeBucket_t& lb = GetLargeBucketRef(pos, depth);
      return lb.buckets[innerIdx];
    }

    void ReadBucketMetadata(const Index& pos, const Index& depth, BucketMetadata_t& ret) {
      const Bucket &bucket = GetBucketRef(pos, depth);
      ret = bucket.md;

      OCI_ONLY(verifier.AssertSameMetadata(pos,depth,ret));
    }
    
    void WriteBucketMetadata(const Index& pos, const Index& depth, const BucketMetadata_t& val) {
      Bucket &bucket = GetBucketRef(pos, depth);
      bucket.md = val;

      OCI_ONLY(verifier.UpdateMetadata(pos,depth,val));
    }

    void ReadBlock(const Index& pos, const Index& depth, const Index& offset, Block_t& ret) {
      Assert(offset < Bucket::BUCKET_SIZE);
      const Bucket &bucket = GetBucketRef(pos, depth);
      ret = bucket.blocks[offset];

      OCI_ONLY(verifier.AssertSameBlock(pos,depth,offset,ret));
    }

    void WriteBlock(const Index& pos, const Index& depth, const Index& offset, const Block_t& val) {
      Assert(offset < Bucket::BUCKET_SIZE);
      Bucket &bucket = GetBucketRef(pos, depth);
      bucket.blocks[offset] = val;

      OCI_ONLY(verifier.UpdateBlock(pos,depth,offset,val));
    }

    void _Dump() {
      for(uint64_t i=0; i<=L_; i++) {
        for(uint64_t j=0; j<(1<<i); j++) {
          Bucket& bkt = GetBucketRef(j*(1<<(L_-i)), i);
          std::cerr << "(" << i << ", " << j << "): " << bkt << std::endl;
        }
      }
    }
  };
  // </struct ORAMClientInterface>
} // namespace _ORAM::ORAMClientInterface