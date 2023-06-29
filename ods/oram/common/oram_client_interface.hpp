#pragma once
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <list>

#ifndef ENCLAVE_MODE
#include <boost/functional/hash.hpp>
#endif

#include "common/defs.hpp"
#include "common/utils.hpp"
#include "common/tracing/tracer.hpp"
#include "common/lrucache.hpp"

#include "oram/common/indexers.hpp"
#include "oram/common/block.hpp"
#include "oram/common/concepts.hpp"
#include "oram/common/bucket.hpp"

#ifndef ENCLAVE_MODE
#include "external_memory/server/fileServer.hpp"
#else
#include "external_memory/server/enclaveFileServer_trusted.hpp"
#endif

#include "external_memory/server/memServer.hpp"

#ifndef NDEBUG
#ifndef ENCLAVE_MODE
#define OCI_ASSERTIONS
#endif
#endif
// This file implements the ORAM server, which turns regular indexes into
// eb-indexes in large buckets, and forward calls the storage class
// FileServer with the eb index.
//

namespace _ORAM::ORAMClientInterface {
  // UNDONE(1): make cached data actually be a perfect forwarding class
  //
  template<typename Block=_ORAM::Block::Block<>
    , typename Bucket=_ORAM::Bucket::Bucket<>
    , bool ENCRYPT_LARGE_BUCKETS=ORAM_SERVER__ENCRYPT_LARGE_BUCKETS
    , unsigned int LEVELS_PER_PACK=ORAM_SERVER__LEVELS_PER_PACK
    , unsigned int DIRECTLY_CACHED_LEVELS=ORAM_SERVER__DIRECTLY_CACHED_LEVELS>
  requires ::Concepts::Encryptable<LargeBucket::LargeBucket<Bucket,ENCRYPT_LARGE_BUCKETS,LEVELS_PER_PACK>>
    && ::Concepts::Encryptable<Block>
    #ifndef ENCLAVE_MODE
    && Concepts::FileServer<
        EM::FileServer::FileServer
      ,  LargeBucket::LargeBucket<Bucket,ENCRYPT_LARGE_BUCKETS,LEVELS_PER_PACK>
      , ORAM_SERVER__BUCKET_CACHE_SIZE>
    #endif
  struct ORAMClientInterface {
    using Block_t = Block;
    using BucketMetadata_t = typename Bucket::BucketMetadata_t;
    using LargeBucket_t = LargeBucket::LargeBucket<typename Bucket::Encrypted_t,ENCRYPT_LARGE_BUCKETS,LEVELS_PER_PACK>;
    using LargeBucket_LastLevel_t = LargeBucket::LargeBucket<typename Bucket::Encrypted_t,ENCRYPT_LARGE_BUCKETS,1>;
    // UNDONE(): paramterize ORAM_SERVER__BUCKET_CACHE_SIZE here
    Cache<Bucket,ORAM_SERVER__BUCKET_CACHE_SIZE> cache;
    
    static_assert( DIRECTLY_CACHED_LEVELS % LEVELS_PER_PACK == 0);

    uint64_t N_; // Number of blocks
    uint64_t L_; // Height of tree (where single root has height 0)
    uint64_t V_; // Number of buckets (well, technically # nodes + 1)
    
    // UNDONE(): generate the key
    uint8_t key_[16] = {0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41,0x41};

    std::vector<Bucket> data;
    #ifndef ENCLAVE_MODE
      #ifndef NO_INMEM_SERVER
        EM::MemServer::MemServer<LargeBucket_t> server;
      #else
        EM::FileServer::FileServer<LargeBucket_t> server;
      #endif
    #else
      EM::MemServer::MemServer<LargeBucket_t> server;
      // EM::EnclaveFileServer::EnclaveFileServer<LargeBucket_t> server;
    #endif

    #ifdef OCI_ASSERTIONS
      typedef std::tuple<Position, Index, Index> Triplet;
      typedef std::tuple<Position, Index> Duplet;
      template<typename K>
      struct KeyHash {
          std::size_t operator()(const K &key) const
          {
              return boost::hash_value(key);
          }
      };

      std::unordered_map<Triplet, Block_t, KeyHash<Triplet> > checkdata_block;
      std::unordered_map<Duplet, BucketMetadata_t, KeyHash<Duplet> > checkdata_metadata;
    #endif
    

    void InitServer() {
      // UNDONE(): make sure this is always correct:
      //
      uint64_t bucketsToWrite = V_;
      uint64_t bytesToWrite = ((V_+LargeBucket_t::BUCKETS_PER_PACK)/LargeBucket_t::BUCKETS_PER_PACK)*sizeof(typename LargeBucket_t::Encrypted_t);
      X_LOG("ORAMClientInterface: Initializing file with ", bytesToWrite, " bytes (", V_, " buckets)");

      uint64_t numberOfLargeBuckets = bucketsToWrite/LargeBucket_t::BUCKETS_PER_PACK + LargeBucket_t::BUCKETS_PER_PACK;
      for (uint64_t i=0; i < numberOfLargeBuckets; i++) {
        LargeBucket_t lb;
        for (int j=0; j<LargeBucket_t::BUCKETS_PER_PACK; j++) {
          BucketMetadata_t md;
          Bucket b;
          for (int z=0; z<Bucket::BUCKET_SIZE; z++) {
            b.blocks[z] = Block_t::DUMMY();
          }
          b.md = Bucket::BucketMetadata_t::DUMMY();
          lb.buckets[j].Encrypt(b);
        }
        server.Write(i, lb);
      }
    }

    explicit ORAMClientInterface(uint64_t N, bool noInit=false) :
      N_(N),
      L_(CeilLog2(N)),
      V_(1ULL << ((L_ + 1 + LEVELS_PER_PACK - 1) / LEVELS_PER_PACK * LEVELS_PER_PACK)),
      server(V_/LargeBucket_t::BUCKETS_PER_PACK + LargeBucket_t::BUCKETS_PER_PACK)
    {
      TRACE_FUNCTION(N);
      PROFILE_F();
      // std::cerr << "oci.N_: " << N_ << std::endl;
      // std::cerr << "oci.L_: " << L_ << std::endl;
      // std::cerr << "oci.V_: " << V_ << std::endl;
      data.resize(1 << DIRECTLY_CACHED_LEVELS, Bucket());
      if (!noInit) {
        InitServer();
      }
    }
    

    // Copy not allowed!
    ORAMClientInterface& operator=(const ORAMClientInterface&) = delete;
    ORAMClientInterface(const ORAMClientInterface&) = delete;

    Bucket& GetBucketRef(Position pos, const Index& depth) {
      Assert(pos < V_);
      Assert(depth <= L_);
      // UNDONE(): hack to make noinit not crash:
      //
      if (pos >= V_) {
        pos = pos % V_;
      }
      Index rootDepth = depth - (depth % LEVELS_PER_PACK);
      Index rootIdx = Indexers::GetHBIndex<LEVELS_PER_PACK>(L_, pos, rootDepth);
      Index innerIdx = Indexers::GetLBIndex<LEVELS_PER_PACK>(L_, pos, depth);
      Index arr_idx = Indexers::GetArrIndex(L_, pos, depth);
     
      if (depth < DIRECTLY_CACHED_LEVELS) {   
        return data[arr_idx];
      } else {
        if (!cache.CheckContains(arr_idx)) {
          if (cache.IsFull()) {
            // UNDONE(): put correct typename here:
            uint64_t evictedIndex;
            {
              Bucket& evicted = cache.GetNextToEvict(evictedIndex).val;
              typename Bucket::Encrypted_t evictedEnc;
              evictedEnc.Encrypt(evicted);
              Index rootIdx, innerIdx;
              _ORAM::Indexers::GetBIndexFromArrIndex<LEVELS_PER_PACK>(L_, evictedIndex, rootIdx, innerIdx);
              LargeBucket_t& lb = server.Access(rootIdx);
              lb.buckets[innerIdx] = evictedEnc;
            }
            cache.EvictLRU(evictedIndex);
          }
          Index rootDepth = depth - (depth % LEVELS_PER_PACK);
          Index rootIdx = Indexers::GetHBIndex<LEVELS_PER_PACK>(L_, pos, rootDepth);
          Index innerIdx = Indexers::GetLBIndex<LEVELS_PER_PACK>(L_, pos, depth);
          LargeBucket_t& lb = server.Access(rootIdx);
          typename Bucket::Encrypted_t& eb = lb.buckets[innerIdx];
          Bucket b;
          eb.Decrypt(b);
          cache.Insert(arr_idx, b);
        }
        return cache.Access(arr_idx);
      }
    }

    void ReadBucketMetadata(const Index& pos, const Index& depth, BucketMetadata_t& ret) {
      const Bucket &bucket = GetBucketRef(pos, depth);
      ret = bucket.md;

      #ifdef OCI_ASSERTIONS
        Position realPos = Indexers::GetArrIndex(L_, pos, depth);

        // TRACE_FUNCTION(pos, realPos, depth, ret, &bucket);
      
        Duplet location = Duplet{realPos, depth};
        if (checkdata_metadata.contains(location)) {
          BucketMetadata_t& expected = checkdata_metadata[location]; 
          Assert(expected == ret, "\n", expected, "\n", ret, pos, realPos, depth, ret, &bucket);
        } else {
          Assert(false);
        }
      #endif
    }
    
    void WriteBucketMetadata(const Index& pos, const Index& depth, const BucketMetadata_t& val) {
      Bucket &bucket = GetBucketRef(pos, depth);
      bucket.md = val;

      #ifdef OCI_ASSERTIONS
        Position realPos = Indexers::GetArrIndex(L_, pos, depth);
        // TRACE_FUNCTION(pos, realPos, depth, val, &bucket);

        Duplet location = Duplet{realPos, depth};
        checkdata_metadata[location] = val;
      #endif
    }

    void ReadBlock(const Index& pos, const Index& depth, const Index& offset, Block_t& ret) {
      Assert(offset < Bucket::BUCKET_SIZE);
      const Bucket &bucket = GetBucketRef(pos, depth);
      ret = bucket.blocks[offset];

      #ifdef OCI_ASSERTIONS
        Position realPos = Indexers::GetArrIndex(L_, pos, depth);

        // if (ret.data.data.k != uint64_t(-1)) {
        //   TRACE_FUNCTION(pos, realPos, depth, offset, ret);
        // }

        Triplet location = Triplet{realPos, depth, offset};
        if (checkdata_block.contains(location)) {
          auto& expected = checkdata_block[location];
          Assert(expected == ret, "\n", expected, "\n", ret);
        } else {
          Assert(false);
        }
      #endif
    }

    void WriteBlock(const Index& pos, const Index& depth, const Index& offset, const Block_t& val) {
      #ifdef OCI_ASSERTIONS
        Position realPos = Indexers::GetArrIndex(L_, pos, depth);
        // if (val.data.data.k != uint64_t(-1)) {
        //   TRACE_FUNCTION(pos, realPos, depth, offset, val);
        // }
      #endif
      Assert(offset < Bucket::BUCKET_SIZE);
      Bucket &bucket = GetBucketRef(pos, depth);
      bucket.blocks[offset] = val;
      #ifdef OCI_ASSERTIONS
        Triplet location = Triplet{realPos, depth, offset};
        checkdata_block[location] = val;
      #endif
    }
  };
  // </struct ORAMClientInterface>
} // namespace _ORAM::ORAMClientInterface