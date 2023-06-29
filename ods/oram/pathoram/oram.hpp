#pragma once
#include <inttypes.h>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>
#include <cstring>
#include <memory>
#include <random>
#include <unordered_map>
#include <algorithm>

#include "common/defs.hpp"
#include "common/utils.hpp"
#include "common/encutils.hpp"
#include "common/tracing/tracer.hpp"
#include "oram/common/block.hpp"
#include "oram/common/bucket.hpp"
#include "oram/common/oram_client_interface.hpp"
#include "oram/pathoram/bucket.hpp"
#include "common/osort.hpp"


// This file implements constant trace PathORAM: https://eprint.iacr.org/2013/280.pdf
//

// Notation:
// N: number of blocks
// L: height of binary tree. Single node has L = 0
// B: block size in bits
// Z: number of blocks in bucket (4 in our case)
namespace _ORAM::PathORAM::ORAMClient {
template<typename T=Block::DefaultBlockData
  , unsigned int Z=ORAM__Z
  , bool ENCRYPT_BUCKETS=true
  , bool ENCRYPT_LARGE_BUCKETS=ORAM_SERVER__ENCRYPT_LARGE_BUCKETS
  , unsigned int LEVELS_PER_PACK=ORAM_SERVER__LEVELS_PER_PACK
  , unsigned int DIRECTLY_CACHED_LEVELS=ORAM_SERVER__DIRECTLY_CACHED_LEVELS
  , bool ObliviousCPUTrace=true>
struct ORAMClient {
  using _T = T;
  using Block_t = typename Block::Block<T, false>;
  using StashedBlock_t = typename StashedBlock::StashedBlock<Block_t>;
  using Bucket_t = typename _ORAM::PathORAM::Bucket::Bucket<Block_t,Z,ENCRYPT_BUCKETS>;
  using BucketMetadata_t = typename Bucket_t::BucketMetadata_t;
  using ORAMClientInterface_t = typename ORAMClientInterface::ORAMClientInterface<Block_t, Bucket_t, ENCRYPT_LARGE_BUCKETS, LEVELS_PER_PACK, DIRECTLY_CACHED_LEVELS>;
  
  ORAMClientInterface_t oramServerClient;

  struct State {
    bool forcedIntoStash;
    Position savedPath;
  };

  uint64_t N_;
  uint64_t L_;
  uint64_t S_;

  // Needs to be in oblivious memory:
  //
  std::vector<StashedBlock_t> stash_;
  State state;

  // We just keep a list of all the buckets we have to reorder,
  // so we can reuse them later on when we do a stash reshufle.
  //
  std::vector<std::tuple<Position, Index> > bucketsToReorder_; 

  bool noInit;

  // This constructor initializes ORAM with an empty AVL tree.
  //
  explicit ORAMClient(uint64_t N, bool _noInit=false) : 
      oramServerClient(N, _noInit)
    , N_(N)
    , L_(CeilLog2(N))
    , S_(Z * (CeilLog2(N)+1) + ORAM__MS)
    , state{false, DUMMY_POSITION}
    , noInit(_noInit)
    {
    PROFILE_F();
      
    // Fill the server with encryptions of dummy block
    // 
    if (!noInit) {
      Bucket_t bucket;
      BucketMetadata_t md;
      for (uint64_t l = 0; l <= L_; ++l) {
        for (uint64_t i = 0; i < (1 << l); ++i) {
          uint64_t pos = i << (L_ - l);
          for (uint64_t j = 0; j < Bucket_t::BUCKET_SIZE; ++j) {
            bucket.blocks[j] = Block_t::DUMMY();
            oramServerClient.WriteBlock(pos, l, j, bucket.blocks[j]);
            md = BucketMetadata_t::DUMMY();
          }
          oramServerClient.WriteBucketMetadata(pos, l, md);
        }
      }
    }
    X_LOG("ORAMClientInterface: Done setting client up!");

  }

  void BeginAccess(const ORAMAddress& oaddress, const Position& newPos, Block_t& ret, const bool markCached=true) {
    // std::cout << "BeginAccess: (" << oaddress.address << ", " << oaddress.position << ")" << std::endl; 
    // TRACE_FUNCTION(oaddress, newPos, markCached);
    PROFILE_F();
    Assert(oaddress.position < N_);
    Assert(oaddress.address < N_);
    Assert(newPos < N_);
    Assert(state.savedPath == DUMMY_POSITION);
    state.savedPath = oaddress.position;


    StashedBlock_t block;
    ReadPathToStash(oaddress.position);
    ReadBlockFromStash(oaddress, block, markCached, /*updatePosition=*/ true, newPos);

    ret = block.block;
  }

  void ReadPathToStash(const Position path) {
    PROFILE_F();
    Assert(path < N_);

    for (Index i=0; i <= L_; i++) {
      BucketMetadata_t md;
      oramServerClient.ReadBucketMetadata(path, i, md);
      ReadBucketToStash(path, i, md);
      oramServerClient.WriteBucketMetadata(path, i, md);
    }  
  }

  void ReadBucketToStash(const Position& pos, const Index& depth, BucketMetadata_t& md) {
    Assert(pos < N_);

    // First count the number of valid non dummies:
    //
    uint32_t validNonDummies = 0;
    for (Index j=0; j < Bucket_t::BUCKET_SIZE; j++) {
      StashedBlock_t sb = DUMMY<StashedBlock_t>();
      oramServerClient.ReadBlock(pos, depth, j, sb.block);
      sb.oaddress = md.priv.addresses[j];
      md.priv.addresses[j] = DUMMY<ORAMAddress>();
      stash_.push_back(sb);
    }
  }

  void WriteBucket(const Position& pos, const Index& depth) {
    PROFILE_F();
    StashedBlock_t blocks[Bucket_t::BUCKET_SIZE];
    BucketMetadata_t md = BucketMetadata_t::DUMMY();
    // DEBUGONLY:
    // oramServerClient.ReadBucketMetadata(pos, i, md);
    // md = GetBucketMetadata(pos, depth);
    // Assert(md.pub.invalidated);

    for (int j=0; j<Bucket_t::BUCKET_SIZE; j++) {
      blocks[j] = StashedBlock_t::DUMMY();
    }

    if constexpr (ObliviousCPUTrace) {
      for (int j=0; j < ORAM__Z; j++) {
        bool hasEvicted = false;
        for (int i=0; i < stash_.size(); i++) {
          StashedBlock_t& b = stash_[i];
          bool shouldEvict = true
          * (!b.cached)
          * (Indexers::PathsIntercept<LEVELS_PER_PACK>(L_, b.oaddress.position, pos, depth))
          * (!hasEvicted)
          ;

          CXCHG(shouldEvict, blocks[j], stash_[i]);
          CMOV(shouldEvict, stash_[i].cached, false);
          CMOV(shouldEvict, hasEvicted, true);
        }
      }
    }
    else 
    {
      // Non double oblivious version:
      //
      int j=0;
      for (int i=0; i < stash_.size(); i++) {
        if (j == ORAM__Z) break;
        StashedBlock_t& b = stash_[i];
        bool shouldEvict = true
        * (!b.cached)
        * (Indexers::PathsIntercept<LEVELS_PER_PACK>(L_, b.oaddress.position, pos, depth))
        ;

        if (shouldEvict) {
          CXCHG(shouldEvict, blocks[j], b);
          j++;
        }
      }
    }


    for (int j=0; j < Bucket_t::BUCKET_SIZE; j++) {
      md.priv.addresses[j] = blocks[j].oaddress;
      oramServerClient.WriteBlock(pos, depth, j, blocks[j].block);
    }
    
    oramServerClient.WriteBucketMetadata(pos, depth, md);
  }

  void ReadBlockFromStash(ORAMAddress oaddress, StashedBlock_t& ret, bool markCached=false, bool updatePosition=false, Position newPos=DUMMY_POSITION) {
    PROFILE_F();
    bool found = false;
    for (uint64_t i=0; i<stash_.size(); i++) {
      StashedBlock_t& stashi = stash_[i];
      bool isRightBlock = stashi.oaddress.address == oaddress.address;
      CMOV(isRightBlock, found, true);
      CMOV(isRightBlock*updatePosition, stashi.oaddress.position, newPos);
      CMOV(isRightBlock*markCached, stashi.cached, true);
      CMOV(isRightBlock, ret, stashi);
    }
    Assert(found, oaddress, "\n", ret, "\n", updatePosition, "\n", newPos);
  }

  void WriteToSomeStashedBlock(const Address& address, const Block_t& block, const bool& markUncached=false) {
    PROFILE_F();
    bool found = false;
    for (uint64_t i=0; i<stash_.size(); i++) {
      StashedBlock_t& stashi = stash_[i];
      bool isRightBlock = stashi.oaddress.address == address;
      CMOV(isRightBlock, found, true);
      CMOV(isRightBlock, stashi.block, block);
      CMOV(isRightBlock*markUncached, stashi.cached, false);
    }
    Assert(found);
  }

  void FinishAccess() {
    PROFILE_F();
    Assert(state.forcedIntoStash || state.savedPath != DUMMY_POSITION);
    if (state.savedPath != DUMMY_POSITION) {
      if constexpr (ObliviousCPUTrace) {
        EvictPath();
      } else {
        EvictPath_ExternalOblivious();
        CleanUpStash_ExternalOblivious();
      }
      if (stash_.size() > ORAM__MS) {
        stash_.resize(ORAM__MS);
      }
    }
    state.forcedIntoStash = false;
  }

  void UncacheAll() {
    for (uint64_t i=0; i<stash_.size(); i++) {
      stash_[i].cached = false;
    }
  }

  void ForceIntoStash(const ORAMAddress oaddr, const bool cached=true) {
    PROFILE_F();
    Assert(oaddr.position < N_);
    Assert(oaddr.address < N_);
    stash_.push_back(StashedBlock_t{cached, oaddr, Block_t::DUMMY()});
    state.forcedIntoStash = true;
  }

  void PurgeFromStash(const Address& address, const bool enable=true) {
    PROFILE_F();
    #ifndef NDEBUG
    bool purged = false;
    #endif


    StashedBlock_t SDUMMY = StashedBlock_t::DUMMY();
    // If enable, purge block from stash:
    //
    for (auto i=0; i < stash_.size(); i++) {
      StashedBlock_t& stashi = stash_[i];

      bool shouldPurge = true
        * (enable)
        * (stashi.oaddress.address == address);
      
      bool shouldClear = true
        * (!enable)
        * (stashi.oaddress.address == address);
      
      #ifndef NDEBUG
        Assert(!shouldPurge || !purged);
        CMOV(shouldPurge, purged, true);
      #endif
      
      CXCHG(shouldPurge, stashi, SDUMMY);
      CMOV(shouldClear, stashi.cached, false);
    }

    #ifndef NDEBUG
    Assert(enable == purged);
    #endif
  }

  void _DumpStash() {
    // for (uint64_t i=0; i<stash_.size(); i++) {
    //   if (stash_[i].oaddress != ORAMAddress::DUMMY()) {
    //     X_LOG_SIMPLE(stash_[i].oaddress, ", ");
    //   }
    // }
    // X_LOG_SIMPLE("\n");
  }

  void _DumpVector(std::vector<TaggedT<StashedBlock_t>>& v, uint64_t x) {
    // for (uint64_t i=0; i<v.size(); i++) {
    //   if (v[i].v.oaddress != ORAMAddress::DUMMY() || !v[i].isDummy || x == 3) {
    //     X_LOG_SIMPLE(i, ": ", v[i].v.oaddress, "[", v[i].tag, (v[i].isDummy ? "D" : ""), "], ");
    //   }
    // }
    // X_LOG_SIMPLE("\n");
  }

private:
  void EvictPath() {
    PROFILE_F();
    Assert(state.savedPath != DUMMY_POSITION);
    PROFILE_V(T_STASH_SIZE, stash_.size());

    std::vector<TaggedT<StashedBlock_t>> v;
    

    // 1) Add extra buckets for padding:
    //
    {
      for (uint64_t i=0; i<stash_.size(); i++) {
        StashedBlock_t& si = stash_[i];
        v.push_back(TaggedT<StashedBlock_t>{
          si.oaddress.address == DUMMY_ADDRESS,
          Indexers::GetDeepestShallowness(L_, si.oaddress.position, state.savedPath),
          si
        });
      }
      for (uint64_t i=0; i<=L_; i++) {
        v.push_back(TaggedT<StashedBlock_t>{
          true,
          i,
          StashedBlock_t::DUMMY()
        });
      }
    }

    uint64_t origSize = v.size();

    // 2) Sort with dummy buckets at the end of each level:
    //
    {
      auto cmp = [](const TaggedT<StashedBlock_t>& b1, const TaggedT<StashedBlock_t>& b2) {
        return 
            (b1.tag < b2.tag)
          + (  b1.tag == b2.tag
              * !b1.isDummy
              * b2.isDummy
            );
      };
      ObliSortBlocks<TaggedT<StashedBlock_t>>(v, cmp);
    }

    // 3) Do a linear pass to count buckets:
    //
    uint64_t count = 0;
    uint64_t currTag = 0;
    {      
      for (uint64_t i=0; i<origSize; i++) {
        bool isEnd = v[i].isDummy                  
                  * (currTag <= L_);
        Assert(!isEnd || (v[i].tag == currTag));      
        CMOV(isEnd, v[i].v.oaddress.address, static_cast<_ORAM::Address>(count));
        CMOV(isEnd, currTag, currTag+1);
        count += 1;
        CMOV(isEnd, count, (uint64_t)0);
      }
      Assert(currTag == L_+1);
    }

    // 4) Do a second sort to make v[0..L] be the counter buckets:
    {
      auto cmp = [L_=L_](const TaggedT<StashedBlock_t>& b1, const TaggedT<StashedBlock_t>& b2) {
        // I swear this is a strict weak ordering.
        //
        return false
          + ((b1.isDummy == b2.isDummy) * (b1.tag < b2.tag))
          + ((b1.isDummy) * (!b2.isDummy) * (b1.tag < L_+1))
          + ((!b1.isDummy) * (b2.isDummy) * (b2.tag >= L_+1));         
      };
      ObliSortBlocks<TaggedT<StashedBlock_t>>(v, cmp);
    }


    count = 0;
    // 5) Do a linear pass to add remaining dummy buckets:
    //
    for (uint64_t i=0; i<=L_; i++) {
      Assert(v[i].isDummy * (v[i].tag == i));
      const uint64_t target = (i+1)*Z;
      count += v[i].v.oaddress.address;
      
      // Turn this block into a dummy block directly:
      //
      v[i].v.oaddress.address = DUMMY_ADDRESS;
      CMOV(count >= target, v[i].tag, L_+1);
      CMOV(count < target, count, count+1);
      for (uint64_t z=0; z<Z-1; z++) {
        bool rem = count < target;
        CMOV(rem, count, count+1);
        {
          TaggedT<StashedBlock_t> newDummy{
            false,
            i,
            StashedBlock_t::DUMMY()
          };
          CMOV(!rem, newDummy.tag, L_+1);
          v.push_back(newDummy);
        }       
      }
    }

    // 6) Sort to put all the buckets in the right position for eviction
    //
    {
      // The second case in the lambda is not used, but we still leave it there:
      //
      auto cmp = [](const TaggedT<StashedBlock_t>& b1, const TaggedT<StashedBlock_t>& b2) {
        return 
            (b1.tag < b2.tag)
          + (  b1.tag == b2.tag
              * !b1.isDummy
              * b2.isDummy
            );
      };
      ObliSortBlocks<TaggedT<StashedBlock_t>>(v, cmp);
    }

    // 7) Evict directly, and then move all the remaining to the stash, padding 
    // to the pad size:
    //
    for (uint64_t i=0; i<=L_; i++) {
      BucketMetadata_t md = BucketMetadata_t::DUMMY();
      StashedBlock_t blocks[Z];
      for (uint64_t j=0; j<Z; j++) {
        md.priv.addresses[j] = v[i*Z+j].v.oaddress;
        oramServerClient.WriteBlock(state.savedPath, L_-i, j, v[i*Z+j].v.block);
      }
      oramServerClient.WriteBucketMetadata(state.savedPath, L_-i, md);
    }
    
    state.savedPath = DUMMY_POSITION;
    stash_.clear();
    for (uint64_t i=Z*(L_+1); i<=std::min(Z*(L_+1) + S_,v.size()); i++) {
      // This if leaks the number of pad blocks, which is public anyways.
      // UNDONE(): put the right number on the loop bound so the if is not needed
      //
      if (!v[i].isDummy) { 
        stash_.push_back(v[i].v);
      } else {
        break;
      }
    }

    PROFILE_V(T_FINAL_STASH_SIZE_1, stash_.size());
    
    uint64_t cachedCount = 0;
    for (uint64_t i=0; i<stash_.size(); i++) {
      CMOV(stash_[i].cached, cachedCount, cachedCount+1);
    }
    PROFILE_V(T_FINAL_STASH_SIZE_2, cachedCount);

  }


  void EvictPath_ExternalOblivious() {
    // PROFILE_F();
    Assert(state.savedPath != DUMMY_POSITION);
    PROFILE_V(T_STASH_SIZE, stash_.size());
    
    // UNDONE(): fast version of this: (lambda log lambda version)
    //
    for (int64_t i=L_; (i+1!=0); i--) {
      BucketMetadata_t md;
      WriteBucket(state.savedPath,i);
    }

    state.savedPath = DUMMY_POSITION;
    CleanUpStash_ExternalOblivious();
    PROFILE_V(T_FINAL_STASH_SIZE_1, stash_.size());
    
    uint64_t cachedCount = 0;
    for (uint64_t i=0; i<stash_.size(); i++) {
      CMOV(stash_[i].cached, cachedCount, cachedCount+1);
    }
    PROFILE_V(T_FINAL_STASH_SIZE_2, cachedCount);
  }


  void CleanUpStash_ExternalOblivious() {
    // PROFILE_F();
    auto cmp = [](const StashedBlock_t& b1, const StashedBlock_t& b2) {
      return 
          (b1.oaddress.address != DUMMY_ADDRESS)
        * (b2.oaddress.address == DUMMY_ADDRESS);
    };
    
    int ctr=0;
    int dummies = 0;
    for (int i=0; i<stash_.size()-dummies; i++) {
      if (stash_[i].oaddress.address != DUMMY_ADDRESS) {        
        ctr++;
      } else {
        dummies++;
        if (i < stash_.size() - dummies) {
          std::swap(stash_[i], stash_[stash_.size()-dummies]);
          i--;
        }
      }
    }
    // ObliSortBlocks<StashedBlock_t>(stash_, cmp);
    auto newSize = std::min(stash_.size(), S_);
    // UNDONE(): recursive ORAM is needing the following, debug why, should not be needed.
    //
    newSize = ctr;
    if (newSize < ctr) {
      Assert(false, "Weak security parameter", ctr, stash_.size());
    }
    if (newSize < stash_.size()) {
      Assert(stash_[newSize].oaddress.address == DUMMY_ADDRESS);
      stash_.resize(newSize);
    }
  }


}; // </struct ORAMClient>
} // namespace _ORAM::PathORAM::ORAMClient
