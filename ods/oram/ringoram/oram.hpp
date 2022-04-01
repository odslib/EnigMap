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
#include "oram/ringoram/bucket.hpp"

// This file implements constant trace RingORAM: https://eprint.iacr.org/2014/997.pdf
//

// Notation:
// N: number of blocks
// L: height of binary tree. Single node has L = 0
// B: block size in bits
// Z: number of blocks in bucket (4 in our case)
namespace _ORAM::RingORAM::ORAMClient {
template<typename T=Block::DefaultData_t
  , unsigned int Z=ORAM__Z
  , unsigned int S=ORAM__S
  , bool ENCRYPT_BLOCKS=ORAM__ENCRYPT_BLOCKS
  , bool ENCRYPT_LARGE_BUCKETS=ORAM_SERVER__ENCRYPT_LARGE_BUCKETS
  , unsigned int LEVELS_PER_PACK=ORAM_SERVER__LEVELS_PER_PACK
  , unsigned int DIRECTLY_CACHED_LEVELS=ORAM_SERVER__DIRECTLY_CACHED_LEVELS
  , bool ObliviousCPUTrace=false>
struct ORAMClient {
  using _T = T;
  using Block_t = typename Block::Block<T, ENCRYPT_BLOCKS>;
  using StashedBlock_t = typename StashedBlock::StashedBlock<Block_t>;
  using Bucket_t = typename _ORAM::RingORAM::Bucket::Bucket<Block_t,Z,S,ENCRYPT_BLOCKS>;
  using BucketMetadata_t = typename Bucket_t::BucketMetadata_t;
  using ORAMClientInterface_t = typename ORAMClientInterface::ORAMClientInterface<Block_t, Bucket_t, ENCRYPT_LARGE_BUCKETS, LEVELS_PER_PACK, DIRECTLY_CACHED_LEVELS>;

  ORAMClientInterface_t oramServerClient;

  uint64_t round = 0;
  uint64_t G = 0;

  struct State {
    Position savedPath;
  };

  uint64_t N_;
  uint64_t L_;

  // Needs to be in oblivious memory:
  //
  std::vector<StashedBlock_t> stash_;

  // We just keep a list of all the buckets we have to reorder,
  // so we can reuse them later on when we do a stash reshufle.
  //
  std::vector<std::tuple<Position, Index> > bucketsToReorder_; 

  // This constructor initializes ORAM with empty data.
  //
  explicit ORAMClient(uint64_t N) : 
      oramServerClient(N)
    , N_(N)
    , L_(CeilLog2(N))
    {
    // Fill the server with encryptions of dummy block
    // 
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
    X_LOG("ORAMClientInterface: Done setting client up!");
  }

  // This is done in a different order from the Access algorithm
  // because it's easier to make oblivious if we go trough the stash first
  // and then just mark we want to access dummy nodes.
  //
  void BeginAccess(const ORAMAddress& oaddress, const Position& newPos, Block_t& ret, const bool markCached=true) {
    // std::cout << "BeginAccess: (" << oaddress.address << ", " << oaddress.position << ")" << std::endl; 
    // TRACE_FUNCTION(oaddress, newPos, markCached);
    PROFILE_F();
    Assert(oaddress.position < N_);
    Assert(oaddress.address < N_);
    Assert(newPos < N_);

    StashedBlock_t block;
    ReadBlockOnPathToStash(oaddress);
    ReadBlockFromStash(oaddress, block, markCached, true, newPos);

    ret = block.block;
  }

  void ReadBlockToStashedBlock(ORAMAddress oaddress, const Index& depth, BucketMetadata_t& md, StashedBlock_t& out) {
    PROFILE_F();
    if (md.pub.invalidated) return;
    if (md.pub.counter >= ORAM__S) {
      // Then we just evict the whole bucket to ORAM,
      // this is what early reshuffle does, but doing it later
      // results in smaller stash size for free.
      //
      ReadBucketToStash(oaddress.position, depth, md);
    } else {
      bool hasBlock = false;
      uint32_t evictIndex = -1;
      uint32_t freeDummies = 0;
      for (uint32_t j=0; j < Bucket_t::BUCKET_SIZE; j++) {
        bool match = true
          * (md.priv.addresses[j] == oaddress)
          * (md.pub.valid[j]);
        
        bool isValidDummy = true
          * (md.priv.addresses[j].address == DUMMY_ADDRESS)
          * (md.pub.valid[j]);
        
        Assert(!match || !hasBlock);
        CMOV(match, evictIndex, j);
        CMOV(match, hasBlock, match);
        CMOV(isValidDummy, freeDummies, freeDummies+1);
      }

      // Chose index to evict:
      //
      uint32_t randomDummyIndex = random() % freeDummies;
      
      for (uint32_t j=0; j < Bucket_t::BUCKET_SIZE; j++) {
        bool isValidDummy = true
          * (md.priv.addresses[j].address == DUMMY_ADDRESS)
          * (md.pub.valid[j]);
        
        bool isNewEvictIndex = true
         * (isValidDummy)
         * (!hasBlock)
         * (randomDummyIndex == 0);
        
        CMOV(isValidDummy, randomDummyIndex, randomDummyIndex-1);
        CMOV(isNewEvictIndex, evictIndex, j);
      }

      Assert(evictIndex <= Bucket_t::BUCKET_SIZE);
      md.pub.valid[evictIndex] = false;
      md.pub.counter += 1;
      
      StashedBlock_t sb;
      oramServerClient.ReadBlock(oaddress.position, depth, evictIndex, sb.block);
      sb.cached = false;
      sb.oaddress = md.priv.addresses[evictIndex];
      // NOTE: Make sure than when we reduce the stash size we don't leak it.
      //

      bool shouldMove = sb.oaddress.address == oaddress.address;
      CXCHG(shouldMove, out, sb);
    }
  }

  void ReadBlockOnPathToStash(ORAMAddress oaddress) {
    PROFILE_F();
    StashedBlock_t newBlock = StashedBlock_t::DUMMY();

    for (Index i=0; i <= L_; i++) {
      BucketMetadata_t md;
      oramServerClient.ReadBucketMetadata(oaddress.position, i, md);
      if (md.pub.invalidated) continue;
      ReadBlockToStashedBlock(oaddress, i, md, newBlock);
      oramServerClient.WriteBucketMetadata(oaddress.position, i, md);
    }
    
    stash_.push_back(newBlock);
  }

  void ReadBucketToStash(const Position& pos, const Index& depth, BucketMetadata_t& md) {
    // UNDONE(): rewrite this to actually follow the logic of only evicting the maximum number of
    // valid non dummies.
    //
    if (md.pub.invalidated) return;
    md.pub.invalidated = true;
    bucketsToReorder_.push_back({pos, depth});

    // First count the number of valid non dummies:
    //
    uint32_t validNonDummies = 0;
    for (Index j=0; j < Bucket_t::BUCKET_SIZE; j++) {
      bool isValid = md.pub.valid[j];
      
      bool isValidDummy = true
        * (isValid)
        * (md.priv.addresses[j].address == DUMMY_ADDRESS);
      
      CMOV(isValid * (!isValidDummy), validNonDummies, validNonDummies+1);
    }
    uint32_t dummiesToEvict = ORAM__Z - validNonDummies;

    for (Index j=0; j < Bucket_t::BUCKET_SIZE; j++) {
      bool isValid = md.pub.valid[j];
      
      bool isValidDummy = true
        * (isValid)
        * (md.priv.addresses[j].address == DUMMY_ADDRESS);
        
      bool shouldEvict = true
        * (isValid)
        * ( (  (isValidDummy)
            * (dummiesToEvict > 0)
            )
          + (!isValidDummy));

      if (shouldEvict) {
        StashedBlock_t sb = StashedBlock_t::DUMMY();
        oramServerClient.ReadBlock(pos, depth, j, sb.block);
        sb.cached = false;
        // This might be a dummy block, but we need to stash it anyways.
        //
        sb.oaddress = md.priv.addresses[j]; 
        // NOTE: Make sure than when we reduce the stash size we don't leak how many dummies are there.
        //
        stash_.push_back(sb);
        CMOV(isValidDummy, dummiesToEvict, dummiesToEvict-1);
      } else {
        if (isValidDummy) {
          StashedBlock_t sb = StashedBlock_t::DUMMY();
          oramServerClient.ReadBlock(pos, depth, j, sb.block);
        }
      }
    }
  }

  // UNDONE(): simplify given new bucket interface
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
          CXCHG(shouldEvict, blocks[j], stash_[i]);
          CMOV(shouldEvict, stash_[i].cached, false);
          j++;
        }
      }
    }

    // UNDONE(): do a random permutation of the blocks.
    //

    for (int j=0; j < Bucket_t::BUCKET_SIZE; j++) {
      md.pub.valid[j] = true;
      md.priv.addresses[j] = blocks[j].oaddress;
      oramServerClient.WriteBlock(pos, depth, j, blocks[j].block);
    }
    md.pub.counter = 0;
    
    oramServerClient.WriteBucketMetadata(pos, depth, md);
  }

  void ReadBlockFromStash(ORAMAddress oaddress, StashedBlock_t& ret, bool markCached=false, bool updatePosition=false, Position newPos=DUMMY_POSITION) {
    PROFILE_F();
    bool found = false;
    for (uint64_t i=0; i<stash_.size(); i++) {
      bool isRightBlock = stash_[i].oaddress.address == oaddress.address;
      if (isRightBlock) {
        found = found + true;
        if (updatePosition) {
          stash_[i].oaddress.position = newPos;
        }
        if (markCached) {
          stash_[i].cached = true;
        }
        ret = stash_[i];
      }
    }
    Assert(found, oaddress, "\n", ret, "\n", updatePosition, "\n", newPos);
  }

  void WriteToSomeStashedBlock(const Address& address, const Block_t& block, const bool& markUncached=false) {
    PROFILE_F();
    bool found = false;
    for (uint64_t i=0; i<stash_.size(); i++) {
      bool isRightBlock = stash_[i].oaddress.address == address;
      if (isRightBlock) {
        found = found + true;
        stash_[i].block = block;
        if (markUncached) {
          stash_[i].cached = false;
        }
      }
    }
    Assert(found);
  }

  void FinishAccess() {
    PROFILE_F();
    round = (round + 1) % ORAM__A;
    if (round == 0) {
      EvictPath();
    }
    // We modified the algorithms to include early reshufle in the 
    // readpath/readblock logic, as it reduces the required stash size 
    // and avoids having to cache metadata.
    //
    // EarlyReshuffle();
    //
    LateReshuffle();
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
      bool shouldPurge = true
        * (enable)
        * (stash_[i].oaddress.address == address);
      
      #ifndef NDEBUG
        Assert(!shouldPurge || !purged);
        purged = purged + shouldPurge;
      #endif
      
      CXCHG(shouldPurge, stash_[i], SDUMMY);
    }

    #ifndef NDEBUG
    Assert(enable == purged);
    #endif
  }

private:
  void EvictPath() {
    PROFILE_F();
    Position pos = ToReverseLexicographicalOrder(G, L_);
    G = (G+1)%N_;
    
    for (int64_t i=0; i<L_; i++) {
      BucketMetadata_t md;
      oramServerClient.ReadBucketMetadata(pos, i, md);
      ReadBucketToStash(pos,i,md);
      oramServerClient.WriteBucketMetadata(pos, i, md);
    }

    // LateReshufle above will remove the buckets from the stash.
    //
  }


  // We are doing the late reshufle heuristic, which does the same
  // as early reshufle, but keeps things on the stash until the buckets
  // are actually accessed, and allows us not to have to cache the metadata
  // decryptions.
  //
  // void EarlyReshuffle() {
  //   for (Index i=0; i<L_; i++) {
  //     ReadBucketToStash(l,i);
  //     WriteBucket(l,i);
  //   }
  // }

  void LateReshuffle() {
    PROFILE_F();
    // UNDONE(): only do this if the public stash size is above some constant?
    // (this should improve even further the maximum stash size, as we have a higher probability of
    // evicting to lower levels).
    //

    // First sort bucketsToReorder from the deepest to the least deep.
    //
    std::sort(bucketsToReorder_.begin(), bucketsToReorder_.end());
    std::reverse(bucketsToReorder_.begin(), bucketsToReorder_.end());

    PROFILE_V(T_STASH_SIZE, stash_.size());
    PROFILE_V(T_INIT_STASH_SIZE_1, stash_.size());
    CleanUpStash();
    PROFILE_V(T_STASH_SIZE, stash_.size());
    PROFILE_V(T_FINAL_STASH_SIZE_1, stash_.size());
    
    // For each bucket to reorder, we call writebucket.
    //
    for (const auto& [pos, depth] : bucketsToReorder_) {
      WriteBucket(pos, depth);
    }
    
    bucketsToReorder_.clear();

    PROFILE_V(T_STASH_SIZE, stash_.size());
    PROFILE_V(T_INIT_STASH_SIZE_2, stash_.size());
    CleanUpStash();
    PROFILE_V(T_STASH_SIZE, stash_.size());
    PROFILE_V(T_FINAL_STASH_SIZE_2, stash_.size());
  }

  void CleanUpStash() {
    PROFILE_F();
    auto cmp = [](const StashedBlock_t& b1, const StashedBlock_t& b2) {
      return 
          (b1.oaddress.address != DUMMY_ADDRESS)
        * (b2.oaddress.address == DUMMY_ADDRESS);
    };
    // UNDONE(): do a random permutation of the stash
    //


    int ctr=0;
    int dummies = 0;
    for (int i=0; i<stash_.size()-dummies; i++) {
      if (stash_[i].oaddress.address != DUMMY_ADDRESS) {
        ctr++;
      } else {
        dummies++;
        std::swap(stash_[i], stash_[stash_.size()-dummies]);
        i--;
      }
    }
    auto newSize = std::min(stash_.size(), (size_t)ORAM__MS);
    newSize = ctr; // For debugging.
    // ObliSortBlocks(stash_, cmp);
    Assert(stash_.size() == newSize || stash_[newSize].oaddress.address == DUMMY_ADDRESS);
    stash_.resize(newSize);
    PROFILE_V(T_STASH_SIZE, newSize);
  }


}; // </struct ORAMClient>
} // namespace _ORAM::RingORAM::ORAMClient
