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
#include "oram/common/indexers.hpp"
#include "oram/common/block.hpp"
#include "oram/common/concepts.hpp"

// This file implements an oram interface that doesn't really implement ORAM
// In debug mode this asserts that a position based ORAM is working as expected.
// It is usefull to debug algorithms that use ORAM.
//

namespace _ORAM::NotORAM::ORAMClient
{
  template <typename T = Block::DefaultBlockData
    , bool ENCRYPT_BLOCKS = ORAM__ENCRYPT_BLOCKS
    , bool NOT_ORAM_ASSERTIONS = true>
  requires true
    && ::Concepts::Encryptable<Block::Block<T, ENCRYPT_BLOCKS> >
    && (IS_POD<T>())
  struct ORAMClient
  {
    using _T = T;
    using Block_t = typename Block::Block<T, ENCRYPT_BLOCKS>;
    using StashedBlock = typename StashedBlock::StashedBlock<Block_t>;

    uint8_t key_[16] = {0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41, 0x41};
    std::random_device rd_; // uses RDRND or /dev/urandom
    std::uniform_int_distribution<uint64_t> dist_;

    uint64_t N_;

    struct Slot
    {
      typename Block_t::Encrypted_t v;
      struct FilledAssertionData {
        Position position;
        bool cached;
      };
      struct EmptyAssertionData { };
      using AssertionData = std::conditional_t<NOT_ORAM_ASSERTIONS, FilledAssertionData, EmptyAssertionData>;

      AssertionData ad;
      
      // UNDONE(): This can't be made consteval?
      //
      static INLINE Slot DUMMY() { 
        Slot ret;
        ret.v.Encrypt(Block_t::DUMMY());
        
        if constexpr (NOT_ORAM_ASSERTIONS) {
          ret.ad.position = DUMMY_POSITION;
          ret.ad.cached = false;
        }

        return ret;
      }
    };


    struct FilledState {
      enum Tp {
        READY,
        PREWRITE
      };

      Tp tp;
      Address addr;
      Position pos;
      
      static FilledState ReadyState() {
        return FilledState{READY, 0, 0};
      }

      static FilledState PreWriteState(Address _addr, Position _pos) {
        return FilledState{PREWRITE, _addr, _pos};
      }
    };
    struct EmptyState { };

    using State = std::conditional_t<NOT_ORAM_ASSERTIONS, FilledState, EmptyState>; 
    State state;

    std::vector<Slot> data_;

    explicit ORAMClient(uint64_t N, bool noInit = false)
    : N_(N)
    {
      data_ = std::vector(N_, Slot::DUMMY());

      if constexpr (NOT_ORAM_ASSERTIONS) {
        state = State::ReadyState();
      }
    }

    // Reads a block at a given position
    //
    void BeginAccess(const ORAMAddress &oaddress, const Position newPos, Block_t &ret, const bool markCached = true)
    {
      if constexpr (NOT_ORAM_ASSERTIONS) { 
        Assert(state.tp == State::Tp::READY);
      }
      Assert(oaddress.address < N_);
      Assert(oaddress.position < N_);
      Assert(newPos < N_);
      if constexpr (NOT_ORAM_ASSERTIONS) {
        state = State::PreWriteState(oaddress.address, newPos);
      }

      Slot &slot = data_[oaddress.address];
      slot.v.Decrypt(ret);

      if constexpr (NOT_ORAM_ASSERTIONS) {
        Assert(slot.ad.cached || slot.ad.position == oaddress.position);
        Assert(slot.ad.cached <= markCached);
        slot.ad.position = newPos;
        if (markCached) slot.ad.cached = true;
      }
    }

    // Writes to a block that is on the stash
    //
    void WriteToSomeStashedBlock(const Address address, const Block_t &block, const bool markUncached = false)
    {
      Assert(address < N_);
      Slot &slot = data_[address];
      if constexpr (NOT_ORAM_ASSERTIONS) {
        Assert(slot.ad.cached);
      }

      slot.v.Encrypt(block);

      if constexpr (NOT_ORAM_ASSERTIONS) {
        if (markUncached)
        {
          slot.ad.cached = false;
        }
      }
    }

    // This is called after every oram operation.
    // At this point communication with oram server might happen.
    //
    void FinishAccess()
    {
      if constexpr (NOT_ORAM_ASSERTIONS) {
        Assert(state.tp == State::Tp::PREWRITE);
        state = State::ReadyState();
      }
    }

    void UncacheAll() {
      // UNDONE(): this
    }

    // Appends a new block to the stash.
    //
    void ForceIntoStash(const ORAMAddress oaddr, const bool cached=true) {
      if constexpr (NOT_ORAM_ASSERTIONS) {
        Assert(state.tp == State::Tp::READY);
        state = State::PreWriteState(oaddr.address, oaddr.position);
      }
      Assert(oaddr.address < N_);
      Assert(oaddr.position < N_);
      data_[oaddr.address].v.Encrypt(Block_t::DUMMY());

      if constexpr (NOT_ORAM_ASSERTIONS) {        
        data_[oaddr.address].ad.position = oaddr.position;
        data_[oaddr.address].ad.cached = cached;
      }
    } 

    // Deletes a block that is on the stash.
    //
    void PurgeFromStash(const Address address, const bool enable=true) {
      if constexpr (NOT_ORAM_ASSERTIONS) {
        Assert(state.tp == State::Tp::READY); 
      }

      if (enable) {
        Assert(address < N_);
      }

      if constexpr (NOT_ORAM_ASSERTIONS) {          
        if (enable) {
          data_[address].ad.cached = false;
        }
      }
    }
  }; /// struct ORAMClient
} // namespace _ORAM::NotORAM::ORAMClient
