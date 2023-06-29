#pragma once
#include "oram/common/types.hpp"
#include "oram/common/block.hpp"
#include "external_memory/emvector.hpp"
#include "external_memory/algorithm/sort.hpp"

namespace _OBST {
namespace OramClient {

using namespace _ORAM;

template <typename T>
concept HasBlockOfNode = std::is_same_v<typename T::_T, Node>;
template<typename ORAMClient /*=_ORAM::ORAMClient::ORAMClient<Node,ORAM__Z,ORAM__S,false,false,4>*/ >
requires HasBlockOfNode<ORAMClient>
struct OramClient {
  using _T = typename ORAMClient::_T;
  using Block_t = typename ORAMClient::Block_t;

  ORAMClient oram;
  ORAMAddress nextFreeBlock;

  #ifndef NDEBUG
  struct State {
    // This structure is used to debug calls to this class follow a valid pattern.
    //
    enum Tp {
      READY,
      PREWRITE
    };

    Tp tp;
    Address addr;
    Position pos;
    
    static State ReadyState() {
      return State{READY, 0, 0};
    }

    static State PreWriteState(Address _addr, Position _pos) {
      return State{PREWRITE, _addr, _pos};
    }
  };
  #endif

  Address N_;
  uint64_t L_;

  #ifndef NDEBUG
  State state;
  #endif

  explicit OramClient(uint64_t N, bool noInit=false) : N_(N), oram(N, noInit) {
    L_ = CeilLog2(N_);

    #ifndef NDEBUG 
    state = State::ReadyState();
    #endif

    nextFreeBlock = {0,0};
  }

  Position GenRandomPosition() {
    // UNDONE(0): This needs to be random
    //
    static Position pos = 0;
    pos = (pos + 1)%N_;
    return pos;
  }

  void BeginReadOp(const Address& i, const Position& oldPos, Block_t& ret, const Position& newPos) {
    PROFILE_F();

    #ifndef NDEBUG
    Assert(state.tp == State::Tp::READY);
    // 
    //
    // std::cout << "R:" << i << " " << oldPos << " ";
    Assert(oldPos < N_);
    Assert(newPos < N_);
    state = State::PreWriteState(i, newPos);
    #endif

    oram.BeginAccess(ORAMAddress{i,oldPos}, newPos, ret, true);
  }

  void FinishOramOp(const Address& address, const Block_t& toWrite, bool markUncached = true) {
    PROFILE_F();
    #ifndef NDEBUG
    const Address i = state.addr;
    Assert(state.tp == State::Tp::PREWRITE);
    Assert(state.addr == address);
    #endif

    // std::cout << "E:" << address 
      // << " (" << toWrite.data.child[0].address << ", " << toWrite.data.child[0].position << ") " 
      // << " (" << toWrite.data.child[1].address << ", " << toWrite.data.child[1].position << ") " 
      // << std::endl;
    oram.WriteToSomeStashedBlock(address, toWrite, markUncached);
    oram.FinishAccess();

    #ifndef NDEBUG
    state = State::ReadyState();
    #endif
  }

  void UncacheAll()
  {
    oram.UncacheAll();
  }

  // UNDONE(): Currently we do not support freeing arbitrary blocks, we want to 
  // add freeing of arbitrary blocks in the future.
  // UNDONE(): Currently we require the blocks to be cached on the stash to be freable.
  //
  // UNDONE(): Remove ret from here
  //
  INLINE void GetNextFreeBlockOp(Address& i, Block_t& ret, const Position& newPos) {
    // Allocates a new block that will reside in newPos.
    // This function, if fully implemented should be very similar to BeginOramOp, except it reads the
    // address and position of the next block on the freelist from the block
    // it just allocated.
    //
    #ifndef NDEBUG
    Assert(state.tp == State::Tp::READY);
    Assert(nextFreeBlock.address < N_);
    Assert(0 <= newPos < N_);
    #endif

    i = nextFreeBlock.address;

    #ifndef NDEBUG
    state = State::PreWriteState(i, newPos);
    #endif
    
    oram.ForceIntoStash(ORAMAddress{i, newPos}, /*cached=*/true);

    ret = Block_t::DUMMY();
    nextFreeBlock.address += 1;
  }

  // UNDONE(): remove currPos from here
  //
  INLINE void FreeBlockMaybeFOp(const Address& address, Position currPos, bool shouldFree) {
    #ifndef NDEBUG
    Assert(state.tp == State::Tp::READY);
    Assert(!shouldFree || ((address+1) == nextFreeBlock.address));
    // state = State::PreWriteState(DUMMY_ADDRESS, DUMMY_POSITION);
    Assert(0 <= currPos < N_);
    #endif
    
    oram.PurgeFromStash(address, /*enable=*/shouldFree);

    CMOV(shouldFree, nextFreeBlock.address, nextFreeBlock.address-1);
  }
/**/
  INLINE void InitViaLargeStashVector(EM::Vector::Vector<_ORAM::StashedBlock::StashedBlock<Block_t>> builtstash) {
    // UNDONE(24): Review and write tests for this code, it was rushed to have
    // experimental numbers for the optimal initialization algorithm, it probably has bugs.
    //
    EM::Vector::Vector<TaggedT<Block_t>> maxStash(builtstash.size(), TaggedT<Block_t>::DUMMY());

    for (uint64_t i=0; i<builtstash.size(); i++) {
      maxStash[i].isDummy = false;
      maxStash[i].v = builtstash[i].block;
      maxStash[i].tag = builtstash[i].oaddress.position;
    }

    uint64_t currSize = builtstash.size();


    // EM::Algorithm::BucketObliviousShuffle(maxStash.begin(), maxStash.end());

    for (uint64_t l=L_+1; l-->0;) {
      // Bellow it's the actuall algorithm, we are just doing K iterations of 
      // random permutations to simulate the largest overhead part of it.
      
      // EM::Algorithm::BucketObliviousShuffle(maxStash);
      // if (maxStash.N > 20) {
        // maxStash.N = GetNextPowerOfTwo(((maxStash.N/4)+1)*2);
      // }

      // auto cmp = [](const TaggedT<Block_t>& b1, const TaggedT<Block_t>& b2) {
      //   return 
      //       (b1.tag < b2.tag)
      //     + ( (b1.tag == b2.tag)
      //       * (b1.isDummy < b2.isDummy)
      //       );
      // };
      // ObliSortBlocks(maxStash.mem, cmp);
      // ObliSortBlocks(maxStash.mem, cmp);


      // 
      currSize /= 2;

      if (currSize == 0) {
        break;
      }
    }
    // Assert(maxStash[0].isDummy);
  }
/**/ 

}; // </struct OramClient>
} // </namespace OramClient>
} // </namespace _OBST>