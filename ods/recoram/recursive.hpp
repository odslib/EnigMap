#pragma once
#include <vector>

#include "common/utils.hpp"
#include "oram/common/block.hpp"

// This file implements a recursive oram that works for position based orams.
// Initialization has  T_access * N log N time and logN memory
// Access has T_Access * log N time.
//
namespace _ORAM::RecursiveORAM {
  struct Data {
    _ORAM::Position pos[2];
    static consteval INLINE Data DUMMY() {return Data{DUMMY_POSITION, DUMMY_POSITION}; }
    bool operator==(const Data &other) const = default;
    friend std::ostream &operator<<(std::ostream &o, const Data &x) { 
      o << "{" << x.pos[0] << ", " << x.pos[1] << "}";
      return o;
    }
  };

  template <typename T, typename BaseORAM, typename PositionsORAM>
  // UNDONE():
  // requires BaseORAM::T == T;
  // requires PositionsORAM::T == Data;
  struct RecursiveORAM {
    using IntermidiateBlock_t = typename PositionsORAM::Block_t;
    using Block_t = typename BaseORAM::Block_t;
    uint64_t L_; // The number of intermediate ORAM's
    uint64_t N_; // Number of positions on the last ORAM.
    std::vector<std::unique_ptr<PositionsORAM> > orams;
    BaseORAM lastORAM;

    explicit RecursiveORAM(uint64_t N) :
      N_(N)
    , L_(CeilLog2(N))
    , lastORAM(N)
    {
      // UNDONE(): Use a single ORAM?
      //
      // UNDONE(): fast recursive ORAM initialization
      //
      // Builds the intermediate ORAMs by first assigning a constant position to every node in every level and then accessing every path (so everything looks random) -> ~O(N log^2 N) initialization.
      //
      orams.reserve(L_);
      for (uint64_t l=0; l<L_; l++) {
        orams.emplace_back(new PositionsORAM(1ULL<<l));
      }

      // Assign position i to every node in every level
      //
      for (uint64_t l=0; l<L_; l++) {
        for (_ORAM::Address n=0; n<(1<<l); n++) {
          const _ORAM::ORAMAddress oa = _ORAM::ORAMAddress{n, n};
          const Position lp = (n<<1) | 0;
          const Position rp = (n<<1) | 1;
          // std::cout << l << " " << n << " " << lp << " " << rp << std::endl;
          // UNDONE(): Actually have a path to evict :)
          //
          orams[l]->ForceIntoStash(oa);
          orams[l]->WriteToSomeStashedBlock(n, IntermidiateBlock_t{Data{{lp, rp}}}, /*uncache=*/true);
          orams[l]->FinishAccess();
        }
      }

      for (_ORAM::Address n=0; n<N_; n++) {
        // UNDONE(): Actually have a path to evict :)
        //
        const _ORAM::ORAMAddress oa = _ORAM::ORAMAddress{n, n};
        lastORAM.ForceIntoStash(oa, /*cache=*/false);
        lastORAM.FinishAccess();
      }

      // Access every path:
      //
      for (_ORAM::Address n=0; n<N_; n++) {
        Block_t tmpBlock;
        BeginAccess(n, tmpBlock, /*markCached=*/false);
        FinishAccess();
      }
    }
    
    _ORAM::Address GetAddress(_ORAM::Address addr, uint64_t l) {
      return (addr >> (L_-l));
    }

    void BeginAccess(_ORAM::Address addr, Block_t& ret, bool markCached = true) {
      _ORAM::Position currPos = 0;
      _ORAM::Position newPos = 0;
      _ORAM::Position newChildPos = 0;
      _ORAM::ORAMAddress oaddr;

      for (uint64_t l=0; l<L_; l++) {
        newPos = newChildPos;
        newChildPos = UniformRandom(std::min(uint64_t{N_-1}, (uint64_t{1}<<(l))-1));

        _ORAM::Address levelAddress = GetAddress(addr, l);
        
        // std::cout << l << " " << levelAddress << " " << currPos << " " << newPos << " " << newChildPos << " " << std::endl;

        IntermidiateBlock_t ib;
        oaddr = _ORAM::ORAMAddress{levelAddress, currPos};
        orams[l]->BeginAccess(oaddr, newPos, ib, /*cache=*/true);
        if ((GetAddress(addr, l+1) & 1) == 0) {
          currPos = ib.data.pos[0];
          ib.data.pos[0] = newChildPos;
        } else {
          currPos = ib.data.pos[1];
          ib.data.pos[1] = newChildPos;
        }
        orams[l]->WriteToSomeStashedBlock(levelAddress, ib, /*markUncached=*/true);
        orams[l]->FinishAccess();
      }

      oaddr = _ORAM::ORAMAddress{addr, currPos};
      lastORAM.BeginAccess(oaddr, newChildPos, ret, markCached);
    }

    void WriteToSomeStashedBlock(const Address address, const Block_t &block, const bool markUncached = true) {
      lastORAM.WriteToSomeStashedBlock(address, block, markUncached);
    }

    void FinishAccess() {
      lastORAM.FinishAccess();
    }

    void UncacheAll() {
      lastORAM.UncacheAll();
    }

    void Read(_ORAM::Address addr, T& out) {
      Block_t aux;
      BeginAccess(addr, aux, /*markCached=*/false);
      out = aux.data;
      lastORAM.FinishAccess();
    }

    void Write(_ORAM::Address addr, const T& in) {
      Block_t aux;
      BeginAccess(addr, aux, /*markCached=*/true);
      lastORAM.WriteToSomeStashedBlock(addr, Block_t{in}, /*markUncached=*/true);
      lastORAM.FinishAccess();
    }
  }; // struct RecursiveORAM
} // namespace _ORAM::RecursiveORAM


template<>
INLINE void CMOV<_ORAM::RecursiveORAM::Data>(const uint64_t& condition, _ORAM::RecursiveORAM::Data& A, const _ORAM::RecursiveORAM::Data& B) {
  CMOV(condition, A.pos[0], B.pos[0]);
  CMOV(condition, A.pos[1], B.pos[1]);
}