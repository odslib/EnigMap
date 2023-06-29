#pragma once
#include <inttypes.h>
#include <map>
#include <cassert>
#include <type_traits>
#include "common/defs.hpp"
#include "common/tracing/tracer.hpp"
#include "common/mov_intrinsics.hpp"
#include "oram/common/types.hpp" //UNDONE(): remove this import?
#include "otree/node.hpp"
#include "otree/oram_interface.hpp"
#include "external_memory/emvector.hpp"
#include "external_memory/algorithm/sort.hpp"

// This file implements https://eprint.iacr.org/2014/185.pdf with our optimizations.
namespace _OBST {

namespace OBST {
  template <typename T>
  concept HasBlockOfNode = std::is_same_v<typename T::_T, Node>;
  template<typename OramClient /*=_ORAM::OramClient::OramClient<Node,ORAM__Z,ORAM__S,false,false,4>*/ > 
  requires HasBlockOfNode<OramClient> 
  struct OBST {
    using Block_t = typename OramClient::Block_t;
    using EMPoint_t = typename std::pair<K,V>;
    using Vector_t = typename EM::Vector::Vector<EMPoint_t>;
    using StashedBlock_t = typename _ORAM::StashedBlock::StashedBlock<typename OramClient::Block_t>;
    #ifndef NDEBUG
    std::map<K,V> m;
    #endif
    uint64_t maxDepth = 0;
    uint64_t maxNodes = 0;
    uint64_t count = 0;

    _ORAM::Index nextPosition = 0; //Index for the next generated position
    _ORAM::ORAMAddress root;
    _ORAM::ORAMAddress FAKE_NODE;
    OramClient oram;

    private:
     explicit OBST(uint64_t _maxNodes, bool noInit, char)
    : maxNodes(GetNextPowerOfTwo(_maxNodes + 2)),
      maxDepth(GetLogBaseTwo(GetNextPowerOfTwo(_maxNodes + 2))), 
      oram(_maxNodes + 2, noInit)
    {
      PROFILE_F(); 
      // UNDONE(): this bound can be improved, it can be 1.44, not 1.5;
      //
      maxDepth += (maxDepth >> 1) + 1;

      count = 0;
    }

    public:
    explicit OBST(uint64_t _maxNodes, bool noInit=false) : OBST(_maxNodes,noInit,'.')
    {
      Block_t  NodeV;

      // The fake node is the node used for dummy accesses, we need to have it
      // cached because a lot of nodes will point to it.
      //
      FAKE_NODE.position = oram.GenRandomPosition();
      oram.GetNextFreeBlockOp(FAKE_NODE.address, NodeV, FAKE_NODE.position);
      NodeV = Block_t {
          Node{INVALID_KEY, INVALID_VALUE, {FAKE_NODE, FAKE_NODE}, B_BALANCED}};
      // used FAKE_NODE.address for assertion
      oram.FinishOramOp(FAKE_NODE.address, NodeV);

      // The root is a node where we always pick the left child to travel to. The
      // value on the root should not be a returnable value.
      //
      root.position = oram.GenRandomPosition();
      oram.GetNextFreeBlockOp(root.address, NodeV, root.position);
      NodeV = Block_t {
          Node{INVALID_KEY, INVALID_VALUE, {FAKE_NODE, FAKE_NODE}, B_BALANCED}};
      // used root.address for assertion
      oram.FinishOramOp(root.address, NodeV);
    }

    
    // Full ORAM initialization
    explicit OBST(uint64_t _maxNodes, Vector_t& points) : OBST(_maxNodes, '.') {
      Assert (points.size() > 0);
      Assert(_maxNodes >= points.size());
      // Stage 1) 
      // Sort the nodes
      Assert(
        EM::Algorithm::IsSorted(points, [](const EMPoint_t& p1, const EMPoint_t& p2) {
          return p1.first < p2.first;
        })
      );

      // Give addresses and positions to nodes:
      EM::Vector::Vector<StashedBlock_t> auxv(points.size() + 2);

      FAKE_NODE.address = 0;
      FAKE_NODE.position = oram.GenRandomPosition();
      for (uint64_t i=0; i<points.size(); i++) {
        StashedBlock_t& sb = auxv[i];
        sb.oaddress.address = i+2;
        sb.oaddress.position = oram.GenRandomPosition();        
        sb.block.data.balance = B_BALANCED;
        sb.block.data.child[0] = FAKE_NODE;
        sb.block.data.child[1] = FAKE_NODE;
        sb.block.data.k = points[i].first;
        sb.block.data.v = points[i].second;
      }

      uint64_t height;
      _ORAM::ORAMAddress rootAddress;
      // Build the vector recursively:
      //
      std::tie(height, rootAddress) = _RecursiveConstruct(auxv, 0, points.size()-1);

      // Setup root block and FAKENODE:
      //
      auxv[points.size()].oaddress = FAKE_NODE;
      auxv[points.size()].block = Block_t {
          Node{INVALID_KEY, INVALID_VALUE, {FAKE_NODE, FAKE_NODE}, B_BALANCED}};
      
      auxv[points.size()+1].oaddress = _ORAM::ORAMAddress{1,oram.GenRandomPosition()};
      auxv[points.size()+1].block = Block_t {
          Node{INVALID_KEY, INVALID_VALUE, {rootAddress, FAKE_NODE}, B_BALANCED}};
      

      // Stage 2) Evict blocks level by level to ORAM
      //
      oram.InitViaLargeStashVector(auxv);
    }


    // It's ok for the recursive construct to leak information from the if, as the access pattern only depends on 
    // the size of the vector, and not on the actual content of the vector.
    //
    std::pair<uint64_t, _ORAM::ORAMAddress> _RecursiveConstruct(EM::Vector::Vector<StashedBlock_t>& points, uint64_t l, uint64_t r) {
      Assert(l <= r);
      uint64_t m = (l+r)/2;
      
      // UNDONE(): add a reference counter the current node to avoid having
      // to do several page reads on tree top pages for very large trees.
      //
      uint64_t leftH = 0;
      uint64_t rightH = 0;
      _ORAM::ORAMAddress childs[2] = {FAKE_NODE,FAKE_NODE};

      if (m > l) {
        std::tie(leftH, childs[0]) = _RecursiveConstruct(points, l, m-1);
      }
      
      if (r > m) {
        std::tie(rightH, childs[1]) = _RecursiveConstruct(points, m+1, r);
      }

      StashedBlock_t& nodeV = points[m];
      nodeV.block.data.child[0] = childs[0];
      nodeV.block.data.child[1] = childs[1];
      nodeV.block.data.balance = B_BALANCED;
      if (rightH > leftH) { nodeV.block.data.balance = B_RIGHT; Assert(rightH-leftH == 1); }
      if (rightH < leftH) { nodeV.block.data.balance = B_LEFT; Assert(leftH-rightH == 1); }        
      std::pair<uint64_t,_ORAM::ORAMAddress> ret = std::pair<uint64_t,_ORAM::ORAMAddress>{1 + std::max(rightH, leftH), nodeV.oaddress};
      return ret;
    }

    const bool Get(const K &k, V &ret) {
      PROFILE_F();

      _ORAM::ORAMAddress curr = root;
      _ORAM::Address prev;
      Block_t currBlock;
      _ORAM::Position newPos, newPos0;
      bool contained = false;
      newPos = newPos0 = oram.GenRandomPosition();
      for (int i = 0; i < maxDepth; i++) {

        // After we hit a fake node, we just keep accessing fake path's:
        //
        prev = curr.address;
        oram.BeginReadOp(curr.address, curr.position, currBlock, newPos);

        CMOV(curr.address == FAKE_NODE.address, FAKE_NODE.position, newPos);

        // Find where to go next:
        //
        Dir_t dir = Direction(currBlock.data, k);
        bool match = dir == B_BALANCED;
        bool right = dir == B_RIGHT;
        newPos = oram.GenRandomPosition();
        TSET(right, curr, currBlock.data.child[0], currBlock.data.child[1]);
        CMOV(!right, currBlock.data.child[0].position, newPos);
        CMOV(right, currBlock.data.child[1].position, newPos);

        // Local updates:
        CMOV(match, contained, true);
        CMOV(match, ret, currBlock.data.v);

        // Commit node value:
        oram.FinishOramOp(prev, currBlock);
        
        CMOV(curr.address == FAKE_NODE.address, curr.position, FAKE_NODE.position);
      }
      Assert(curr.address == FAKE_NODE.address);
      root.position = newPos0;

      // Testing code:
      //
      #ifndef NDEBUG
      if (m.count(k) > 0) {
        ret = m[k];
        return true;
      } else {
        return false;
      }
      #endif

      return contained;
    }

    // For insertions we need to do two passes:
    // On the first pass we do the actual insertion and find out what kind of
    // rotation we need to do. On the second pass we do the actual rotation and fix
    // the balances. On both passes, we only need to read nodes on the path to the
    // insertion point.
    //

    #define GO_NEXT(curr, currBlock, newPos)                                       \
    do {                                                                         \
      dir = Direction(currBlock.data, k);                                        \
      right = dir == B_RIGHT;                                                    \
      match = dir == B_BALANCED;                                                 \
      TSET(right, curr, currBlock.data.child[0], currBlock.data.child[1]);     \
      TSET(right, other, currBlock.data.child[1], currBlock.data.child[0]);    \
      CMOV(!right, currBlock.data.child[0].position, newPos);                    \
      CMOV(right, currBlock.data.child[1].position, newPos);                     \
      TSET(right, nCurr, currBlock.data.child[0], currBlock.data.child[1]);    \
    } while (false);

    #define GO_NEXT_COND(cond, curr, currBlock, newPos)                                       \
    do {                                                                         \
      CMOV(cond,        dir, Direction(currBlock.data, k));                             \
      CMOV(cond,        right, dir == B_RIGHT);                                         \
      CMOV(cond,        match, dir == B_BALANCED);                                      \
      CTSET(cond, right,  curr, currBlock.data.child[0], currBlock.data.child[1]);     \
      CTSET(cond, right,  other, currBlock.data.child[1], currBlock.data.child[0]);    \
      CMOV(cond*!right, currBlock.data.child[0].position, newPos);                    \
      CMOV(cond*right,  currBlock.data.child[1].position, newPos);                     \
      CTSET(cond, right,  nCurr, currBlock.data.child[0], currBlock.data.child[1]);    \
    } while (false);

    #define NO_REBALANCE 0
    #define DIRECT_REBALANCE 1
    #define ROT2 2
    #define ROT3 3

    void Insert(const K &k, const V &v) {
      PROFILE_F();
      #ifndef NDEBUG
      m[k] = v;
      #endif

      _ORAM::ORAMAddress currP;
      _ORAM::ORAMAddress curr = root;
      _ORAM::ORAMAddress nCurr;
      _ORAM::Address prev;
      Block_t currBlock;
      _ORAM::Position newPos, newPos0, newPos1;
      _ORAM::ORAMAddress newBlock;
      Block_t newBlockV;

      Dir_t dir = B_BALANCED;
      Dir_t pDir1 = B_BALANCED;
      Dir_t pDir2;

      _ORAM::ORAMAddress A, B, C, D, E, F, G, X;
      _OBST::V vB, vF, vD;
      _OBST::K kB, kF, kD;

      // Values saved for second pass:
      //
      bool match = false;
      int rotType = NO_REBALANCE;
      int unState = 10;
      Dir_t dirArg1 = B_BALANCED;
      Dir_t dirArg2 = B_BALANCED;
      _ORAM::ORAMAddress other;
      //

      bool contained = false;
      bool inserted = false;
      bool right = false;



      // First pass:
      //
      newBlock.position = oram.GenRandomPosition();
      oram.GetNextFreeBlockOp(newBlock.address, newBlockV, newBlock.position);
      newBlockV.data.balance = B_BALANCED;
      newBlockV.data.k = k;
      newBlockV.data.v = v;
      newBlockV.data.child[0] = FAKE_NODE;
      newBlockV.data.child[1] = FAKE_NODE;
      oram.FinishOramOp(newBlock.address, newBlockV, false);

      newPos = newPos0 = oram.GenRandomPosition();
      newPos1 = oram.GenRandomPosition();
      for (int i = 0; i < maxDepth; i++) {
        // After we hit a fake node, we just keep accessing fake path's:
        //
        prev = curr.address;
        oram.BeginReadOp(curr.address, curr.position, currBlock, newPos);
        CMOV(curr.address == FAKE_NODE.address, FAKE_NODE.position, newPos);

        // Find where to go next:
        //
        dir = Direction(currBlock.data, k);
        pDir2 = pDir1;
        pDir1 = dir;
        currP = curr;
        currP.position = newPos;
        newPos = oram.GenRandomPosition();
        GO_NEXT(curr, currBlock, newPos);

        // Compute the parent information for rebalancing:
        //
        if (i > 0) {
          bool case0 = 
            (currBlock.data.balance != B_BALANCED)
          * (currP.address != FAKE_NODE.address);
          bool case0_0 = case0 
                * (currBlock.data.balance != pDir1);
          bool case1 = !case0 * (unState == 1);
          bool case1_1 = case1 * (pDir1 == pDir2);
          bool case1_1_1 = case1_1 * (nCurr.address == FAKE_NODE.address);
          bool case2 = !case0 * !case1 * (unState == 2);
          bool case2_1 = case2 * (D.address != FAKE_NODE.address);
          bool case2_1_1 = case2_1 * (nCurr.address == FAKE_NODE.address);
          bool case2_1_2 = case2_1 * (dir != dirArg1);
          bool case2_2 = case2 * !case2_1;
          // if (case0) {
            CMOV(case0,   unState, 0);
            CMOV(case0,   B,       currP);
            CMOV(case0,   A,       other);
            CMOV(case0,   F,       nCurr);
            CMOV(case0,   vB,      currBlock.data.v);
            CMOV(case0,   kB,      currBlock.data.k);
            CMOV(case0_0, rotType, DIRECT_REBALANCE);
            CMOV(case0_0, unState, 3);
          // } else {
            // if (unState == 1) {
              CMOV(case1, G,       other);
              CMOV(case1, D,       nCurr);
              CMOV(case1, vF,      currBlock.data.v);
              CMOV(case1, kF,      currBlock.data.k);
              CMOV(case1, dirArg1, pDir2);
              // if (pDir1 == pDir2) {
                CMOV(case1_1, rotType, ROT2);
                CMOV(case1_1, unState, 3);
                // if (D.address == FAKE_NODE.address) {
                  CMOV(case1_1_1, D, newBlock);
                // }
              // }
            // } else if (unState == 2) {
              // if (D.address != FAKE_NODE.address) {
                CMOV(case2_1, C, other);
                CMOV(case2_1, E, nCurr);

                // if (nCurr.address == FAKE_NODE.address) {
                  CMOV(case2_1_1, E, newBlock);
                // }

                // if (dir != dirArg1) {
                  CXCHG(case2_1_2, C, E);
                // }
                CMOV(case2_1, vD,      currBlock.data.v);
                CMOV(case2_1, kD,      currBlock.data.k);
                CMOV(case2_1, rotType, ROT3);
                CMOV(case2_1, dirArg2, dir);
              // } else {
                CMOV(case2_2, E,       FAKE_NODE);
                CMOV(case2_2, C,       FAKE_NODE);
                CMOV(case2_2, D,       newBlock);
                CMOV(case2_2, vD,      v);
                CMOV(case2_2, kD,      k);
                CMOV(case2_2, rotType, ROT3);
                CMOV(case2_2, dirArg2, B_BALANCED);
              // }
            // }
          // }
          unState++;
        }
        CMOV(match, contained, true);
        CMOV(match, currBlock.data.v, v);
        bool cond0 = !contained
                   * !inserted
                   * (curr.address == FAKE_NODE.address);
        // if (!contained && !inserted && curr.address == FAKE_NODE.address) {
          // We don't access the new node on this pass,
          // as it is not used.
          //
          CMOV(cond0*!right, currBlock.data.child[0], newBlock);
          CMOV(cond0*right,  currBlock.data.child[1], newBlock);
          CMOV(cond0,        curr,                    FAKE_NODE);
          CMOV(cond0,        inserted,                true);
        // }

        // Commit node value:
        oram.FinishOramOp(prev, currBlock);
        CMOV(curr.address == FAKE_NODE.address, curr.position,
            FAKE_NODE.position);
      }
      Assert(curr.address == FAKE_NODE.address);
      root.position = newPos0;
      Assert(unState >= 3 || !inserted);
      oram.FreeBlockMaybeFOp(newBlock.address, newBlock.position, !inserted);

      // Second pass:
      //

      // UNDONE(): Split this code into functions?
      //
      newPos = oram.GenRandomPosition();
      bool rebalancing = false;
      CMOV(inserted, curr, root);
      CMOV(inserted, root.position, newPos);
      CMOV(!inserted, curr, FAKE_NODE);

      for (int i = 0; i < maxDepth; i++) {
        // After we hit a fake node, we just keep accessing fake path's:
        //
        prev = curr.address;
        oram.BeginReadOp(curr.address, curr.position, currBlock, newPos);
        CMOV(curr.address == FAKE_NODE.address, FAKE_NODE.position, newPos);
        CMOV(currBlock.data.v == v, rebalancing, false);

        bool cond_dr = true
          * (rotType == DIRECT_REBALANCE)
          * (curr.address == B.address);
        bool cond_r2_1 = !cond_dr
          * (rotType == ROT2)
          * (curr.address == B.address);
        bool cond_r2_2 = !cond_dr * !cond_r2_1
          * (rotType == ROT2)
          * (curr.address == F.address);
        bool cond_r3_1 = !cond_dr * !cond_r2_1 * !cond_r2_2
          * (rotType == ROT3)
          * (curr.address == B.address);
        bool cond_r3_2 = !cond_dr * !cond_r2_1 * !cond_r2_2 * !cond_r3_1
          * (rotType == ROT3)
          * (curr.address == F.address);
        bool cond_r3_2_1 = cond_r3_2
          * (dirArg2 != B_BALANCED);
        bool cond_r3_2_1_1 = cond_r3_2_1
          * (dirArg1 == dirArg2);
        bool cond_r3_2_1_2 = cond_r3_2_1
          * !cond_r3_2_1_1;
        bool cond_r3_2_2 = cond_r3_2
          * !cond_r3_2_1;
        bool cond_r3_3 = !cond_dr * !cond_r2_1 * !cond_r2_2 * !cond_r3_1 * !cond_r3_2
          * (rotType == ROT3)
          * (curr.address == D.address);
        bool cond_r3_3_1 = cond_r3_3
          * (dirArg2 != dirArg1)
          * (dirArg2 != B_BALANCED);
        bool cond_r3_3_2 = cond_r3_3 * (!cond_r3_3_1);
        bool cond_update = !cond_dr * !cond_r2_1 * !cond_r2_2 * !cond_r3_1 * !cond_r3_2 * !cond_r3_3;
        bool cond_update_1 = cond_update * rebalancing;
        bool cond_update_2 = cond_update * !cond_update_1
          * (i == 0)
          * (rotType == NO_REBALANCE);
        // if (rotType == DIRECT_REBALANCE && curr.address == B.address) {
          CMOV(cond_dr, currBlock.data.balance, B_BALANCED);
          CMOV(cond_dr, rebalancing,            true);
          CMOV(cond_dr, newPos,                 oram.GenRandomPosition());
          GO_NEXT_COND(cond_dr, curr, currBlock, newPos);
          // currBlock.data.balance = dir;
        // } else if (rotType == ROT2 && curr.address == B.address) {
          CMOV(cond_r2_1, currBlock.data.v, vF);
          CMOV(cond_r2_1, currBlock.data.k, kF);
          CMOV(cond_r2_1, currBlock.data.balance, B_BALANCED);
          CMOV(cond_r2_1, curr, F);
          CMOV(cond_r2_1, newPos, oram.GenRandomPosition());
          CMOV(cond_r2_1, F.position, newPos);
          CMOV(cond_r2_1, newPos0, oram.GenRandomPosition());
          CMOV(cond_r2_1, newPos1, D.position);
          CMOV(cond_r2_1, D.position, newPos0);
          CTSET(cond_r2_1, !dirArg1, currBlock.data.child[0], F, D);
          CTSET(cond_r2_1, dirArg1, currBlock.data.child[1], F, D);
        // } else if (rotType == ROT2 && curr.address == F.address) {
          CMOV(cond_r2_2, currBlock.data.v, vB);
          CMOV(cond_r2_2, currBlock.data.k, kB);
          CMOV(cond_r2_2, currBlock.data.balance, B_BALANCED);
          CMOV(cond_r2_2, curr, D);
          CMOV(cond_r2_2, curr.position, newPos1);
          CMOV(cond_r2_2, newPos, newPos0);
          CTSET(cond_r2_2, !dirArg1, currBlock.data.child[0], A, G);
          CTSET(cond_r2_2, dirArg1, currBlock.data.child[1], A, G);
          CMOV(cond_r2_2, rebalancing, true);
        // } else if (rotType == ROT3 && curr.address == B.address) {
          CMOV(cond_r3_1, currBlock.data.v, vD);
          CMOV(cond_r3_1, currBlock.data.k, kD);
          CMOV(cond_r3_1, currBlock.data.balance, B_BALANCED);
          CMOV(cond_r3_1, curr, F);
          CMOV(cond_r3_1, newPos, oram.GenRandomPosition());
          CMOV(cond_r3_1, F.position, newPos);
          CMOV(cond_r3_1, newPos0, D.position);
          CMOV(cond_r3_1, D.position, oram.GenRandomPosition());
          CTSET(cond_r3_1, dirArg1, currBlock.data.child[0], D, F);
          CTSET(cond_r3_1, !dirArg1, currBlock.data.child[1], D, F);
        // } else if (rotType == ROT3 && curr.address == F.address) {
          CMOV(cond_r3_2, currBlock.data.v, vB);
          CMOV(cond_r3_2, currBlock.data.k, kB);
          // if (dirArg2 != B_BALANCED) {
            // if (dirArg1 == dirArg2) {
              CMOV(cond_r3_2_1_1, currBlock.data.balance, 1 - dirArg1);
              CMOV(cond_r3_2_1_1, X, E);
              CMOV(cond_r3_2_1_1, newPos1, oram.GenRandomPosition());
              CMOV(cond_r3_2_1_1, E.position, newPos1);
            // } else {
              CMOV(cond_r3_2_1_2, currBlock.data.balance, B_BALANCED);
              CMOV(cond_r3_2_1_2, X, C);
              CMOV(cond_r3_2_1_2, newPos1, oram.GenRandomPosition());
              CMOV(cond_r3_2_1_2, C.position, newPos1);
            // }
          // } else {
            CMOV(cond_r3_2_2, currBlock.data.balance, B_BALANCED);
            CMOV(cond_r3_2_2, X, FAKE_NODE);
          // }
          CMOV(cond_r3_2, newPos, D.position);
          CMOV(cond_r3_2, curr, D);
          CMOV(cond_r3_2, curr.position, newPos0);
          CTSET(cond_r3_2, !dirArg1, currBlock.data.child[0], A, C);
          CTSET(cond_r3_2, dirArg1, currBlock.data.child[1], A, C);
        // } else if (rotType == ROT3 && curr.address == D.address) {
          CMOV(cond_r3_3, currBlock.data.v, vF);
          CMOV(cond_r3_3, currBlock.data.k, kF);
          // if (dirArg2 != dirArg1 && dirArg2 != B_BALANCED) {
            CMOV(cond_r3_3_1, currBlock.data.balance, dirArg1);
          // } else {
            CMOV(cond_r3_3_2, currBlock.data.balance, B_BALANCED);
          // }
          CMOV(cond_r3_3, curr, X);
          CMOV(cond_r3_3, newPos, newPos1);
          CTSET(cond_r3_3, !dirArg1, currBlock.data.child[0], E, G);
          CTSET(cond_r3_3, dirArg1, currBlock.data.child[1], E, G);
          CMOV(cond_r3_3, rebalancing, true);
        // } else {
          CMOV(cond_update, newPos, oram.GenRandomPosition());
          GO_NEXT_COND(cond_update, curr, currBlock, newPos);
          
          // if (rebalancing) {
            CMOV(cond_update_1, currBlock.data.balance, dir);
          // } else if (i == 0 && rotType == NO_REBALANCE) {
            CMOV(cond_update_2, rebalancing, true);
          // }
        // }

        // Commit node value:
        oram.FinishOramOp(prev, currBlock);
        CMOV(curr.address == FAKE_NODE.address, curr.position,
            FAKE_NODE.position);
      }


      oram.UncacheAll();
    }

    struct DeletionMetadata {
      #define DELTP_NOTFOUND 0
      #define DELTP_LEAF 1
      #define DELTP_SINGLE_CHILD 2
      #define DELTP_TWO_CHILDREN 3

      uint64_t deletionType = DELTP_NOTFOUND; 

      _ORAM::ORAMAddress writeAddress = _ORAM::ORAMAddress::DUMMY();
      Node writeValue = Node::DUMMY();

      _ORAM::ORAMAddress deletionAddress = _ORAM::ORAMAddress::DUMMY();
      
      _ORAM::ORAMAddress replacementChildren = _ORAM::ORAMAddress::DUMMY();
    };

    // void _Delete_Step1(bool enabled, const K& k, DeletionMetadata& md, _ORAM::ORAMAddress& curr) {
    //   _ORAM::Position newPos = oram.GenRandomPosition();
    //   Dir_t nextDir;
    //   if (!enabled) { Assert(curr.address == FAKE_NODE.address); }
    //   Block_t currBlock;
    //   _ORAM::ORAMAddress& next;

    //   oram.BeginReadOp(curr.address, curr.position, currBlock, newPos);
    //   curr.position = newPos;
    //   CMOV(curr.address == FAKE_NODE.address, FAKE_NODE.position, newPos);
    //   nextDir = Direction(currBlock, k);
    //   if (nextDir == 2) {
    //     // This is the node to delete:
    //     //
    //     if (currBlock.child[0].address == FAKE_NODE.address && currBlock.child[1].address == FAKE_NODE.address) {
    //       md.deleteType = DELTP_LEAF;
    //       md.deletionAddress = curr;
    //     } else if (currBlock.child[0].address == FAKE_NODE.address || currBlock.child[1].address == FAKE_NODE.address) {
    //       md.deleteType = DELTP_SINGLE_CHILD;
    //       md.deletionAddress = curr;
    //     } else {
    //       md.writeAddress = DELTP_SINGLE_CHILD;
    //       md.deletionAddress = curr;
    //     }
    //   } else {
    //     next = currBlock.data.child[nextDir];
    //   }
    //   oram.FinishOramOp(curr.address, currBlock);

    //   // Recurse:
    //   //
    //   _Delete_Step1(enabled, k, next);

    //   oram.BeginReadOp(curr.address, curr.position, currBlock, newPos);
    //   curr.position = newPos;
    //   CMOV(curr.address == FAKE_NODE.address, FAKE_NODE.position, newPos);
    //   // Do some computation.
    //   oram.FinishOramOp(curr.address, currBlock);
    // }

    // void _Delete_Step2(bool enabled, const K& k) {

    // }

    // void Delete(const K &k) {
    //   // UNDONE(): make the allocator allow to do optional deletions.
    //   // Also, make the allocator use the free list, which it currently isn't using.
    //   // 

    //   // First iteration finds the deletion point, the type of deletion,
    //   // and possibly also the successor of the deletion point.
    //   //
    //   _Delete_Step1(true, k);

    //   // Second iteration writes to the deletion point, and removes the deleted
    //   // node.
    //   //
    //   _Delete_Step2(true, k);
    // }

    void RecursivePrint() {
      std::string prefix = "-";
      _ORAM::Position newPos = oram.GenRandomPosition();
      _ORAM::ORAMAddress addr;
      Block_t currV;
      oram.BeginReadOp(root.address, root.position, currV, newPos);
      root.position = newPos;
      newPos = oram.GenRandomPosition();
      addr = currV.data.child[0];
      currV.data.child[0].position = newPos;
      oram.FinishOramOp(root.address, currV);
      _RecursivePrint(addr, newPos, prefix, " > ", B_BALANCED);
      X_LOG_SIMPLE("\n");
    }

  private:
    INLINE Dir_t Direction(const Node& head, const K& k) {
      Dir_t ret = B_LEFT;
      CMOV(k == head.k, ret, B_LEFT);
      CMOV(k > head.k, ret, B_RIGHT);
      return ret;
    }
    

    int _RecursivePrint(_ORAM::ORAMAddress curr, _ORAM::Position newPos, const std::string& prefix, const char* x, Dir_t direction) {
      if (curr.address == FAKE_NODE.address) return 0;

      _ORAM::Position newPosL = oram.GenRandomPosition();
      _ORAM::Position newPosR = oram.GenRandomPosition();
      Block_t currV;
      oram.BeginReadOp(curr.address, curr.position, currV, newPos);
      _ORAM::ORAMAddress left = currV.data.child[0];
      _ORAM::ORAMAddress right = currV.data.child[1];
      currV.data.child[0].position = newPosL;
      currV.data.child[1].position = newPosR;
      oram.FinishOramOp(curr.address, currV);

      int l1 = _RecursivePrint(left, newPosL, prefix + "  " + x[0], " ,|", 0);
      X_LOG_SIMPLE(prefix, "  ", x[1], "--+", ("/\\:"[currV.data.balance]), " ", currV.data.k,  " @ ", curr.address, ",", curr.position , " (", currV.data.v, ")\n");
      int l2 = _RecursivePrint(right, newPosR, prefix + "  " + x[2], "|` ", 1);
      if (l1 > l2) {
        Assert(currV.data.balance == B_LEFT);
      }
      if (l1 < l2) {
        Assert(currV.data.balance == B_RIGHT);
      }
      if (l1 == l2) {
        Assert(currV.data.balance == B_BALANCED);
      }
      Assert(abs(l1-l2) <= 1);
      return 1 + std::max(l1,l2);
    }
  }; // </struct OBST>
} // </namespace OBST>
}; // </namespace _OBST>