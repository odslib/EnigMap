#pragma once
#include <inttypes.h>
#include <sstream>
#include "common/defs.hpp"
#include "oram/common/types.hpp"
#include "oram/common/block.hpp"
namespace _OBST {
  typedef uint64_t K;
  typedef uint64_t V;
  constexpr inline K INVALID_KEY = -1;
  constexpr inline V INVALID_VALUE = -1;

  typedef uint32_t Dir_t;
  constexpr inline Dir_t B_LEFT = 0;
  constexpr inline Dir_t B_RIGHT = 1;
  constexpr inline Dir_t B_BALANCED = 2;

  inline std::string Dir_t_to_string(const Dir_t& t) {
    if (t == B_LEFT) {
      return "B_LEFT";
    }
    if (t == B_RIGHT) {
      return "B_RIGHT";
    }
    if (t == B_BALANCED) {
      return "B_BALANCED";
    }
    Assert(false);
    return "B_ERRORBALANCE";
  }

  struct Node;
  struct Node
  {
    K k;
    V v;
    _ORAM::ORAMAddress child[2];
    Dir_t balance;

    inline bool IsBalanced() { return balance == B_BALANCED; }
    auto operator==(const Node &other) const
    {
      return true * (k == other.k) * (v == other.v) * (child[0] == other.child[0]) * (child[1] == other.child[1]) * (balance == other.balance);
    }

    friend std::ostream &operator<<(std::ostream &o, const Node &x)
    {
      // UNDONE(): c++14 libraries with c++20 compiler don't have << working for non string types.
      //
      o << "(" << std::to_string(x.k) << ", " << std::to_string(x.v) << ", " << x.child[0] << ", " << x.child[1] << "," << Dir_t_to_string(x.balance) <<  ")";
      return o;
    }

    static consteval inline Node DUMMY() { return Node
      {
        K{INVALID_KEY},
        V{INVALID_VALUE},
        {_ORAM::ORAMAddress::DUMMY(), _ORAM::ORAMAddress::DUMMY()},
        Dir_t{B_BALANCED}
      };
    }

  };
  static_assert(IS_POD<Node>());
} // </namespace _OBST>

template<>
INLINE void CMOV<_OBST::Node>(const uint64_t& condition, _OBST::Node& A, const _OBST::Node& B) {
  CMOV(condition, A.v, B.v);
  CMOV(condition, A.k, B.k);
  CMOV(condition, A.v, B.v);
  CMOV(condition, A.child[0], B.child[0]);
  CMOV(condition, A.child[1], B.child[1]);
  CMOV(condition, A.balance, B.balance);
}