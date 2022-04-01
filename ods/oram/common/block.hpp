#pragma once
#include "common/defs.hpp"
#include "common/mov_intrinsics.hpp"
#include "common/encrypted.hpp"
#include "oram/common/types.hpp"
#include "oram/common/concepts.hpp"

namespace _ORAM::Block {
  struct DefaultBlockData {
    uint64_t v;
    static consteval inline DefaultBlockData DUMMY() {return DefaultBlockData{static_cast<uint64_t>(-1)}; }
    bool operator==(const DefaultBlockData &other) const = default;
    friend std::ostream &operator<<(std::ostream &o, const DefaultBlockData &x) { o << x.v; return o; }
  };
  using DefaultData_t = DefaultBlockData;


  template<typename T=DefaultData_t, bool ENCRYPT_BLOCKS=ORAM__ENCRYPT_BLOCKS>
  requires (IS_POD<T>())
  struct Block
  {
    using _T = T;

    T data;

    auto operator==(const Block &other) const
    {
      return data == other.data;
    }
    friend std::ostream &operator<<(std::ostream &o, const Block &x)
    {
      o << "(" << x.data << ")";
      return o;
    }
    static consteval inline Block DUMMY() { return Block{MakeDummy<T>()}; }

    using Encrypted_t = std::conditional_t<ENCRYPT_BLOCKS, Encrypted<Block>, NonEncrypted<Block>>;
  }; // struct Block

  static_assert(IS_POD<Block<DefaultData_t> >());
  static_assert(IS_POD<Block<DefaultData_t>::Encrypted_t >());
} // namespace _ORAM::Block

namespace _ORAM::StashedBlock
{  
  template<typename Block=Block::Block<Block::DefaultData_t> >
  requires ::Concepts::Encryptable<Block>
  struct StashedBlock
  {
    using Block_t = Block;

    bool cached;
    ORAMAddress oaddress;
    Block block;

    friend std::ostream &operator<<(std::ostream &o, const StashedBlock &x)
    {
      o << "(";
      if (x.cached) {
        o << "cached, ";
      }
      o << x.oaddress << ", " << x.block << ")";
      return o;
    }

    static consteval inline StashedBlock DUMMY() { return 
      StashedBlock{
        /*cached=*/false,
        ORAMAddress::DUMMY(),
        Block::DUMMY()
      };
    }
  }; // struct StashedBlock

  static_assert(IS_POD<StashedBlock<Block::Block<Block::DefaultData_t> > >());

} // namespace _ORAM::StashedBlock

template<>
INLINE void CMOV<_ORAM::Block::DefaultBlockData>(const uint64_t& condition, _ORAM::Block::DefaultBlockData& A, const _ORAM::Block::DefaultBlockData& B) {
  CMOV(condition, A.v, B.v);
}

template<typename T, bool W>
INLINE void CMOV(const uint64_t& condition, _ORAM::Block::Block<T,W>& A, const _ORAM::Block::Block<T,W>& B) {
  CMOV(condition, A.data, B.data);
}

template<typename T>
INLINE void CMOV(const uint64_t& condition, _ORAM::StashedBlock::StashedBlock<T>& A, const _ORAM::StashedBlock::StashedBlock<T>& B) {
  CMOV(condition, A.cached, B.cached);
  CMOV(condition, A.oaddress, B.oaddress);
  CMOV(condition, A.block, B.block);
}

// UNDONE(): figure out how to not have to overload CXCHG and TSET here:
//
OVERLOAD_TSET_CXCHG(_ORAM::StashedBlock::StashedBlock<T>, typename T)

OVERLOAD_TSET_CXCHG(_ORAM::Block::Block<T COMMA W>, typename T, bool W)


// UNDONE: move this to tthe right file, we just moved here beause we know 
// the include order won't break.
//
template <typename T>
struct TaggedT {
  bool isDummy;
  uint64_t tag;
  T v;

  static consteval TaggedT DUMMY() {
    TaggedT ret;
    ret.isDummy = true;
    ret.tag = static_cast<uint64_t>(-1);
    ret.v = T::DUMMY();
    return ret;
  }
};


template<typename T>
INLINE void CMOV(const uint64_t& condition, TaggedT<T>& A, const TaggedT<T>& B) {
  CMOV(condition, A.isDummy, B.isDummy);
  CMOV(condition, A.tag, B.tag);
  CMOV(condition, A.v, B.v);
}

OVERLOAD_TSET_CXCHG(TaggedT<T>, typename T);