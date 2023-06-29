#pragma once
#include "common/defs.hpp"
#include "common/dummy.hpp"
#include "common/encrypted.hpp"
#include "common/mov_intrinsics.hpp"
#include "oram/common/concepts.hpp"
#include "oram/common/types.hpp"
#ifndef ELEMENT_SIZE
#define ELEMENT_SIZE 128
#endif

namespace EM::Algorithm {
template <bool perf = true>
INLINE void condSwap(const auto& cond, auto& v1, auto& v2) {
  if constexpr (perf) {
    PERFCTR_INCREMENT(swapCount);
  }

  obliSwap(cond, v1, v2);
  // CXCHG(cond, v1, v2);
}

INLINE void swap(auto& v1, auto& v2) {
  // PERFCTR_INCREMENT(swapCount);
  std::swap(v1, v2);
}
// enum DummyFlag { NOT_DUMMY, DUMMY_DEFAULT, DUMMY_LESS, DUMMY_GE };
struct DefaultBlockData {
  uint64_t v;

  static consteval inline DefaultBlockData DUMMY() {
    return DefaultBlockData{static_cast<uint64_t>(-1)};
  }
  bool operator==(const DefaultBlockData& other) const = default;
  bool operator<(const DefaultBlockData& other) const { return v < other.v; }
#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const DefaultBlockData& x) {
    o << x.v;
    return o;
  }
#endif
};
using DefaultData_t = DefaultBlockData;

template <typename T = DefaultData_t,
          bool ENCRYPT_BLOCKS = ORAM__ENCRYPT_BLOCKS>
  requires(IS_POD<T>())
struct Block {
  using _T = T;
  static constexpr size_t paddingSize = sizeof(T) % 32 == 16 ? 8 : 0;
  T data;
  uint32_t tag;  // tie breaker
  bool dummyFlag;
  bool lessFlag;
  char padding[paddingSize];
  auto operator==(const Block& other) const { return data == other.data; }
  bool operator<(const Block& other) const {
    return (data < other.data) | ((data == other.data) & (tag < other.tag));
  }
#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const Block& x) {
    o << "(" << x.data << ")";
    return o;
  }
#endif
  inline void setData(const T& _data) {
    this->data = _data;
    this->tag = UniformRandom32();
    this->dummyFlag = false;
  }

  inline const T& getData() { return data; }

  static consteval inline Block DUMMY() {
    return Block{T::DUMMY(), 0, true, false};
  }
  inline bool isDummy() const { return dummyFlag; }
  inline bool setAndGetMarked(const Block& pivot) {
    return lessFlag = *this < pivot;
  }
  inline bool isMarked(const Block& unused) const { return lessFlag; }

  inline bool isLess() const { return lessFlag; }
  inline void setLessFlag(bool flag) { this->lessFlag = flag; }
  inline void condChangeMark(bool cond, const Block& unused) {
    this->lessFlag ^= cond;
  }
  inline void setDummy() { setDummyFlag(true); }
  inline void setDummyFlag(bool flag) { this->dummyFlag = flag; }
  inline void setDummyFlagCond(bool cond, bool flag) {
    CMOV(cond, this->dummyFlag, flag);
  }

  using Encrypted_t =
      std::conditional_t<ENCRYPT_BLOCKS, Encrypted<Block>, NonEncrypted<Block>>;
};  // struct Block
template <typename T>
struct TaggedT {
  static constexpr size_t paddingSize = sizeof(T) % 32 == 16 ? 8 : 0;
  uint64_t tag;
  T v;
  char padding[paddingSize];

  // UNDONE(): this function needs to stop being a member function and be a
  // template function on it's own. The main reason is that T::DUMMY() for int
  // is not callable.
  //
  static consteval TaggedT DUMMY() {
    TaggedT ret;
    ret.tag = static_cast<uint64_t>(-1);
    ret.v = T();
    return ret;
  }

  inline void setData(const T& _data) { v = _data; }

  inline const T& getData() { return v; }

  inline bool isDummy() const { return tag >> 63; }

  inline void setDummy() { tag |= 0x8000'0000'0000'0000UL; }

  inline void setTag(uint64_t _tag) { tag = _tag & 0x7fff'ffff'ffff'ffffUL; }

  inline bool setAndGetMarked(uint64_t bitMask) const {  // same as isMarked
    return isMarked(bitMask);
  }

  inline bool isMarked(uint64_t bitMask) const { return !(tag & bitMask); }

  inline void condChangeMark(bool cond, uint64_t bitMask) {
    CMOV(cond, tag, tag ^ bitMask);
  }

  inline uint8_t getMarkAndUpdate(uint64_t k) {
    uint64_t realTag = tag & 0x7fff'ffff'ffff'ffffUL;
    tag &= 0x8000'0000'0000'0000UL;
    tag |= realTag / k;
    uint8_t mark = realTag % k;
    return mark;
  }
};

static_assert(IS_POD<Block<DefaultData_t>>());
static_assert(IS_POD<Block<DefaultData_t>::Encrypted_t>());
}  // namespace EM::Algorithm

struct SortElement {
  uint64_t key;
  char payload[ELEMENT_SIZE - sizeof(key)];
  static consteval inline SortElement DUMMY() {
    return SortElement{static_cast<uint64_t>(-1)};
  }
  bool operator==(const SortElement& other) const { return key == other.key; }
  bool operator<(const SortElement& other) const { return key < other.key; }
#ifndef ENCLAVE_MODE
  friend std::ostream& operator<<(std::ostream& o, const SortElement& x) {
    o << x.key;
    return o;
  }
#endif
};

INLINE void CMOV(const uint64_t& condition, SortElement& A,
                 const SortElement& B) {
  CMOV(condition, A.key, B.key);
  for (uint64_t i = 0; i < sizeof(SortElement::payload); i += 8) {
    CMOV(condition, *(uint64_t*)&(A.payload[i]), *(uint64_t*)&(B.payload[i]));
  }
}

template <>
INLINE void CMOV<EM::Algorithm::DefaultBlockData>(
    const uint64_t& condition, EM::Algorithm::DefaultBlockData& A,
    const EM::Algorithm::DefaultBlockData& B) {
  CMOV(condition, A.v, B.v);
}

template <typename T, bool W>
INLINE void CMOV(const uint64_t& condition, EM::Algorithm::Block<T, W>& A,
                 const EM::Algorithm::Block<T, W>& B) {
  CMOV(condition, A.data, B.data);
  CMOV(condition, A.tag, B.tag);
  CMOV(condition, A.dummyFlag, B.dummyFlag);
  CMOV(condition, A.lessFlag, B.lessFlag);
}

template <typename T>
INLINE void CMOV(const uint64_t& condition, EM::Algorithm::TaggedT<T>& A,
                 const EM::Algorithm::TaggedT<T>& B) {
  CMOV(condition, A.tag, B.tag);
  CMOV(condition, A.v, B.v);
}

// UNDONE(): figure out how to not have to overload CXCHG and TSET here:

OVERLOAD_TSET_CXCHG(EM::Algorithm::Block<T COMMA W>, typename T, bool W)
OVERLOAD_TSET_CXCHG(EM::Algorithm::TaggedT<T>, typename T);