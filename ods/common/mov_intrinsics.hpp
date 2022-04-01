#pragma once
#include <utility>
#include <inttypes.h>
#include <typeinfo>
#include "cpp_extended.hpp"

INLINE void CSWAP8(const uint64_t swap, uint64_t *guy1, uint64_t *guy2) {
  asm volatile("test %%rdi, %%rdi\n\t"
               "mov (%%rsi), %%r10\n\t"
               "mov (%%rdx), %%r9\n\t"
               "mov %%r9, %%r11\n\t"
               "cmovnz %%r10, %%r9\n\t"
               "cmovnz %%r11, %%r10\n\t"
               "mov %%r9, (%%rdx)\n\t"
               "mov %%r10, (%%rsi)\n\t"
               : // no out
               : // no in
               : "r9", "r10", "r11");
}

INLINE void CMOV8_internal(const uint64_t cond, uint64_t& guy1, const uint64_t& guy2) {
  asm volatile("test %[mcond], %[mcond]\n\t"
               "cmovnz %[i2], %[i1]\n\t"
               : [i1] "=r"(guy1)
               : [mcond] "r"(cond), "[i1]" (guy1), [i2] "r"(guy2)
               : );
}

INLINE void CMOV4_internal(const uint64_t cond, uint32_t& guy1, const uint32_t& guy2) {
  asm volatile("test %[mcond], %[mcond]\n\t"
               "cmovnz %[i2], %[i1]\n\t"
               : [i1] "=r"(guy1)
               : [mcond] "r"(cond), "[i1]" (guy1), [i2] "r"(guy2)
               : );
}


INLINE void CMOV1(const bool& cond, uint8_t& val1, const uint8_t& val2) {
  uint32_t r1 = 0 | val1;
  uint32_t r2 = 0 | val2;
  CMOV4_internal(cond, r1, r2);
  val1 = r1 & 0xff;
}

INLINE void CMOV2(const bool& cond, uint16_t& val1, const uint16_t& val2) {
  uint32_t r1 = 0 | val1;
  uint32_t r2 = 0 | val2;
  CMOV4_internal(cond, r1, r2);
  val1 = r1 & 0xffff;
}

INLINE void CMOV4(const bool& cond, uint32_t& val1, const uint32_t& val2) {
  CMOV4_internal(cond,val1,val2);
}

INLINE void CMOV8(const bool& cond, uint64_t& val1, const uint64_t& val2) {
  CMOV8_internal(cond,val1,val2);
}

INLINE void CMOV_BOOL(const uint64_t& cond, bool& val1, const bool& val2) {
  uint32_t v1 = val1;
  CMOV4_internal(cond,v1,val2);
  val1 = v1;
}

// Note: Possible generic implementation of CMOV, that we don't use.
// template<typename T>
// void ObliMov(bool mov, T* guy1, T* guy2) {
//   static_assert(sizeof(T)%8 == 0);
//   uint64_t* curr1 = (uint64_t*)guy1;
//   uint64_t* curr2 = (uint64_t*)guy2;
//   for (uint64_t i = 0; i < sizeof(T) / 8; ++i) {
//     CMOV(mov, *curr1, *curr2);
//     curr1++;
//     curr2++;
//   }
// }

template<typename T>
INLINE void CMOV(const uint64_t& cond, T& val1, const T& val2) {
  Assert(false, "This should mov not be compiled check that you implemented CMOV"
    "and that if you used overloading that you called OVERLOAD_TSET_CXCHG. For type: ", typeid(T).name());
  if (cond) {
    val1 = val2;
  }
}

template<typename T>
INLINE void TSET(bool selector, T &A, const T &B, const T &C) {
  CMOV(selector, A, C);
  CMOV(!selector, A, B);
}

template<typename T>
INLINE void CTSET(bool condition, bool selector, T &A, const T &B, const T &C) {
  CMOV(condition*selector, A, C);
  CMOV(condition*!selector, A, B);
}

template<typename T>
INLINE void CXCHG(const uint64_t& cond, T& A, T& B) {
  const T C = A;
  CMOV(cond, A, B);
  CMOV(cond, B, C);
}


// Some CMOV specializations:
//
template<>
INLINE void CMOV<uint64_t>(const uint64_t& cond, uint64_t& val1, const uint64_t& val2) {
  CMOV8(cond, val1, val2);
}

template<>
INLINE void CMOV<uint32_t>(const uint64_t& cond, uint32_t& val1, const uint32_t& val2) {
  CMOV4(cond, val1, val2);
}

template<>
INLINE void CMOV<uint16_t>(const uint64_t& cond, uint16_t& val1, const uint16_t& val2) {
  CMOV2(cond, val1, val2);
}

template<>
INLINE void CMOV<uint8_t>(const uint64_t& cond, uint8_t& val1, const uint8_t& val2) {
  CMOV1(cond, val1, val2);
}

template<>
INLINE void CMOV<bool>(const uint64_t& cond, bool& val1, const bool& val2) {
  CMOV_BOOL(cond, val1, val2);
}

template<>
INLINE void CMOV<int>(const uint64_t& cond, int& val1, const int& val2) {
  // UNDONE(): Make this a reinterpret cast?
  //
  CMOV4(cond, (uint32_t&) val1, val2);
}

// Other cmov specializations are inside the respective headers for the types.
// Make sure to always include this file first so that overloads get resolved,
// correctly.
//
// Aditionally, TSET and CXCHG need to be overloaded if CMOV is overloaded
// for a given type also:
//

#define OVERLOAD_TSET_CXCHG(TYPE, ...) \
 \
template<__VA_ARGS__> \
INLINE void TSET(bool selector, TYPE& A, const TYPE& B, const TYPE& C) { \
  CMOV(selector, A, C); \
  CMOV(!selector, A, B); \
} \
 \
template<__VA_ARGS__> \
INLINE void CTSET(bool condition, bool selector, TYPE& A, const TYPE& B, const TYPE& C) { \
  CMOV(condition*selector, A, C); \
  CMOV(condition*!selector, A, B); \
} \
\
template<__VA_ARGS__> \
INLINE void CXCHG(const uint64_t& cond, TYPE& A, TYPE& B) { \
  const TYPE C = A; \
  CMOV(cond, A, B); \
  CMOV(cond, B, C); \
} \

// Overload for CMOV of pair
template<typename A, typename B>
INLINE void CMOV(const uint64_t& cond, std::pair<A,B>& val1, const std::pair<A,B>& val2) {
  CMOV(cond, val1.first, val2.first);
  CMOV(cond, val1.second, val2.second);
}

OVERLOAD_TSET_CXCHG(std::pair<X COMMA Y>, typename X, typename Y);