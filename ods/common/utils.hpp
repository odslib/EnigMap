#pragma once

#include <inttypes.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include "common/encutils.hpp"
#include "common/mov_intrinsics.hpp"
#include "common/defs.hpp"

constexpr INLINE uint64_t GetNextPowerOfTwo(uint64_t n)  {
  Assert(n <= 0x8000'0000'0000'0000);
  n = n-1;
  n |= n >> 1;
  n |= n >> 2;
  n |= n >> 4;
  n |= n >> 8;
  n |= n >> 16;
  n |= n >> 32;
  n = n + 1;
  return n;
}

INLINE uint64_t GetLogBaseTwo(uint64_t n) {
  const uint64_t masks[6] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF'0000, 0xFFFF'FFFF'0000'0000};
  uint64_t c = 32;
  uint64_t r = 0;
  for (int32_t i=5; i>=0; i--) {
    const bool cond = n & masks[i];
    CMOV(cond, n, n >> c);
    CMOV(cond, r, r | c);
    c >>= 1;
  }
  return r;
}

INLINE uint64_t CeilLog2(uint64_t x)
{
  static const uint64_t t[6] = {
    0xFFFFFFFF00000000ull,
    0x00000000FFFF0000ull,
    0x000000000000FF00ull,
    0x00000000000000F0ull,
    0x000000000000000Cull,
    0x0000000000000002ull
  };

  uint64_t y = (((x & (x - 1)) == 0) ? 0 : 1);
  int j = 32;
  int i;

  for (i = 0; i < 6; i++) {
    int k = (((x & t[i]) == 0) ? 0 : j);
    y += k;
    x >>= k;
    j >>= 1;
  }

  return y;
}


INLINE uint64_t reverseLowest32Bits(uint64_t x) {
    uint64_t cells;

    const uint64_t filter16 = 0x0000ffff0000fffful;
    const uint64_t filter8  = 0x00ff00ff00ff00fful;
    const uint64_t filter4  = 0x0f0f0f0f0f0f0f0ful;
    const uint64_t filter2  = 0x3333333333333333ul;
    const uint64_t filter1  = 0x5555555555555555ul;

    // swap neighboring bits in bundle of 16 and shift left by 16 bits
    cells = x & filter16;
    x = (x ^ cells) | (cells << 32);
    // swap neighboring bits in bundle of 8 and shift left by 8 bits
    cells = x & filter8;
    x = (x ^ cells) | (cells << 16);
    // swap neighboring bits in bundle of 4 and shift left by 4 bits
    cells = x & filter4;
    x = (x ^ cells) | (cells << 8);
    // swap neighboring bits in bundle of 2 and shift left by 2 bits
    cells = x & filter2;
    x = (x ^ cells) | (cells << 4);
    // swap neighboring bits in bundle of 1 and shift left by 1 bit
    cells = x & filter1;
    x = (x ^ cells) | (cells << 2);

    return x >> 31;

}


// Turns a number in [0, (1<<bits)-1] to Reverse Lexicographical Order
// in [0, (1<<bits)-1].
INLINE uint64_t ToReverseLexicographicalOrder(uint64_t n, uint64_t bits) {
  uint64_t reversedNumber = 0;
  for (int i=0; i<(int)bits; i++) {
    reversedNumber |= (((n >> i) & 1) << (bits - i - 1));
  }
  return ((1<<bits) - reversedNumber) - 1;
}

INLINE void GetRand16(uint8_t* out) {
  // UNDONE: we can make this random again,
  // but we don't actually need it for aes-gcm,
  // as long as iv's are not repeated.
  // we should not use /dev/urandom, because it is too
  // slow.
  //
  static uint64_t counter = 0;
  std::fill(out, out+16, 0);
  *(uint64_t*)out = counter;
  counter += 1;
  // static std::ifstream f("/dev/urandom");
	// f.read((char*)out, 16);
	// f.close();
}

// calculates mod(N, p) where N is expressed as a vector of uint64_t, i.e., base 2^64
// @pre: p should not exceed INT32_MAX
inline uint64_t large_num_mod(const std::vector<uint64_t>& nums, uint64_t p) {
  uint64_t res = 0;
  const uint64_t pow2_64_mod_p = (UINT64_MAX % p + 1) % p;
  for (uint64_t num: nums) {
    res = (res * pow2_64_mod_p + num % p) % p;
  }
  return res;
}

extern RandGen default_rand;

// [left,right]
inline uint64_t UniformRandom(uint64_t left, uint64_t right) {
  return default_rand.rand64() % (right - left + 1) + left;
}

inline uint32_t UniformRandom32(uint32_t left, uint32_t right) {
  return default_rand.rand32() % (right - left + 1) + left;
}

inline bool UniformRandomBit() {
  return default_rand.rand1();
}

// [0,right]
INLINE uint64_t UniformRandom(uint64_t right) {
  return UniformRandom(0,right);
}

// [0,right]
INLINE uint64_t UniformRandom() {
  return default_rand.rand64();
}

INLINE uint64_t UniformRandom32(uint32_t right) {
  return UniformRandom32(0,right);
}

// [0,right]
INLINE uint64_t UniformRandom32() {
  return default_rand.rand32();
}

// x/y round up
INLINE uint64_t divRoundUp(size_t x, size_t y) {
  return (x+y-1) / y;
}

/**
  * Note: this function is not oblivious
*/
template<typename Iterator>
void fisherYatesShuffle(Iterator begin, Iterator end) {
    size_t N = end - begin;
    for (size_t n = N - 1; n; --n) {
        size_t randPos = UniformRandom(n);
        std::swap(*(begin + randPos), *(--end));
    }
}