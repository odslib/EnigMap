#pragma once

#include <inttypes.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>

#include "common/mov_intrinsics.hpp"
#include "common/defs.hpp"

inline uint64_t GetNextPowerOfTwo(uint64_t n)  {
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

inline uint64_t GetLogBaseTwo(uint64_t n) {
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

inline uint64_t CeilLog2(uint64_t x)
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

// Turns a number in [0, (1<<bits)-1] to Reverse Lexicographical Order
// in [0, (1<<bits)-1].
inline uint64_t ToReverseLexicographicalOrder(uint64_t n, uint64_t bits) {
  uint64_t reversedNumber = 0;
  for (int i=0; i<bits; i++) {
    reversedNumber |= (((n >> i) & 1) << (bits - i - 1));
  }
  return ((1<<bits) - reversedNumber) - 1;
}

inline void GetRand16(uint8_t* out) {
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

// [left,right]
inline uint64_t UniformRandom(uint64_t left, uint64_t right) {
  // UNDONE(0): This needs to be random
  //
  static uint64_t cache = 0;
  cache = (cache + 1);
  return ((cache) % (right-left+1)) + left;
}

// [0,right]
inline uint64_t UniformRandom(uint64_t right) {
  return UniformRandom(0,right);
}