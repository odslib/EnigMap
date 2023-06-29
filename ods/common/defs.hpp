#pragma once
#include "common/cpp_extended.hpp"

// ORAM bucket sizes:
// For PathORAM bucket size = ORAM__Z, for RingORAM = ORAM__Z + ORAM__S
#define ORAM__Z 4
#define ORAM__S 2

// Parameter A from Ring ORAM (number of accesses between reshufles)
#define ORAM__A 3

// ORAM doubly oblivious stash size:
#define ORAM__MS 80

// At what levels should encryption be used:
#define ORAM_SERVER__ENCRYPT_LARGE_BUCKETS false
#define ORAM__ENCRYPT_BUCKETS false
#define ORAM__ENCRYPT_BLOCKS false

// Info regarding emdas boas packing:
#ifndef ORAM_SERVER__DIRECTLY_CACHED_LEVELS
#define ORAM_SERVER__DIRECTLY_CACHED_LEVELS 16
#endif
#ifndef ORAM_SERVER__BUCKET_CACHE_SIZE
#define ORAM_SERVER__BUCKET_CACHE_SIZE 8192
#endif
#ifndef ORAM_SERVER__CACHE_SIZE
#define ORAM_SERVER__CACHE_SIZE 8192
#endif
#ifndef ORAM_SERVER__LEVELS_PER_PACK
#define ORAM_SERVER__LEVELS_PER_PACK 4
#endif

#define ORAM_USE_INRAM_SERVER true

#ifndef ENCLAVE_SIZE
#define ENCLAVE_SIZE 128
#endif
#ifndef DEFAULT_HEAP_SIZE
#define DEFAULT_HEAP_SIZE ((uint64_t)ENCLAVE_SIZE * 0xE0000UL)
#endif

// ORAM block size:
// #define BLOCK_SIZE 4096
// ORAM plaintext size (4096 - 16 * 2)
// #define PT_SIZE 4064
// ORAM data size (4096 - 16 * 2 - 8 * 2)
// #define DATA_SIZE 4048
// Multiplier for stash limit
#define STASH_MULTIPLIER 8

#define FLAMEGRAPHS_BASE_FOLDER "./quality/flamegraphs/"
// #define FLAMEGRAPHS_BASE_FOLDER "./"