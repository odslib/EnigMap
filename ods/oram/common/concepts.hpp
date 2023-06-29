#pragma once
#include <cinttypes>

#ifndef ENCLAVE_MODE
#include <concepts>
#endif

namespace _ORAM::Concepts {
  typedef uint64_t IndexType;

  template<template <typename T, uint64_t CACHE_SIZE> class FS
    , typename T
    , uint64_t CACHE_SIZE>
  concept FileServer = requires (FS<T, CACHE_SIZE> fs, T& t, const T& ct, IndexType indexType) {
    { FS<T, CACHE_SIZE>(indexType) };
    #ifndef ENCLAVE_MODE
    // concepts library is not supported:
    //
    { fs.Access(indexType) } -> std::same_as<T&>; 
    #endif
    { fs.Write(indexType, ct) };
    { fs.Read(indexType, t) };
  };  

} // namespace _ORAM::Concepts