#pragma once
#include <cinttypes>
#ifndef ENCLAVE_MODE
#warning Defining enclave mode
#define ENCLAVE_MODE
// #define NDEBUG
#endif

#include "common/cpp_extended.hpp"
#include <fstream>


#define DEFAULT_FILENAME "./storage.bin"

static uint64_t N;
static std::unique_ptr<std::fstream> lbios(nullptr);

void ocall_InitServer(uint64_t sizeOfT, uint64_t N_) {
  lbios.reset(new std::fstream());
  N = N_;
  lbios->open(DEFAULT_FILENAME, std::fstream::in | std::fstream::out | std::fstream::binary | std::fstream::trunc);
  Assert(lbios->is_open());
  for (uint64_t i=0; i<N; i++) {
    uint8_t zeros[sizeOfT];
    memset(zeros, 0, sizeOfT);
    lbios->write((char*)zeros, sizeOfT);
  }
  Assert(lbios->is_open());
  lbios->flush();
  X_LOG("ORAMClientInterface: Done allocating file (", lbios->tellp(), " bytes written)");
}

void ocall_SwapPage(uint64_t sizeOfT, uint64_t index_in, uint64_t index_out, uint8_t* page_in,  uint8_t* page_out) {
  Assert(index_in < N);
  Assert(index_out < N);
  std::streampos filePos = index_out * sizeOfT;
  lbios->seekg(filePos);
  lbios->read((char*)&page_out, sizeOfT);
  filePos = index_in * sizeOfT;
  lbios->seekp(filePos);
  lbios->write((char*)&page_in, sizeOfT);
}

void ocall_ReadPage(uint64_t sizeOfT, uint64_t i, uint8_t* page) {
  Assert(i < N);
  std::streampos filePos = i * sizeOfT;
  lbios->seekg(filePos);
  lbios->read((char*)&page, sizeOfT);
}

void ocall_WritePage(uint64_t sizeOfT, uint64_t i, uint8_t* page) {
  Assert(i < N);
  std::streampos filePos = i * sizeOfT;
  lbios->seekp(filePos);
  lbios->write((char*)&page, sizeOfT);
}