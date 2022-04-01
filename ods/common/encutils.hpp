#pragma once
#include <cassert>
#include "common/defs.hpp"

#ifndef AES_BLOCK_SIZE
#define AES_BLOCK_SIZE 32
#endif
void aes_init();
void aes_256_gcm_encrypt(uint64_t plaintextSize, uint8_t* plaintext, const uint8_t key[AES_BLOCK_SIZE], uint8_t iv[AES_BLOCK_SIZE], uint8_t tag[AES_BLOCK_SIZE], uint8_t* ciphertext);
bool aes_256_gcm_decrypt(uint64_t ciphertextSize, uint8_t* ciphertext, const uint8_t key[AES_BLOCK_SIZE], uint8_t iv[AES_BLOCK_SIZE], uint8_t tag[AES_BLOCK_SIZE], uint8_t* plaintext);
