#include <utility>


#ifndef NOOPENSSL
#include <openssl/aes.h>
#include <openssl/evp.h>
#include <openssl/rand.h>
#include <openssl/crypto.h>
#include "common/encutils.hpp"


void aes_init() {
  static int init = 0;
  if (init == 0) {
    OpenSSL_add_all_ciphers();
    // int rv = RAND_load_file("/dev/urandom", 32);
    init = 1;
  }
}

void aes_256_gcm_encrypt(uint64_t plaintextSize, uint8_t *plaintext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *ciphertext) {

  aes_init();

  size_t enc_length = ((plaintextSize + 15) / 16) * 16;

  int actual_size = 0, final_size = 0;
  EVP_CIPHER_CTX *e_ctx = EVP_CIPHER_CTX_new();
  EVP_CIPHER_CTX_ctrl(e_ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL);
  EVP_CIPHER_CTX_set_padding(e_ctx, 0);
  EVP_EncryptInit(e_ctx, EVP_aes_256_gcm(), key, iv);

  EVP_EncryptUpdate(e_ctx, ciphertext, &actual_size, plaintext, plaintextSize);
  int ok = EVP_EncryptFinal(e_ctx, &ciphertext[actual_size], &final_size);
  Assert(final_size <= enc_length);
  EVP_CIPHER_CTX_ctrl(e_ctx, EVP_CTRL_GCM_GET_TAG, 16, tag);
  EVP_CIPHER_CTX_free(e_ctx);
  Assert(ok == 1);
  IGNORE_UNUSED(enc_length);
  IGNORE_UNUSED(ok);
}

bool aes_256_gcm_decrypt(uint64_t ciphertextSize, uint8_t *ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *plaintext) {
  aes_init();

  // UNDONE(): Make sure this is using aesni
  //

  // Assert(ciphertextSize % 16 == 0);
  int actual_size = 0, final_size = 0;
  EVP_CIPHER_CTX *d_ctx = EVP_CIPHER_CTX_new();
  EVP_CIPHER_CTX_ctrl(d_ctx, EVP_CTRL_GCM_SET_IVLEN, 16, NULL);
  EVP_CIPHER_CTX_set_padding(d_ctx, 0);
  EVP_DecryptInit(d_ctx, EVP_aes_256_gcm(), key, iv);
  EVP_DecryptUpdate(d_ctx, plaintext, &actual_size, ciphertext, ciphertextSize);
  EVP_CIPHER_CTX_ctrl(d_ctx, EVP_CTRL_GCM_SET_TAG, 16, tag);
  int ok;
  ok = EVP_DecryptFinal(d_ctx, &plaintext[actual_size], &final_size);
  EVP_CIPHER_CTX_free(d_ctx);

  Assert(ok == 1);
  return ok == 1;
}

#else
#include "sodium.h"
#include "common/encutils.hpp"
#define MAX_ENC_SIZE 8192

static uint8_t* sodium_buffer = nullptr;
// In enclave mode we can't use openssl, so we use libsodium, which should be packed in the enclave directly.
void aes_init() {
  static int init = 0;
  if (init == 0) {
    Assert(crypto_aead_aes256gcm_is_available());
    sodium_buffer = sodium_malloc(MAX_ENC_SIZE + crypto_aead_aes256gcm_ABYTES);
    // OpenSSL_add_all_ciphers();
    // int rv = RAND_load_file("/dev/urandom", 32);
    init = 1;
  }
}

void aes_256_gcm_encrypt(uint64_t plaintextSize, uint8_t *plaintext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *ciphertext) {
  aes_init();

  memset(tag, 0, AES_BLOCK_SIZE);
  size_t enc_length = ((plaintextSize + 15) / 16) * 16;
  uint64_t reportedTagLen;
  crypto_aead_aes256gcm_encrypt_detached(
    ciphertext, 
    tag, &reportedTagLen,
    plaintext, plaintextSize,
    /*ad=*/nullptr, 0,
    /*nsec=*/nullptr,
    iv,
    key
  );
  Assert(reportedTagLen == crypto_aead_aes256gcm_ABYTES);
}

bool aes_256_gcm_decrypt(uint64_t ciphertextSize, uint8_t *ciphertext,
                         const uint8_t key[AES_BLOCK_SIZE],
                         uint8_t iv[AES_BLOCK_SIZE],
                         uint8_t tag[AES_BLOCK_SIZE], uint8_t *plaintext) {
  aes_init();

  int ok;
  ok = crypto_aead_aes256gcm_decrypt_detached(
    plaintext,
    /*nsec=*/nullptr,    
    ciphertext, ciphertextSize,
    tag,
    /*ad=*/nullptr, 0,
    iv,
    key    
  );
  Assert(ok==0);
  return ok == 0;
}


#endif