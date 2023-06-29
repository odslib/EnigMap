#include "../Enclave.h"
#include "Enclave_t.h"
#include "common/encutils.hpp"


void ecall_cxx14_quoted()
{
    printf("Hello!\n");
}

void ecall_swap_page_4096(uint8_t* page_in, uint8_t* page_out)
{
    uint8_t aux = 0x42;
    for (uint64_t i=0; i<4096; i++) {
        page_out[i] = page_in[i] ^ aux;
    }
}

void ecall_swap_page_512(uint8_t* page_in, uint8_t* page_out)
{
    uint8_t aux = 0x42;
    for (uint64_t i=0; i<512; i++) {
        page_out[i] = page_in[i] ^ aux;
    }
}

void ecall_swap_page_32(uint8_t* page_in, uint8_t* page_out)
{
    uint8_t aux = 0x42;
    for (uint64_t i=0; i<32; i++) {
        page_out[i] = page_in[i] ^ aux;
    }
}

void ecall_movs_1(uint8_t* page_in, uint8_t* page_out)
{
    volatile uint64_t idx = 12345;
    for(uint64_t i=0; i<1'000'000; i++) {
        for (uint64_t j=0; j<4096; j++) {
            page_out[j] = idx;
        }
    }
}

void ecall_movs_8(uint8_t* page_in, uint8_t* page_out)
{
    volatile uint64_t idx = 12345;
    for(uint64_t i=0; i<1'000'000; i++) {
        for (uint64_t j=0; j<4096/8; j++) {
            ((uint64_t*)page_out)[j] = idx;
        }
    }
}

void ecall_nops_4096(uint8_t * page_in, uint8_t * page_out)
{
    return;
}

void ecall_nops_512(uint8_t * page_in, uint8_t * page_out)
{
    return;
}

void ecall_nops_32(uint8_t * page_in, uint8_t * page_out)
{
    return;
}

// ARR_SZ needs to be at least 2x the ammount of enclave memory so we ensure the pages
// are very likely not cached when we access them.
//
#define ARR_SZ (1ULL<<32)

uint8_t* arr = nullptr;

void ecall_bm_ewb(uint64_t times)
{
    if (arr == nullptr) {
        arr = new uint8_t[ARR_SZ];
    }
    // UNDONE():
    //
    volatile uint8_t vol = 0x42;
    for (uint64_t t=0; t<times; t++) {
        for (uint64_t i=0; i<ARR_SZ; i++) {
            arr[i] = vol;
        }
    }
}

volatile uint8_t plaintext[1<<16] = {0};
volatile uint8_t ciphertext[1<<16] = {0};
uint8_t key[AES_BLOCK_SIZE];
uint8_t iv[AES_BLOCK_SIZE];
uint8_t tag[AES_BLOCK_SIZE];

void ecall_bm_encrypt(uint64_t size) {
    aes_init();
    
    // 
    // aes_256_gcm_decrypt(64+16, ciphertext, key, iv, tag, plaintext);

    for (uint64_t i=0; i<4'000'000; i++) {
        uint8_t* pt = (uint8_t*) plaintext;
        uint8_t* ct = (uint8_t*) ciphertext;
        aes_256_gcm_encrypt(size, pt, key, iv, tag, ct);
        // aes_256_ctr_encrypt(size, pt, key, iv, ct);
    }
}

void ecall_bm_decrypt(uint64_t size) {
    aes_init();
    
    {
        uint8_t* pt = (uint8_t*) plaintext;
        uint8_t* ct = (uint8_t*) ciphertext;
        aes_256_gcm_encrypt(size, pt, key, iv, tag, ct);
        // aes_256_ctr_encrypt(size, pt, key, iv, ct);
    }

    for (uint64_t i=0; i<4'000'000; i++) {
        uint8_t* pt = (uint8_t*) plaintext;
        uint8_t* ct = (uint8_t*) ciphertext;
        // aes_256_ctr_decrypt(size, pt, key, iv, ct);
        aes_256_gcm_decrypt(size+16, ct, key, iv, tag, pt);
    }    
}

