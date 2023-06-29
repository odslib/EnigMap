#include "../Enclave.h"
#include "Enclave_t.h"

void ecall_cxx14_quoted()
{
    printf("Hello!\n");
}

void ecall_swap_page_4096(uint8_t* page_in)
{
    uint8_t aux = 0x42;
    for (uint64_t i=0; i<4096; i++) {
        page_in[i] = page_in[i] ^ aux;
    }
}

void ecall_swap_page_512(uint8_t* page_in)
{
    uint8_t aux = 0x42;
    for (uint64_t i=0; i<512; i++) {
        page_in[i] = page_in[i] ^ aux;
    }
}

void ecall_swap_page_32(uint8_t* page_in)
{
    uint8_t aux = 0x42;
    for (uint64_t i=0; i<32; i++) {
        page_in[i] = page_in[i] ^ aux;
    }
}

void ecall_movs_1(uint8_t* page_in)
{
    volatile uint64_t idx = 12345;
    for(uint64_t i=0; i<1'000'000; i++) {
        for (uint64_t j=0; j<4096; j++) {
            page_in[j] = idx;
        }
    }
}

void ecall_movs_8(uint8_t* page_in)
{
    volatile uint64_t idx = 12345;
    for(uint64_t i=0; i<1'000'000; i++) {
        for (uint64_t j=0; j<4096/8; j++) {
            ((uint64_t*)page_in)[j] = idx;
        }
    }
}

void ecall_nops_4096(uint8_t * page_in)
{
    return;
}

void ecall_nops_512(uint8_t * page_in)
{
    return;
}

void ecall_nops_32(uint8_t * page_in)
{
    return;
}

// ARR_SZ needs to be at least 2x the ammount of enclave memory so we ensure the pages
// are very likely not cached when we access them.
//
#define ARR_SZ (1ULL<<28)

volatile uint8_t* arr = nullptr;

static inline uint64_t fix_offset(uint64_t idx, uint64_t stride) {
    return (idx^(idx&(stride-1)));
}

void ecall_bm_ewb(uint64_t times, uint64_t stride)
{
    uint64_t prng = 123213;
    if (arr == nullptr) {
        arr = new uint8_t[ARR_SZ+0x10000];
        arr = (uint8_t*) (((uint64_t)(arr + 0x10000)) ^ (((uint64_t)arr) & 0xffff));
    }
    
    volatile uint8_t vol = 0x42;
    // for (uint64_t t=0; t<times; t++) {
    //     for (uint64_t i=0; i<(1<<24); i++) {
    //         prng = (prng*6364136223846793005ULL + 1442695040888963407ULL);
    //         uint64_t idx = fix_offset(prng%ARR_SZ, stride);
    //         arr[idx] = vol;
    //         arr[idx+stride] = vol;
    //     }
    // }
    uint64_t curr = 0;
    for (uint64_t t=0; t<times; t++) {
        for (uint64_t i=0; i<(1<<20); i++) {
            uint64_t idx = curr + (i%stride) ;
            curr = (curr + 2*stride) % ARR_SZ;
            arr[idx] = vol;
            arr[idx+stride] = vol;
        }
    }
}