#include "../Enclave.h"
#include "Enclave_t.h"


void ecall_cxx14_quoted()
{
    printf("Hello!\n");
}

void static_call_swap_page_4096(uint8_t* page_in, uint8_t* page_out) {
    uint8_t aux = 0x42;
    #pragma omp parallel for
    for (uint64_t i=0; i<4096; i++) {
        page_out[i] = page_in[i] ^ aux;
    }
}

void ecall_swap_page_4096(uint8_t* page_in, uint8_t* page_out)
{
    static_call_swap_page_4096(page_in, page_out);
}

void static_call_swap_page_512(uint8_t* page_in, uint8_t* page_out) {
    uint8_t aux = 0x42;
    #pragma omp parallel for
    for (uint64_t i=0; i<512; i++) {
        page_out[i] = page_in[i] ^ aux;
    }
}

void ecall_swap_page_512(uint8_t* page_in, uint8_t* page_out)
{
    static_call_swap_page_512(page_in, page_out);
}

void static_call_swap_page_32(uint8_t* page_in, uint8_t* page_out) {
    uint8_t aux = 0x42;
    #pragma omp parallel for
    for (uint64_t i=0; i<32; i++) {
        page_out[i] = page_in[i] ^ aux;
    }
}

void ecall_swap_page_32(uint8_t* page_in, uint8_t* page_out)
{
    static_call_swap_page_32(page_in, page_out);
}