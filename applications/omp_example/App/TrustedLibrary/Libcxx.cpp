#include <stdio.h>
#include <iostream>
#include "../App.h"
#include "Enclave_u.h"
#include <thread>

#include <chrono>
using namespace std::chrono;


void BenchmarkCall(void)
{
    {
        uint8_t page[4096];
        uint8_t page2[4096];
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000; i++) {
            ecall_swap_page_4096(global_eid, page, page2);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "4096: " << endValue_.count()/1000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    {
        uint8_t page[512];
        uint8_t page2[512];
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000; i++) {
            ecall_swap_page_512(global_eid, page, page2);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "512: " << endValue_.count()/1000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    {
        uint8_t page[32];
        uint8_t page2[32];
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000; i++) {
            ecall_swap_page_32(global_eid, page, page2);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "32: " << endValue_.count()/1000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }
}

void ActualMain(void)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;

    BenchmarkCall();

    ret = ecall_cxx14_quoted(global_eid);
    if (ret != SGX_SUCCESS)
        abort();
}

