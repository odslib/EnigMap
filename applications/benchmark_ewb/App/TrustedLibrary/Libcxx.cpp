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
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000000; i++) {
            ecall_swap_page_4096(global_eid, page);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "4096: " << endValue_.count()/1000000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    {
        uint8_t page[512];
        uint8_t page2[512];
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000000; i++) {
            ecall_swap_page_512(global_eid, page);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "512: " << endValue_.count()/1000000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    {
        uint8_t page[32];
        uint8_t page2[32];
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000000; i++) {
            ecall_swap_page_32(global_eid, page);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "32: " << endValue_.count()/1000000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

        {
        uint8_t page[4096];
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000000; i++) {
            ecall_nops_4096(global_eid, page);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "nops-4096: " << endValue_.count()/1000000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    {
        uint8_t page[512];
        uint8_t page2[512];
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000000; i++) {
            ecall_nops_512(global_eid, page);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "nops-512: " << endValue_.count()/1000000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    {
        uint8_t page[32];
        uint8_t page2[32];
        auto start = high_resolution_clock::now();    
        for(int i=0; i<1000000; i++) {
            ecall_nops_32(global_eid, page);
        }
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "nops-32: " << endValue_.count()/1000000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    
    {
        uint8_t page[4096];
        auto start = high_resolution_clock::now();
        ecall_movs_1(global_eid, page);
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "mov-1: " << endValue_.count()/1'000'000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    {
        uint8_t page[4096];
        auto start = high_resolution_clock::now();
        ecall_movs_8(global_eid, page);
        auto stop = high_resolution_clock::now();
        std::chrono::nanoseconds endValue_ = stop - start;
        std::cout << "[Results]" << std::endl;
        std::cout << "mov-8: " << endValue_.count()/1'000'000 << "ns " << std::endl;
        std::cout << "[/Results]" << std::endl;
    }

    // 3358 * 8 = 26864
    // 7995 / (1-28/64) = 14213
    // 10466 / (1-28/128) = 13396.48
    // 11867 / (1-28/256) = 13324.35
    //
    {
        for (uint64_t times : {1,2,3}) {
            for (uint64_t stride : {512,1024,2048,4096,8192,2*8192,4*8192}) {
                auto start = high_resolution_clock::now();
                ecall_bm_ewb(global_eid, times, stride);
                auto stop = high_resolution_clock::now();
                std::chrono::nanoseconds endValue_ = stop - start;
                std::cout << "[Results]" << std::endl;
                std::cout << "bm_ewb-" << times << "-" << stride << ": " << endValue_.count()/(times*(1<<20)) << "ns " << std::endl;
                // std::cout << "bm_ewb-" << times << "-" << stride << ": " << endValue_.count()/(times*((1ULL<<28)/stride)) << "ns " << std::endl;
                std::cout << "[/Results]" << std::endl;
            }
        }
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

