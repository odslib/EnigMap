#include "../Enclave.h"
#include "Enclave_t.h"
#include "signal.hpp"


void ecall_signal_search(uint64_t sz)
{
    uint64_t t0; ocall_measure_time(&t0);
    Signal s(sz);
    uint64_t t1; ocall_measure_time(&t1);
    s.RegisterUser(103213123);
    s.RegisterUser(1);
    s.RegisterUser(2);
    s.RegisterUser(3);
    s.RegisterUser(4);
    s.RegisterUser(6);
    printf("Beginning test");
    uint64_t t2; ocall_measure_time(&t2);
    for (uint64_t i=0; i<10'000; i++) {
        bool v1 = s.QueryUser(5);
        Assert(!v1);
        bool v2 = s.QueryUser(103213123);
        Assert(v2);
        // printf(v1, v2);
    }
    uint64_t t3; ocall_measure_time(&t3);
    printf("Time to init: %lld\n", t1-t0);
    printf("Time to test: %lld\n", t3-t2);
    printf("Average time: %lld\n", (t3-t2)/20'000);
    printf("Test done");
}

void ecall_signal_insert(uint64_t sz)
{
    uint64_t t0; ocall_measure_time(&t0);
    Signal s(sz);
    uint64_t t1; ocall_measure_time(&t1);
    s.RegisterUser(103213123);
    s.RegisterUser(1);
    s.RegisterUser(2);
    s.RegisterUser(3);
    s.RegisterUser(4);
    s.RegisterUser(6);
    printf("Beginning test");
    uint64_t t2; ocall_measure_time(&t2);
    for (uint64_t i=0; i<10'000; i++) {
        s.RegisterUser(2*i+7);
        s.RegisterUser(2*i+1+7);
    }
    uint64_t t3; ocall_measure_time(&t3);
    printf("Time to init: %lld\n", t1-t0);
    printf("Time to test: %lld\n", t3-t2);
    printf("Average time: %lld\n", (t3-t2)/20'000);
    printf("Test done");
}


using EMPoint_t = typename std::pair<_OBST::K,_OBST::V>;
using Vector_t = typename EM::Vector::Vector<EMPoint_t>;

void ecall_signal_initialization(uint64_t sz)
{
    uint64_t t0; ocall_measure_time(&t0);
    Vector_t v(sz);
    OBST_t* client = new OBST_t(sz, v);
    uint64_t t1; ocall_measure_time(&t1);
    printf("Time to init: %lld\n", t1-t0);
}