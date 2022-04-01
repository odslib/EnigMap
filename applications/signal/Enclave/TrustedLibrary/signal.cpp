#include "../Enclave.h"
#include "Enclave_t.h"
#include "signal.hpp"


void ecall_signal_sanity()
{
    X_LOG("Testing bucket constructor");
    // volatile _ORAM::PathORAM::Bucket::Bucket<> b = _ORAM::PathORAM::Bucket::Bucket<>();
    X_LOG("Done Testing bucket constructor");
    Signal s(10000);
    s.RegisterUser(103213123);
    s.RegisterUser(1);
    s.RegisterUser(2);
    s.RegisterUser(3);
    s.RegisterUser(4);
    s.RegisterUser(6);
    X_LOG("Beginning test");
    for (uint64_t i=0; i<10'000; i++) {
        bool v1 = s.QueryUser(5);
        Assert(!v1);
        bool v2 = s.QueryUser(103213123);
        Assert(v2);
        // X_LOG(v1, v2);
    }
    X_LOG("Test done");
}

