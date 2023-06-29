#include <stdio.h>

#include "../App.h"
#include "Enclave_u.h"
#include <thread>
#ifdef DISK_IO
#include "external_memory/server/enclaveFileServer_untrusted.hpp"
#else
#include "external_memory/server/enclaveMemServer_untrusted.hpp"
#endif


void ActualMain(void)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    
    printf("Doing ecall\n");
    
    uint64_t TEST_SELECTOR = 0;

    for (uint64_t i=1<<8; i<1<<30; i*=2) {
        if (TEST_SELECTOR == 0) { // Search
            ret = ecall_signal_search(global_eid, i);
        } else if (TEST_SELECTOR == 1) { // Insert
            ret = ecall_signal_insert(global_eid, i);
        } else { // Initialization
            ret = ecall_signal_initialization(global_eid, i);
        }
        if (ret != SGX_SUCCESS)
            abort();
    }
    printf("Done\n");

    if (ret != SGX_SUCCESS)
        abort();    
}