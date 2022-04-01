#include <stdio.h>

#include "../App.h"
#include "Enclave_u.h"
#include <thread>
#include "external_memory/server/enclaveFileServer_untrusted.hpp"


void ActualMain(void)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    
    printf("Doing ecall\n");

    ret = ecall_signal_sanity(global_eid);

    printf("Done\n");

    if (ret != SGX_SUCCESS)
        abort();
        
    printf("Did not abort\n");
}

