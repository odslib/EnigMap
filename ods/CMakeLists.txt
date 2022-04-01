searchfiles()

include_directories (${OPENSSL_INCLUDE_DIR})
add_library (common common/encutils.cpp common/tracing/tracer.cpp)
target_compile_definitions(common PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_include_directories(common PUBLIC .)
target_link_libraries (common dl OpenSSL::SSL OpenSSL::Crypto)

add_library (oram INTERFACE)
target_include_directories(oram INTERFACE .)
target_link_libraries (oram INTERFACE common)

add_library (otree INTERFACE)
target_include_directories(otree INTERFACE .)
target_link_libraries (otree INTERFACE oram common)

add_executable (main main.cpp)
target_compile_definitions(main PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries (main dl otree)

# set_property(TARGET encutils PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
# set_property(TARGET oram PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
# set_property(TARGET tracer PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
# set_property(TARGET utils PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
# set_property(TARGET main PROPERTY INTERPROCEDURAL_OPTIMIZATION True)
