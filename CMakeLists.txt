cmake_minimum_required(VERSION 3.18)
project("Private Contact Discovery" LANGUAGES C CXX)

include(cmake/buildtype.cmake)
include(cmake/faster.cmake)

set(PROJECT_TRDPARTY_DIR "thirdparty")

# gtest is currently buggy with clang and interprocess optimization, we needded to remove this flag: -Werror
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -fno-inline -no-pie")
# This doesn't work if the machine doesn't support icelake-server march.
# I'm removing it for now and replacing with -march=native, -mtune=native
# UNDONE(): Make sure that on evaluation machine we validate that AVX512VL is being set.
#
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -march=native -mtune=native")
include(cmake/colors.cmake)

# Includes:
#
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/3ea587050da9447536d0b55fece0a240273d9927.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# If you don't want to install openssl, uncommenting the following lines should work:
# FetchContent_Declare(
#   ssl_content
#   GIT_REPOSITORY https://github.com/janbar/openssl-cmake.git
#   GIT_TAG "1.1.1k-20210430"
# )
# FetchContent_MakeAvailable(ssl_content)
# add_library(OpenSSL::SSL ALIAS ssl)
# add_library(OpenSSL::Crypto ALIAS crypto)

find_package(OpenSSL REQUIRED)
find_package (Boost REQUIRED)

include(cmake/searchfiles.cmake)
include(cmake/boost.cmake)



add_subdirectory(ods)


enable_testing()
add_subdirectory(tests)



include(cmake/cppcheck.cmake)
include(cmake/clangformat.cmake)
