searchfiles()

include_directories (${PROJECT_SOURCE_DIR}/ods)
include_directories (${OPENSSL_INCLUDE_DIR})


add_executable(
  test_sanity
  indexers.cpp
  util_test.cpp
)
target_compile_definitions(test_sanity PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_sanity
  dl
  gtest_main
  oram
  common
  otree
  recoram
)

add_executable(
  test1
  test.cpp
  oram.cpp
  otree.cpp
  recursive.cpp
)
target_compile_definitions(test1 PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test1
  dl
  gtest_main
  oram
  common
  otree
  recoram
)

add_executable(
  test_algo
  algorithms.cpp
  sort.cpp
)
target_compile_definitions(test_algo PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_algo
  dl
  gtest_main
  common
)

add_executable(
  stagingtests
  otree_init.cpp
)
target_compile_definitions(stagingtests PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  stagingtests
  dl
  gtest_main
  oram
  common
  otree
  recoram
)



add_executable(
  test_basic_perf
  basic_perf.cpp
  perf_ods.cpp
  perf_sort.cpp
)
target_compile_definitions(test_basic_perf PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_basic_perf
  dl
  gtest_main
  oram
  common
  otree
)

add_executable(
  test_large_perf
  perf_ods_large.cpp
)
target_compile_definitions(test_basic_perf PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  test_large_perf
  dl
  gtest_main
  oram
  common
  otree
)


add_executable(
  profiling_test
  profiling_test.cpp
)
target_compile_definitions(profiling_test PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  profiling_test
  dl
  gtest_main
  oram
  common
  otree
)

add_executable(
  improvements_none
  improvements_none.cpp
)
target_compile_definitions(improvements_none PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  improvements_none 
  dl gtest_main oram common otree
)

add_executable(
  improvements_packing
  improvements_packing.cpp
)
target_compile_definitions(improvements_packing PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  improvements_packing 
  dl gtest_main oram common otree
)

add_executable(
  improvements_filecache
  improvements_filecache.cpp
)
target_compile_definitions(improvements_filecache PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  improvements_filecache 
  dl gtest_main oram common otree
)

add_executable(
  improvements_bucketcache
  improvements_bucketcache.cpp
)
target_compile_definitions(improvements_bucketcache PRIVATE BOOST_STACKTRACE_USE_ADDR2LINE)
target_link_libraries(
  improvements_bucketcache 
  dl gtest_main oram common otree
)

# set_property(TARGET test1 PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

include(GoogleTest)
gtest_discover_tests(test_sanity)
gtest_discover_tests(test1)
gtest_discover_tests(test_algo)
# gtest_discover_tests(basic_perf)
# gtest_discover_tests(stagingtests)
