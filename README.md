# Oblivious Data Structure Library

A library providing external memory efficient, cpu instruction and memory access trace oblivious algorithms.

# Prerequisites:
Install cmake, ninja and intel sgx sdk, or use the cppbuilder docker image.

## How to build the builder docker image:
```bash
cd ./tools/docker/cppbuilder
docker build -t cppbuilder:latest .
```

## How to enter the docker environment to build under the builder:
```bash
docker run -it --rm -v $PWD:/builder -u $(id -u) cppbuilder
```

## An example for SGX
```bash
docker run --device=/dev/sgx_enclave -v /tmp/bucketfile:/ssdmount -it --rm -v $PWD:/builder cppbuilder
```

## How to enter the docker environment to build under the builder and run in SGX:
docker run --device=/dev/sgx_enclave -it --rm -v $PWD:/builder cppbuilder

## How to run the unit tests
```bash
rm -rf build
cmake -B build -G Ninja
ninja -C build
ninja -C build test
ninja -C build cppcheck
```

## How to do a release compilation:

```bash
rm -rf build # Needed after the CC/CXX export or after changing the CMAKE_BUILD_TYPE
export CC=/usr/bin/clang
export CXX=/usr/bin/clang++
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
ninja -C build
```


## Build the signal example enclave
```bash
source /startsgxenv.sh
cd applications/signal
make
```

## Performance graph on latest build

<!-- Commented out because due to not public yet:
### Map
![map perf graph](https://git.xtrm0.com/dsf/odsl/-/jobs/artifacts/main/raw/quality/graphs/perf.png?job=graphs)

![map perf graph](https://git.xtrm0.com/dsf/odsl/-/jobs/artifacts/main/raw/quality/graphs/perf-guidelines.png?job=graphs)


### External Memory Sorting
![sorting perf graph time](https://git.xtrm0.com/dsf/odsl/-/jobs/artifacts/main/raw/quality/graphs/sort_time.png?job=graphs)
![sorting perf graph page swaps](https://git.xtrm0.com/dsf/odsl/-/jobs/artifacts/main/raw/quality/graphs/sort_pages.png?job=graphs)
![sorting perf graph instructions](https://git.xtrm0.com/dsf/odsl/-/jobs/artifacts/main/raw/quality/graphs/sort_instr.png?job=graphs)

## Intrinsics performance graphs
![latest perf graph](https://git.xtrm0.com/dsf/odsl/-/jobs/artifacts/main/raw/quality/graphs/intrinsics.png?job=graphs)
-->

## Folder structure high level details

ods - C++ odsl library code
tests - C++ tests modules
applications - Enclaves using odsl
tools - tools used to generate graphs or test sets
tools/docker - dockerfiles used for reproducible builds

### Ods folder structure:

common - common c++ utilies, cpu abstractions, cryptography abstractions and tracing code
external_memory - external memory abstraction and sorting algorithms
external_memory/server - server abstraction for different external memory scenarios (sgx, file system, ram)
oram - oram implementations
otree - oblivious binary search tree
recoram - recursive oram implementations


## Links to view flamegraph files:

https://www.speedscope.app/

chrome://tracing/

It's also possible to view in qt creator as described here: https://doc.qt.io/qtcreator/creator-ctf-visualizer.html

At last: https://ui.perfetto.dev/ seems ok also



### Profiling code

1) Compile with ENABLE_PROFILING

2) For the functions that need profiling, add PROFILE_F(); at as the first line of the function code. Additionally add the function name to trace_events.hxx

3) Use PROFILER_SET(false); to disable profiling, use PROFILER_RESET() to write the profile to the log file (see profiling related functions in profiling_test to confirm).

4) Use any of the tools in "Links to view flamegraph files" above to look at the profiling, adjust uncached IO time based on the results of the benchmarks enclave (benchmark_sgx).

