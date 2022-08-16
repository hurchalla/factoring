This directory is used by the main factoring CMakeLists.txt (via add_subdirectory), in order to optionally facilitate benchmarking from github actions.  See the option BENCH_HURCHALLA_FACTORING within the main factoring CMakeLists.txt file.

However, you can also build and run the benchmark from this directory by executing ./benchmark.sh.

Benchmark.cpp contains a basic speed test of the API function factorize().  It is not a comparative benchmark with other libraries and so it may not be particularly interesting.  Though if you wish to peformance tune the factorize() function, it might be useful as a starting point.
