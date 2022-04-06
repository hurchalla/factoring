# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET hurchalla_factoring)
    return()
endif()

cmake_minimum_required(VERSION 3.14)

project(hurchalla_factoring VERSION 1.0.0 LANGUAGES CXX)


# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()


# TODO: this section seems slightly messy for detecting/setting up testing
# --------------------
option(TEST_HURCHALLA_LIBS
       "Build the tests for all Hurchalla library projects."
       OFF)

# if this is the top level CMakeLists.txt, add testing options, and enable
# testing when testing options have been set to ON.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(TEST_HURCHALLA_FACTORING
        "Build the tests for the Hurchalla modular arithmetic library project."
        ON)
    if(TEST_HURCHALLA_FACTORING OR TEST_HURCHALLA_LIBS)
        enable_testing()
        # include(CTest)
    endif()
elseif(TEST_HURCHALLA_LIBS)
    # If TEST_HURCHALLA_LIBS is set to ON, enable_testing() should have been
    # called either directly or indirectly by the top level project. (Note that
    # if a project calls include(CTest), the included CTest.cmake defines a
    # BUILT_TESTING option and calls enable_testing if BUILD_TESTING is ON.)
    if (NOT CMAKE_TESTING_ENABLED)
        message(FATAL_ERROR "Fatal error: TEST_HURCHALLA_LIBS was set, but enable_testing() was never called")
    endif()
endif()


# an option to facilitate benchmarking from github actions
option(BENCH_HURCHALLA_FACTORING
       "Build a benchmarking app for the factoring project."
       OFF)



if(NOT DEFINED CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
if(NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()



add_library(hurchalla_factoring INTERFACE)


target_sources(hurchalla_factoring INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/factorize.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/greatest_common_divisor.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/is_prime.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/resource_intensive_api/factorize_intensive_uint32.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/resource_intensive_api/IsPrimeIntensive.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/factorize_pollard_rho.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/factorize_trialdivision.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/factorize_wheel210.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/impl_factorize.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/impl_greatest_common_divisor.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/impl_is_prime.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/ImplIsPrimeIntensive.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/is_prime_bruteforce.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/is_prime_miller_rabin.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/is_prime_trialdivision.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/OddPrimes.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PollardRhoBrentTrial.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PollardRhoIsPrime.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PrimeTrialDivisionMayer.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PrimeTrialDivisionWarren.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/SieveOfEratosthenes.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/trial_divide_mayer.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/experimental/PollardRhoBrentMontgomeryTrial.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/experimental/PollardRhoBrentMontgomeryTrialParallel.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/experimental/PollardRhoBrentSwitchingTrial.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/experimental/PollardRhoTrial.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases16_1.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases16_2.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases30_1.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases30_2.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases31_1.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases31_2.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases32_1.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases32_2.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases32_3.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases62_2.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases62_3.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases62_4.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases62_5.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases62_6.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases63_2.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases63_3.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases63_4.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases63_5.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases63_6.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_2.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_3.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_4.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_5.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_6.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_7.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/MillerRabinProbabilisticBases128.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/alternative_tables/MillerRabinBases63_5_alt.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/miller_rabin_bases/alternative_tables/MillerRabinBases64_3_alt.h>
    )

install(DIRECTORY
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/>
    DESTINATION include)


target_include_directories(hurchalla_factoring INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)


target_compile_features(hurchalla_factoring INTERFACE cxx_std_17)


include(FetchContent)

FetchContent_Declare(
    hurchalla_util
    GIT_REPOSITORY https://github.com/hurchalla/util.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(hurchalla_util)
target_link_libraries(hurchalla_factoring
                      INTERFACE hurchalla_util)


FetchContent_Declare(
    modular_arithmetic
    GIT_REPOSITORY https://github.com/hurchalla/modular_arithmetic.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(modular_arithmetic)
target_link_libraries(hurchalla_factoring
                      INTERFACE hurchalla_modular_arithmetic)



# ***Tests***

# if this is the top level CMakeLists.txt
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    if(TEST_HURCHALLA_FACTORING OR TEST_HURCHALLA_LIBS)
        add_subdirectory(test)
    endif()
elseif(TEST_HURCHALLA_LIBS)
    add_subdirectory(test)
endif()


# ***Benchmarking***
if(BENCH_HURCHALLA_FACTORING)
    add_subdirectory(benchmark)
endif()
