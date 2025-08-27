# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET hurchalla_factoring)
    return()
endif()

# later versions are probably fine, but are untested
cmake_minimum_required(VERSION 3.14...4.03)

project(hurchalla_factoring VERSION 1.0.0 LANGUAGES CXX)


# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()


# TODO: this section seems slightly messy for detecting/setting up testing
#
# if this is the top level CMakeLists.txt, add testing options, and enable
# testing when testing options have been set to ON.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(TEST_HURCHALLA_FACTORING
        "Build the tests for the Hurchalla modular arithmetic library project."
        ON)
    if(TEST_HURCHALLA_FACTORING)
        enable_testing()
        # include(CTest)
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
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/resource_intensive_api/FactorByTable32.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/resource_intensive_api/factorize_intensive_uint32.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/resource_intensive_api/is_prime_intensive.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/resource_intensive_api/IsPrimeTable.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/resource_intensive_api/is_prime_ultimate.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/DualMontgomeryForm.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/FactorizeStage2.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/factorize_trialdivision.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/factorize_wheel210.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/impl_factorize.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/impl_greatest_common_divisor.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/impl_is_prime.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/impl_is_prime_intensive.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/ImplFactorByTable32.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/is_prime_bruteforce.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/is_prime_miller_rabin.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/is_prime_trialdivision.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/IsPrimeFactor.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/microecm.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/OddPrimes.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PollardRhoBrentTrial.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PollardRhoBrentSwitchingTrial.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PollardRhoBrentTrialParallel.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PrimeTrialDivisionMayer.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/PrimeTrialDivisionWarren.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/SieveOfEratosthenes.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/trial_divide_mayer.h>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/detail/experimental/get_single_factor.h>
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
    GIT_TAG        9163344ee69f21a21cb8928dd30fc2ff15e94f0c
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
    if(TEST_HURCHALLA_FACTORING)
        add_subdirectory(test)
    endif()
endif()


# ***Benchmarking***
if(BENCH_HURCHALLA_FACTORING)
    add_subdirectory(benchmark)
endif()
