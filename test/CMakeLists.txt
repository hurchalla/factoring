# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET test_hurchalla_factoring)
    return()
endif()

cmake_minimum_required(VERSION 3.14)


include(FetchGoogleTest.cmake)
include(EnableMaxWarnings.cmake)


#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/tests)
#set(CTEST_BINARY_DIRECTORY ${PROJECT_BINARY_DIR}/tests)


# needed for gtest_discover_tests()
include(GoogleTest)


add_executable(test_hurchalla_factoring)

target_sources(test_hurchalla_factoring PRIVATE
               factorize_bruteforce.h
               integer_sqrt.h
               sample_primes_and_nonprimes.h
               test_FactorByTable32.cpp
               test_factorize.cpp
               test_factorize_intensive_uint32.cpp
               test_greatest_common_divisor.cpp
               test_greatest_common_divisor_euclid.cpp
               test_is_prime.cpp
               test_IsPrimeIntensive.cpp
               detail_tests/is_prime_wheel210.h
               detail_tests/test_FactorizeStage2.cpp
               detail_tests/test_factorize_trialdivision.cpp
               detail_tests/test_factorize_wheel210.cpp
               detail_tests/test_is_prime_miller_rabin.cpp
               detail_tests/test_is_prime_trialdivision.cpp
               detail_tests/test_PrimeTrialDivisionWarrenAndMayer.cpp
               detail_tests/test_SieveOfEratosthenes.cpp
               detail_tests/test_trial_divide_mayer.cpp
               detail_tests/test_trial_divide_mayer_standard_division.cpp
               )

EnableMaxWarnings(test_hurchalla_factoring)


set_target_properties(test_hurchalla_factoring
                      PROPERTIES FOLDER "Tests")

target_link_libraries(test_hurchalla_factoring
                      hurchalla_factoring
                      gtest_main)

#add_test(test_hurchalla_factoring  test_hurchalla_factoring)
gtest_discover_tests(test_hurchalla_factoring)
