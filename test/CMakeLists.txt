# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET test_hurchalla_factoring)
    return()
endif()

cmake_minimum_required(VERSION 3.13)


include(FetchGoogleTest.cmake)
include(EnableMaxWarnings.cmake)


#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/tests)
#set(CTEST_BINARY_DIRECTORY ${PROJECT_BINARY_DIR}/tests)


# needed for gtest_discover_tests()
include(GoogleTest)


add_executable(test_hurchalla_factoring)

target_sources(test_hurchalla_factoring PRIVATE
               test_factoring.cpp
               is_prime_wheel210.h
               is_prime_bruteforce.h
               integer_sqrt.h
               )

EnableMaxWarnings(test_hurchalla_factoring)


set_target_properties(test_hurchalla_factoring
                      PROPERTIES FOLDER "Tests")

target_link_libraries(test_hurchalla_factoring
                      hurchalla_factoring
                      gtest_main)

#add_test(test_hurchalla_factoring  test_hurchalla_factoring)
gtest_discover_tests(test_hurchalla_factoring)
