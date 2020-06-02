# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET hurchalla_factoring)
    return()
endif()

cmake_minimum_required(VERSION 3.13)

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
        OFF)
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



# TODO: I'm not sure whether these lines provide any benefit:
# --------------------
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


# Check for required C++ features only if the C++ compiler ID is non-empty.
#target_compile_features(hurchalla_factoring INTERFACE 
#          $<$<NOT:$<CXX_COMPILER_ID:>>:
#               cxx_static_assert
#           >)


target_sources(hurchalla_factoring INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/hurchalla/factoring/is_prime_monty.h>
    )


target_include_directories(hurchalla_factoring INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)



# TODO: Currently we assume the programming_by_contract dependency is in a
# folder at the relative location specified below.  This is a hack, since users
# usually won't know this is a requirement, and it obligates the user to perform
# a manual step (even if only once) to get the dependency from git and place it
# in the correct location.  The dependency should be completely invisible to the
# user, and this cmake script should "just work" without any extra steps.
# Among the possible solutions might be git submodules, or FetchContent.  Git
# submodules often seems to be warned against, making me more optimistic about
# FetchContent or perhaps some other solution.  Note: the tests for this project
# currently use FetchContent to get the Google Tests dependency.
# ---------------------------
add_subdirectory(../programming_by_contract
            ${CMAKE_CURRENT_BINARY_DIR}/programming_by_contract)
target_link_libraries(hurchalla_factoring
                      INTERFACE programming_by_contract)

# TODO: address the same issue as above, but for modular_arithmetic dependency
# ---------------------------
add_subdirectory(../modular_arithmetic
            ${CMAKE_CURRENT_BINARY_DIR}/modular_arithmetic)
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
