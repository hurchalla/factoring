# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET bench_hurchalla_factoring)
    return()
endif()

# later versions are probably fine, but are untested
cmake_minimum_required(VERSION 3.14...4.03)


#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/tests)
#set(CTEST_BINARY_DIRECTORY ${PROJECT_BINARY_DIR}/tests)


add_executable(bench_hurchalla_factoring)

target_sources(bench_hurchalla_factoring PRIVATE
               benchmark.cpp
               )


set_target_properties(bench_hurchalla_factoring
                      PROPERTIES FOLDER "Benchmarks")


if(NOT (TARGET hurchalla_factoring))
  # this CMakeLists.txt is usually invoked via add_subdirectory by the
  # factoring library's CMakeLists.txt.  In such a case, hurchalla_factoring
  # should be defined, and we would not enter this section.  However we can also
  # invoke this CMakeLists.txt through this directory's benchmark.sh.  In such a
  # case, hurchalla_factoring will not be defined, and we must add it.
  # In add_subdirectory below, we use ".." because (unless you moved it)
  # this CMakeLists.txt file has a path that is within the factoring repo, and
  # that path is one directory level below the repo's root.
  add_subdirectory(..  ${CMAKE_CURRENT_BINARY_DIR}/factoring)
endif()


target_link_libraries(bench_hurchalla_factoring
                      hurchalla_factoring)

if(MSVC)
  target_compile_options(bench_hurchalla_factoring PRIVATE /W4 /WX)
else()
  target_compile_options(bench_hurchalla_factoring PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
