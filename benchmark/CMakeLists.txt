# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


if(TARGET bench_hurchalla_factoring)
    return()
endif()

cmake_minimum_required(VERSION 3.14)


#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/tests)
#set(CTEST_BINARY_DIRECTORY ${PROJECT_BINARY_DIR}/tests)


add_executable(bench_hurchalla_factoring)

target_sources(bench_hurchalla_factoring PRIVATE
               benchmark.cpp
               )


set_target_properties(bench_hurchalla_factoring
                      PROPERTIES FOLDER "Benchmarks")

target_link_libraries(bench_hurchalla_factoring
                      hurchalla_factoring)

if(MSVC)
  target_compile_options(bench_hurchalla_factoring PRIVATE /W4 /WX)
else()
  target_compile_options(bench_hurchalla_factoring PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
