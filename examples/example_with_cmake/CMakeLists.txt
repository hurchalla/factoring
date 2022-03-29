# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


cmake_minimum_required(VERSION 3.14)

project(factoring_example VERSION 1.0.0 LANGUAGES CXX)

add_executable(factoring_example)

target_sources(factoring_example PRIVATE
    example.cpp
    )

# in add_subdirectory below, the first argument must be the path on your system
# to the root of the cloned factoring repository.  In this example, "../.." is
# used because (unless you move it) this CMakeLists.txt file has a path that is
# within the repository, and that path is two directory levels below the repo's
# root.
add_subdirectory(../..  ${CMAKE_CURRENT_BINARY_DIR}/factoring)

target_link_libraries(factoring_example
                      hurchalla_factoring)


# To build this example, see the file example.sh.
