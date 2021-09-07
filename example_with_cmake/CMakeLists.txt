# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---


cmake_minimum_required(VERSION 3.14)

project(factoring_example VERSION 1.0.0 LANGUAGES CXX)

add_executable(factoring_example)

target_sources(factoring_example PRIVATE
    example.cpp
    )

# in add_subdirectory below, the first argument must be the path on your system
# to the root of the cloned factoring repository.  In this example, ".." is used
# because (unless you move it) this example's CMakeLists.txt has a path that is
# within the repository, and that path is one directory level below the repo's
# root.
add_subdirectory(..  ${CMAKE_CURRENT_BINARY_DIR}/factoring)

target_link_libraries(factoring_example
                      hurchalla_factoring)
