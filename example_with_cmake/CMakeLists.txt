# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---


cmake_minimum_required(VERSION 3.14)

project(factoring_example VERSION 1.0.0 LANGUAGES CXX)

add_executable(factoring_example)

target_sources(factoring_example PRIVATE
    example.cpp
    )

add_subdirectory(..  ${CMAKE_CURRENT_BINARY_DIR}/factoring)

target_link_libraries(factoring_example
                      hurchalla_factoring)
