#!/bin/bash

# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Note that using CMake with -DCMAKE_BUILD_TYPE=Release will ensure that the
# standard macro NDEBUG (see <cassert>) is defined, which is essential for best
# performance.


cpp_standard="-std=c++17"

build_dir=tmp
mkdir -p $build_dir
cmake -S. -B./$build_dir  \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS="$cpp_standard"  \

cmake --build ./tmp --config Release

echo
echo Running example...
echo
./$build_dir/bench_hurchalla_factoring
