#!/bin/bash

# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Note that using CMake with -DCMAKE_BUILD_TYPE=Release will ensure that the
# standard macro NDEBUG (see <cassert>) is defined, which is essential for best
# performance.


allow_ecm="-DHURCHALLA_FACTORING_ALLOW_ECM_EXPERIMENTAL"
#expect_large_factors="-DHURCHALLA_FACTORING_EXPECT_LARGE_FACTORS"

use_all_inline_asm="-DHURCHALLA_ALLOW_INLINE_ASM_ALL=1"


cpp_standard="-std=c++1z"

build_dir=tmp
mkdir -p $build_dir
cmake -S. -B./$build_dir  \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS="$allow_ecm  \
        $expect_large_factors  \
        $cpp_standard  \
        $use_all_inline_asm"  \

cmake --build ./tmp --config Release

echo
echo Running example...
echo
./$build_dir/bench_hurchalla_factoring
