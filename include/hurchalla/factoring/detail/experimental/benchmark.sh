#!/bin/bash

# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

# Note that using CMake with -DCMAKE_BUILD_TYPE=Release will ensure that the
# standard macro NDEBUG (see <cassert>) is defined, which is essential for best
# performance.


cpp_standard="-std=c++17"
# pr_trial="-DHURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME=PollardRhoBrentTrialParallelAlt"
ecm_bits="-DHURCHALLA_FACTORING_ECM_THRESHOLD_BITS=70"
gcd_thresh="-DHURCHALLA_PRB_GCD_THRESHOLD=500"


build_dir=tmp
mkdir -p $build_dir
cmake -S. -B./$build_dir  \
        -DCMAKE_CXX_COMPILER=/opt/local/bin/clang++-mp-17 \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_CXX_FLAGS="$cpp_standard $pr_trial $ecm_bits $gcd_thresh"

cmake --build ./tmp --config Release

echo
echo Running example...
echo
./$build_dir/bench_hurchalla_factoring
