#!/bin/bash

# Copyright (c) 2020-2022 Jeffrey Hurchalla.
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


exit_on_failure () {
  if [ $? -ne 0 ]; then
    exit 1
  fi
}


#allow_microecm_dual_monty="-DHURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM"

expect_arbitrary_factors="-DEXPECT_ARBITRARY_SIZE_FACTORS"

#comment out the following macro, in order to test the C++ version of microecm
use_c_interface="-DUSE_ECM_C_INTERFACE"


ccompiler=clang
cppcompiler=clang++
cpp_standard="-std=c++17"


# SET THIS TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT REPOSITORIES.
# (or otherwise ensure the compiler /I flags correctly specify the needed
# hurchalla include directories)
repo_directory=/home/jeff/Desktop



$ccompiler -O3  -DNDEBUG  -fomit-frame-pointer -march=native  -c microecm_c.c


$cppcompiler  \
        -O3  -DNDEBUG \
        $cpp_standard \
        -I${repo_directory}/factoring/include \
        -I${repo_directory}/modular_arithmetic/modular_arithmetic/include \
        -I${repo_directory}/modular_arithmetic/montgomery_arithmetic/include \
        -I${repo_directory}/util/include \
        $expect_arbitrary_factors   $allow_microecm_dual_monty  \
        $use_c_interface  \
        -c test_single_factor.cpp

$cppcompiler  -O3  -std="c++17"  -o test_single_factor_c  test_single_factor.o microecm_c.o -lm



exit_on_failure

echo "compilation finished, now executing"


./test_single_factor_c
