#!/bin/bash

# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

exit_on_failure () {
  if [ $? -ne 0 ]; then
    exit 1
  fi
}


# using a macro to index the semiprimes file number we want to use is a hack...
# ...but it's easy and it works
semiprime_filenum="-DUECM_SEMIPRIME_FILENUM=13"

input_type="-DUECM_INPUT_TYPE=uint64_t"
#input_type="-DUECM_INPUT_TYPE=__uint128_t"

num_tries="-DUECM_NUM_TRIES=5"

start=0
length=100000



# this macro can increase the C++ version's speed at 63 and 64 bits, but otherwise slows things down
#allow_microecm_dual_monty="-DHURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM"


microecm_c_filename="microecm_c"

ccompiler=clang
cppcompiler=clang++
#ccompiler=gcc
#cppcompiler=g++
cpp_standard="-std=c++17"

# SET THIS TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT REPOSITORIES.
# (or otherwise ensure the compiler /I flags correctly specify the needed
# hurchalla include directories)
repo_directory=/home/jeff/Desktop

#yafu_include="-I/home/jeff/Desktop/yafu/include"



$ccompiler  -O3  -DNDEBUG  -fomit-frame-pointer -march=native  \
        -pthread $yafu_include  $semiprime_filenum  -c ${microecm_c_filename}.c
#yafu/arith/arith.c


$cppcompiler  \
        -O3  -DNDEBUG  -fomit-frame-pointer -march=native \
        $cpp_standard \
        $yafu_include \
        -I${repo_directory}/factoring/include \
        -I${repo_directory}/modular_arithmetic/modular_arithmetic/include \
        -I${repo_directory}/modular_arithmetic/montgomery_arithmetic/include \
        -I${repo_directory}/util/include \
        $semiprime_filenum  $allow_microecm_dual_monty  $num_tries  $input_type \
        -pthread -c ben_test.cpp

#$cppcompiler  -O3  -std="c++17"  -o ben_test   ben_test.o ${microecm_c_filename}.o arith.o -lm
$cppcompiler  -pthread -O3  -std="c++17"  -o ben_test   ben_test.o ${microecm_c_filename}.o  -lm


exit_on_failure

echo "compilation finished, now executing"


./ben_test  $start  $length
