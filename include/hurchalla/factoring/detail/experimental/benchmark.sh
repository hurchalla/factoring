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


cpp_standard="-std=c++17"
pr_trial="-DHURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME=PollardRhoBrentTrialParallelOpt"
ecm_bits="-DHURCHALLA_FACTORING_ECM_THRESHOLD_BITS=50"
# 34 or 40 depending on whether arbitary size factors

gcd_thresh="-DHURCHALLA_PRB_GCD_THRESHOLD=100"
prb_startlen="-DHURCHALLA_PRB_STARTING_LENGTH=25"

td_size="-DHURCHALLA_TRIAL_DIVISION_SIZE=139"

use_lehman="-DHURCHALLA_USE_EXPERIMENTAL_LEHMAN"



#use_inline_asm_redc="-DHURCHALLA_ALLOW_INLINE_ASM_REDC=1"
use_all_inline_asm="-DHURCHALLA_ALLOW_INLINE_ASM_ALL=1"
#use_inline_asm_add="-DHURCHALLA_ALLOW_INLINE_ASM_MODADD"
#use_inline_asm_sub="-DHURCHALLA_ALLOW_INLINE_ASM_MODSUB"
#use_alt_hr_addsubs="-DHURCHALLA_MONTYHALFRANGE_USE_ALT_ADDSUBS"

microecm_expect_large_factors="-DHURCHALLA_FACTORING_EXPECT_LARGE_FACTORS"
#allow_microecm_dual_monty="-DHURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM"



# SET THIS TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT REPOSITORIES.
# (or otherwise ensure the compiler /I flags correctly specify the needed
# hurchalla include directories)
repo_directory=/home/jeff/Desktop

cppcompiler=clang++

source_file=benchmark.cpp


$cppcompiler  \
        -O3  -DNDEBUG \
        $cpp_standard \
        -I${repo_directory}/factoring/include \
        -I${repo_directory}/modular_arithmetic/modular_arithmetic/include \
        -I${repo_directory}/modular_arithmetic/montgomery_arithmetic/include \
        -I${repo_directory}/util/include \
        $microecm_expect_large_factors   $allow_microecm_dual_monty  \
        $use_inline_asm_add  $use_inline_asm_sub  $use_alt_hr_addsubs  \
        $use_inline_asm_redc  $use_all_inline_asm  \
        \
        $pr_trial $ecm_bits $gcd_thresh $prb_startlen \
        $td_size \
        $use_lehman \
        \
        $source_file -o benchmark

exit_on_failure

echo "compilation finished, now executing"


./benchmark
