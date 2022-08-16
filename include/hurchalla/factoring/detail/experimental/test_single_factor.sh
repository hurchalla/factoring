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


# --------
# Various macros for performance experimentation.
# For details, see the file  ../../../../../../macros_for_performance.md

ecm_threshold="-DHURCHALLA_FACTORING_ECM_THRESHOLD_BITS=40"

trialdiv_size_crossover_bits="-DHURCHALLA_TRIAL_DIVISION_CROSSOVER_BITS=45"
trialdiv_size_small="-DHURCHALLA_TRIAL_DIVISION_SIZE_SMALL=109"
trialdiv_size_large="-DHURCHALLA_TRIAL_DIVISION_SIZE_LARGE=139"

prbst_gcd_threshold="-DHURCHALLA_PRBST_GCD_THRESHOLD=711"
prbst_starting_length="-DHURCHALLA_PRBST_STARTING_LENGTH=19"

#trial_type=-DHURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME=PollardRhoBrentTrial

#use_inline_asm_redc="-DHURCHALLA_ALLOW_INLINE_ASM_REDC=1"
use_all_inline_asm="-DHURCHALLA_ALLOW_INLINE_ASM_ALL=1" 
#use_inline_asm_add="-DHURCHALLA_ALLOW_INLINE_ASM_MODADD"
#use_inline_asm_sub="-DHURCHALLA_ALLOW_INLINE_ASM_MODSUB"

#use_alt_hr_addsubs="-DHURCHALLA_MONTYHALFRANGE_USE_ALT_ADDSUBS"

#microecm_expect_large_factors="-DHURCHALLA_FACTORING_EXPECT_LARGE_FACTORS"

#allow_microecm_dual_monty="-DHURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM"

# -----------



# SET THIS TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT REPOSITORIES.
# (or otherwise ensure the compiler -I flags correctly specify the needed
# hurchalla include directories)
repo_directory=/Users/jean/Downloads


cppcompiler=clang++
cpp_standard="-std=c++1z"


$cppcompiler  \
        -O3  -DNDEBUG  -march=haswell -mtune=haswell \
        $cpp_standard \
        -I${repo_directory}/factoring/include \
        -I${repo_directory}/modular_arithmetic/modular_arithmetic/include \
        -I${repo_directory}/modular_arithmetic/montgomery_arithmetic/include \
        -I${repo_directory}/util/include \
        $ecm_threshold  \
        $trialdiv_size_small  $trialdiv_size_large  $trialdiv_size_crossover_bits \
        $prbst_gcd_threshold  $prbst_starting_length \
        $microecm_expect_large_factors   $allow_microecm_dual_monty  \
        $cpp_standard  $use_inline_asm_add  $use_inline_asm_sub  $use_alt_hr_addsubs  \
        $use_inline_asm_redc  $use_all_inline_asm  $trial_type \
        test_single_factor.cpp -o test_single_factor

exit_on_failure

echo "compilation finished, now executing"


./test_single_factor



