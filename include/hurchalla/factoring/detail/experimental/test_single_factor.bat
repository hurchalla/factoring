@echo off


set ecm_threshold="-DHURCHALLA_FACTORING_ECM_THRESHOLD_BITS=40"

set trialdiv_size="-DHURCHALLA_TRIAL_DIVISION_SIZE=139"

set prbst_gcd_threshold="-DHURCHALLA_PRBST_GCD_THRESHOLD=711"
set prbst_starting_length="-DHURCHALLA_PRBST_STARTING_LENGTH=19"

REM set trial_type=-DHURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME=PollardRhoBrentTrial

REM set use_inline_asm_redc="-DHURCHALLA_ALLOW_INLINE_ASM_REDC=1"
set use_all_inline_asm="-DHURCHALLA_ALLOW_INLINE_ASM_ALL=1"
REM set use_inline_asm_add="-DHURCHALLA_ALLOW_INLINE_ASM_MODADD"
REM set use_inline_asm_sub="-DHURCHALLA_ALLOW_INLINE_ASM_MODSUB"

REM set use_alt_hr_addsubs="-DHURCHALLA_MONTYHALFRANGE_USE_ALT_ADDSUBS"

REM set microecm_expect_large_factors="-DHURCHALLA_FACTORING_EXPECT_LARGE_FACTORS"

REM set allow_microecm_dual_monty="-DHURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM"



REM SET THIS TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT REPOSITORIES.
REM (or otherwise ensure the compiler /I flags correctly specify the needed
REM hurchalla include directories)
set repo_directory=C:\Users\Jeff\source\repos\github



set cppcompiler=cl
set cpp_standard="/std:c++17"


%cppcompiler%  ^
        /O2  /DNDEBUG  /EHsc ^
        %cpp_standard% ^
        /I%repo_directory%\util\include\ ^
        /I%repo_directory%\factoring\include\ ^
        /I%repo_directory%\modular_arithmetic\modular_arithmetic\include\ ^
        /I%repo_directory%\modular_arithmetic\montgomery_arithmetic\include\ ^
        %ecm_threshold% ^
        %trialdiv_size% ^
        %prbst_gcd_threshold%  %prbst_starting_length% ^
        %microecm_expect_large_factors%   %allow_microecm_dual_monty%  ^
        %cpp_standard%  %use_inline_asm_add%  %use_inline_asm_sub%  %use_alt_hr_addsubs%  ^
        %use_inline_asm_redc%  %use_all_inline_asm%  %trial_type% ^
        /Fe"test_single_factor" ^
        test_single_factor.cpp  microecm_c.c


echo "compilation finished, now executing"


test_single_factor.exe

