@echo off


REM set expect_arbitrary_factors="-DEXPECT_ARBITRARY_SIZE_FACTORS"

REM set allow_microecm_dual_monty="-DHURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM"

REM set use_c_interface="-DUSE_ECM_C_INTERFACE"


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
        %expect_arbitrary_factors%   %allow_microecm_dual_monty%  ^
        %use_c_interface% ^
        /Fe"test_single_factor" ^
        test_single_factor.cpp  microecm_c.c


echo "compilation finished, now executing"


test_single_factor.exe

