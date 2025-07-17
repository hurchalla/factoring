@echo off

REM Copyright (c) 2025 Jeffrey Hurchalla.
REM This Source Code Form is subject to the terms of the Mozilla Public
REM License, v. 2.0. If a copy of the MPL was not distributed with this
REM file, You can obtain one at https://mozilla.org/MPL/2.0/.



REM You need to clone the util, factoring, and modular_arithmetic repos
REM from https://github.com/hurchalla


REM SET repo_directory TO THE DIRECTORY WHERE YOU CLONED THE HURCHALLA GIT
REM REPOSITORIES.  (or otherwise ensure the compiler /I flags correctly specify
REM the needed hurchalla include directories)
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
        /Fe"benchmark_bpsw" ^
        benchmark_bpsw.cpp


echo "compilation finished, now executing"


benchmark_bpsw.exe

