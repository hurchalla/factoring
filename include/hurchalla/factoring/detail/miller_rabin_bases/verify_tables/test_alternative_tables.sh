#!/bin/bash

clang++ -stdlib=libc++  -Wall -Wextra -Wpedantic  -ferror-limit=3  \
        -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-old-style-cast  \
        -O2  -DNDEBUG  -march=haswell -mtune=haswell \
        -std="c++17" \
        -I/home/jeff/Desktop/factoring/include \
        -I/home/jeff/Desktop/modular_arithmetic/modular_arithmetic/include \
        -I/home/jeff/Desktop/modular_arithmetic/montgomery_arithmetic/include \
        -I/home/jeff/Desktop/util/include \
        test_alternative_tables.cpp  -o test_alternative_tables

./test_alternative_tables  strong_psps_to_2_64.txt
