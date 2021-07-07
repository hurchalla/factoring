#!/bin/bash

./build_tests.sh -t -cgcc -mdebug

cov-build --dir cov-int cmake --build ./build/debug_gcc0 --clean-first --config Debug

tar czvf factoring.tgz cov-int

# upload the tar archive at https://scan.coverity.com/projects/hurchalla-factoring/builds/new
