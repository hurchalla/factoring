#!/bin/bash

# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


./build_tests.sh -t -cgcc -mdebug

cov-build --dir cov-int cmake --build ./build/debug_gcc0 --clean-first --config Debug

tar czvf factoring.tgz cov-int

# upload the tar archive at https://scan.coverity.com/projects/hurchalla-factoring/builds/new
