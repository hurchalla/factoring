// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// This example is intended for the case that you are not using CMake.
// If you haven't already, you need to follow the steps in the README.md
// for "How to use the library" | "Without CMake"
#include "hurchalla/factoring/factorize.h"
#include <iostream>

#ifndef NDEBUG
// remove this if you want to allow asserts
// (they're great for testing and debugging).
#error "asserts are enabled and will slow performance"
#endif

int main()
{
   unsigned int x = 322u;

   int num_factors;
   auto array = hurchalla::factorize(x, num_factors);
   std::cout << "the factors of " << x << " are:" << "\n";
   for (int i = 0; i < num_factors; ++i)
      std::cout << array[i] << "\n";
   return 0;
}
