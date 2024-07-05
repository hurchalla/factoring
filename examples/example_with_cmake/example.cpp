// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// This example is intended for the case that you are using CMake.
// If you haven't already, you need to follow the steps in the README.md
// for "How to use the library" | "With CMake"
#include "hurchalla/factoring/factorize.h"
#include <iostream>

#ifndef NDEBUG
// remove this if you want to allow asserts
// (they're very good for testing and debugging but may drastically slow down
// the library).
#error "Performance warning: asserts are enabled and will slow performance"
#endif

int main()
{
   unsigned int x = 322u;

   unsigned int num_factors;
   auto array = hurchalla::factorize(x, num_factors);
   std::cout << "the factors of " << x << " are:" << "\n";
   for (unsigned int i = 0; i < num_factors; ++i)
      std::cout << array[i] << "\n";
   return 0;
}
