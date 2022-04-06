// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/factoring/factorize.h"
#include <iostream>
#include <chrono>

#ifndef NDEBUG
#error "asserts are enabled and will slow performance"
#endif


int main()
{
   using T = std::int64_t;

   using namespace std::chrono;
   using dsec = duration<double>;

   bool impossible_happened = false;

   std::cout << "---started---\n";
   auto t0 = steady_clock::now();

   using U = typename hurchalla::extensible_make_unsigned<T>::type;
   T start = static_cast<T>((static_cast<U>(1) << 63) - 1);
   for (T x = start; x > start - 800000; x = x-2) {
      int num_factors;
      auto arr = hurchalla::factorize(x, num_factors);
      // We need to prevent the compiler from completely removing
      // the factorize calls due to arr never being used.
      // So we'll check arr[0] (which is never 0) just so it's used.
      if (arr[0] == 0) {
         impossible_happened = true;
         break;
      }

//      std::cout << "the factors of " << x << " are:" << "\n";
//      for (int i = 0; i < num_factors; ++i)
//         std::cout << arr[i] << "\n";
   }
   if (impossible_happened)
      std::cout << "impossible\n";

   auto t1 = steady_clock::now();
   std::cout << dsec(t1-t0).count() << "\n";
}
