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
#if defined(__GNUC__)
#  include <cpuid.h>
#  include <string>
#  include <cstring>
#endif

#ifndef NDEBUG
#error "asserts are enabled and will slow performance"
#endif


#define STRINGIFY(x) #x
#define STRINGIFYMACRO(y) STRINGIFY(y)


#if defined(__GNUC__)
std::string displayCPU()
{
    // this code is copied from https://stackoverflow.com/a/50021699
    // licensed CC-BY-SA-3.0  https://creativecommons.org/licenses/by-sa/3.0/
    char CPUBrandString[0x40];
    unsigned int CPUInfo[4] = {0,0,0,0};

    __get_cpuid(0x80000000, &CPUInfo[0], &CPUInfo[1], &CPUInfo[2], &CPUInfo[3]);
    unsigned int nExIds = CPUInfo[0];

    std::memset(CPUBrandString, 0, sizeof(CPUBrandString));

    for (unsigned int i = 0x80000000; i <= nExIds; ++i)
    {
        __get_cpuid(i, &CPUInfo[0], &CPUInfo[1], &CPUInfo[2], &CPUInfo[3]);

        if (i == 0x80000002)
            memcpy(CPUBrandString, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000003)
            memcpy(CPUBrandString + 16, CPUInfo, sizeof(CPUInfo));
        else if (i == 0x80000004)
            memcpy(CPUBrandString + 32, CPUInfo, sizeof(CPUInfo));
    }
    std::string result("CPU Type: ");
    result += CPUBrandString;
    return result;
}
#endif


// this function benchmarks factoring of numbers just under
// 1 << (bit_width_T - 1).  For example, if T is int64_t, T has a bit
// width of 64 and this function would benchmark numbers just under
// 1 << 63.
template <typename T>
void bench_halfrange()
{
   using namespace std::chrono;
   using dsec = duration<double>;
   bool impossible_happened = false;
   auto t0 = steady_clock::now();

   using U = typename hurchalla::extensible_make_unsigned<T>::type;
   T start = static_cast<T>((static_cast<U>(1) << 63) - 1);
   for (T x = start; x > start - 400000; x = x-2) {
      int num_factors;
      auto arr = hurchalla::factorize(x, num_factors);
      // We need to prevent the compiler from completely removing
      // the factorize calls due to arr never being used.
      // So we'll check arr[0] (which is never 0) just so it's used.
      if (arr[0] == 0) {
         impossible_happened = true;
         break;
      }
//         std::cout << "the factors of " << x << " are:" << "\n";
//         for (int j = 0; j < num_factors; ++j)
//            std::cout << arr[j] << "\n";
   }
   if (impossible_happened)
      std::cout << "impossible\n";

   auto t1 = steady_clock::now();
   std::cout << dsec(t1-t0).count() << "\n";
}


int main()
{
#if defined(__GNUC__)
   std::cout << displayCPU() << "\n";
#endif
    
   std::cout << "HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME is ";
   std::cout << STRINGIFYMACRO(HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME) << "\n";

   std::cout << "---started---\n";

   std::cout << "(using type int64_t)\n";
   for (int i=0; i<5; ++i) {
      bench_halfrange<std::int64_t>();
   }
   std::cout << "(using type uint64_t)\n";
   for (int i=0; i<5; ++i) {
      bench_halfrange<std::uint64_t>();
   }
}
