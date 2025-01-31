// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
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



// this function benchmarks factoring of odd numbers between min and max.
template <typename T>
void bench_range(T min, T max)
{
   if (max % 2 == 0)
      --max;
   if (min == 0)
      min = 1;
   using namespace std::chrono;
   using dsec = duration<double>;
   bool impossible_happened = false;
   auto t0 = steady_clock::now();

   for (T x = max; x > min; x = x-2) {
      unsigned int num_factors;
      auto arr = hurchalla::factorize(x, num_factors);
      // We need to prevent the compiler from completely removing
      // the factorize calls due to the array never being used.
      // So we'll check arr[0] (which is never 0) just so it's used.
      if (arr[0] == 0) {
         impossible_happened = true;
         break;
      }
#if 0
      std::cout << "the factors of " << x << " are:" << "\n";
      for (unsigned int j = 0; j < num_factors; ++j)
         std::cout << arr[j] << "\n";
#endif
   }
   if (impossible_happened)
      std::cout << "impossible\n";

   auto t1 = steady_clock::now();
   std::cout << dsec(t1-t0).count() << "\n";
}


template <typename T>
void print_int_type()
{
   if constexpr (hurchalla::ut_numeric_limits<T>::is_signed) {
        std::cout << "signed int" <<
                hurchalla::ut_numeric_limits<T>::digits + 1;
   } else {
        std::cout << "unsigned int" <<
                hurchalla::ut_numeric_limits<T>::digits;
   }
}


int main()
{
#if defined(__GNUC__)
   std::cout << displayCPU() << "\n";
#endif
    
   std::cout << "HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME is ";
   std::cout << STRINGIFYMACRO(HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME) << "\n";

   std::cout << "---started---\n";


   int num_test_runs = 5;
   using T = int64_t;
   T span = 400000;


   T max = hurchalla::ut_numeric_limits<T>::max();
   if (max < span) {
       std::cout << "Error: max < span\n";
       return 1;
   }
   T min = max - span;

   std::cout << "using ";  print_int_type<T>();  std::cout << "\n";
   for (int i=0; i<num_test_runs; ++i) {
      bench_range(min, max);
   }
}

