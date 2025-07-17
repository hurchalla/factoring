// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/factoring/detail/is_prime_trialdivision.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/util/traits/extensible_make_signed.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
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
#define GCD_FUNCTION_NAME hurchalla::greatest_common_divisor
#define LUCAS_FORCE_INLINE HURCHALLA_FORCE_INLINE
#define NUM_LIMITS hurchalla::ut_numeric_limits
#define TRAIT_MAKE_SIGNED hurchalla::extensible_make_signed
#define TRAIT_MAKE_UNSIGNED hurchalla::extensible_make_unsigned
//
// I'm not sure how the Perl5 license from dana_lucas.h interacts with the
// MPL2 license of this file.  (This file is MPL2 because it uses header
// libraries (factoring, util, montgomery_arithmetic) that are MPL2.)
//
// So by default the following #include
// is disabled.
// For purposes of your own experiments, you can enable it...
#if 0
# include "dana_lucas.h"
#else
# error "in order to compile, you must enable this section"
#endif



// --- The macros below here change what is being benchmarked ---

// You can define one of these two macros to isolate the performance of the Lucas portion of
// BPSW, or the Miller-Rabin portion of BPSW.  Note: the perf results also include the Trial
// division and Montgomery initialization cost.
// If you define both macros, you get perf results for only the trial division and mont init.
//
//#define SKIP_MILLER_RABIN
//#define SKIP_LUCAS


// You can define this macro to measure perf for hashed Miller Rabin, which is an alternative
// to BPSW (for uint64_t and smaller integer types).
// Hashed MR uses a 448kb hash table, but can potentially be faster than BPSW.
//
//#define MEASURE_HASHED_MR64_PLUS_TD_AND_MONT_INIT


// These settings for TRIAL_DIVISION_SIZE roughly produced the best measured
// performance for me on Haswell(x64 CPU), MSVC 64bit compiler Windows, testing
// a range of 64 bit numbers that was essentially arbitrary (numbers with no
// particular likelihood to be composite or prime).
// Different systems will have different best settings of TRIAL_DIVISION_SIZE...
#ifndef MEASURE_HASHED_MR64_PLUS_TD_AND_MONT_INIT
# define TRIAL_DIVISION_SIZE (110)
//34bit
//115 5032
//50  5041
//80  5012
//150 5089
//95  5005
//
//64bit
// 50  5897
// 75  5730
// 100 5636
// 125 5615
// 150 5647
// 200 5682
#else
# define TRIAL_DIVISION_SIZE (150)
//34bit
//160 4403
//100 4329 4331 4317
//75  4317
//50  4399
//130 4313 4301
//115 4310
//150 4363
//
//64bit
// 100  5380 5349
// 125  5327
// 150  5327 5306
// 160  5285
// 175  5314
// 200  5315
// 250  5351
// 300  5397


//34bit (*120*)

//32 bit   (*88*)
//150 3136
//125 3095 (3083)
//100 3068
//75  3070 (3066)
//50  3100

//28bit  (*75*)
//125 2957
//100 2922
// 90 2915 (2906) (2872)
// 75 2859 (2890)
// 60 2914 (2891)
// 50 2928
// 35 2928 2916
#endif










// For BPSW's lucas test, LUCAS_STRENGTH can be set from 0 to 4.
// The *much* preferred (fastest) lucas test is the almost-extra-strong version,
// which is LUCAS_STRENGTH of 4.  Though measured performance differs from system to
// system, it would surprise me if you measure any LUCAS_STRENGTH != 4 to be fastest.
#define LUCAS_STRENGTH (4)
//
#if !defined(LUCAS_STRENGTH) || (LUCAS_STRENGTH < 0 || LUCAS_STRENGTH > 4)
# error "LUCAS_STRENGTH must be defined as 0 1 2 3 or 4."
#endif


// You can define this macro to (approximately) measure how long it takes BPSW to test
// prime numbers, regardless of whether the numbers you give it are prime.  Normally BPSW
// would not run the lucas test when it determines a number is composite via Miller-Rabin.
// If you define this macro, it will still run the Lucas test even in those cases.
//
//#define ALWAYS_USE_LUCAS_IN_BPSW
//
#if defined(SKIP_LUCAS) && defined(ALWAYS_USE_LUCAS_IN_BPSW)
# error "You can not define both SKIP_LUCAS and ALWAYS_USE_LUCAS_IN_BPSW"
#endif







template <typename MontType>
static bool bpsw_montgomery(const MontType& mf)
{
  using T = typename MontType::IntegerType;
  static_assert(hurchalla::ut_numeric_limits<T>::is_integer, "");
  HPBC_PRECONDITION2(mf.getModulus() > 1);

  T modulus = mf.getModulus();

#if defined(TRIAL_DIVISION_SIZE) && (TRIAL_DIVISION_SIZE > 0)
  {
    namespace hcd = hurchalla::detail;
    bool success;
    bool isPrime = hcd::is_prime_trialdivision::
             call<hcd::PrimeTrialDivisionWarren, TRIAL_DIVISION_SIZE>(modulus, success);
    if (success)
      return isPrime;
    // is_prime_trialdivision::call should have successfully handled any even
    // modulus, and any modulus < 2.
    HPBC_ASSERT2(modulus % 2 != 0);
    HPBC_ASSERT2(modulus >= 2);
  }
#else
  {
    if ((modulus % 2) == 0) return (modulus == 2);
  }
#endif

  {
    // Ensure preconditions are met for  is_lucas_probable_prime()
    HPBC_ASSERT2(modulus % 2 != 0);
    T n = modulus;
    if (n < 13) return (n == 2 || n == 3 || n == 5 || n == 7 || n == 11);
    using UT = typename hurchalla::extensible_make_unsigned<T>::type;
    auto UT_MAX = hurchalla::ut_numeric_limits<UT>::max();
    // since (1<<4) == 16, squaring that value gives a result that ends in 6.
    // Squaring that result again gives a result ending in 6...
    // And if we subtract 1 from any of those, the ending digit is 5, and
    // therefore the difference is composite.
    // So (1<<8) - 1 is composite.  As is (1<<16) - 1, and (1<<32) - 1, etc.
    // These are the possible values of UT_MAX, and are all composite.
    if (n == UT_MAX) return false;
  }


#ifdef MEASURE_HASHED_MR64_PLUS_TD_AND_MONT_INIT
  {
    constexpr std::size_t TOTAL_BASES = 2;
    constexpr std::size_t TRIAL_SIZE = 2;
    const auto bases = hurchalla::detail::MillerRabinBases<64, TOTAL_BASES>::get(modulus);
    return hurchalla::detail::IPMR_internal::miller_rabin_trials<TRIAL_SIZE>(mf, bases);
  }
#endif


#if defined(SKIP_MILLER_RABIN) && defined(SKIP_LUCAS)
// Measure only the time to perform the trial division and montgomery initialization.
#  ifdef MEASURE_HASHED_MR64_PLUS_TD_AND_MONT_INIT
#    error "do not define MEASURE_HASHED_MR64_PLUS_TD_AND_MONT_INIT if you define both SKIP_LUCAS and SKIP_MILLER_RABIN"
#  endif
  // we want to force montgomery initialization to happen, but
  // almost nothing more than that.
  auto val = mf.convertIn(1);
  // this should prevent the compiler optimizer from omitting the calculation of val.
  // we expect this to always be true, but the compiler doesn't know that.
  return (mf.getCanonicalValue(val) != mf.getZeroValue());
#endif


  using C = typename MontType::CanonicalValue;
  const C one = mf.getUnityValue();
  const C two = mf.add(one, one);
  bool is_probable_prime = true;


#ifndef SKIP_MILLER_RABIN
  constexpr std::size_t TOTAL_BASES = 1;
  constexpr std::size_t TRIAL_SIZE = 1;
  // Base 2 MR test
  const std::array<std::uint8_t, 1> bases = { UINT8_C(2) };
  is_probable_prime = hurchalla::detail::IPMR_internal::miller_rabin_trials<TRIAL_SIZE>(mf, bases);
#  ifndef ALWAYS_USE_LUCAS_IN_BPSW
  if (!is_probable_prime)
    return false;
#  endif
#endif

#ifdef SKIP_LUCAS
  return is_probable_prime;
#endif


  bool islpp = is_lucas_probable_prime<LUCAS_STRENGTH>(mf, two); // !!!! You'll get a compile error unless you enable #include "dana_lucas.h" at the top of this file

  return (islpp && is_probable_prime);
}




template <typename T>
static bool is_prime_bpsw(T x)
{
  static_assert(hurchalla::ut_numeric_limits<T>::is_integer, "");
  constexpr int TDIGITS = hurchalla::ut_numeric_limits<T>::digits;
  static_assert(TDIGITS == 63 || TDIGITS == 64);
  constexpr T Rdiv4 = static_cast<T>(1) << 62;

  HPBC_PRECONDITION2(x % 2 == 1);
  HPBC_PRECONDITION2(x > 1);

  if (x < Rdiv4)
      return bpsw_montgomery(hurchalla::MontgomeryQuarter<T,true>(static_cast<T>(x)));
  else
      return bpsw_montgomery(hurchalla::MontgomeryForm<T,true>(static_cast<T>(x)));
}







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



// this function benchmarks primality testing of odd numbers between min and max.
template <typename T>
void bench_range(T min, T max)
{
   if (max % 2 == 0)
      --max;
   if (min == 0)
      min = 1;
   using namespace std::chrono;
   using dsec = duration<double>;
   auto t0 = steady_clock::now();

   T total_primes = 0;

   for (T x = max; x > min; x = x-2) {
      bool isprime = is_prime_bpsw(x);
      if (isprime)
         total_primes++;
#if 0
      std::cout << "isprime(" << x << ") == " << isprime << "\n";
#endif
   }

   auto t1 = steady_clock::now();
   std::cout << "time: " << dsec(t1-t0).count() << "\n";
   std::cout << "total_primes == " << total_primes << "\n";
}




int main()
{
#if defined(__GNUC__)
   std::cout << displayCPU() << "\n";
#endif
   std::cout << "---started---\n";


   int num_test_runs = 5;
   using T = uint64_t;
   T span = 10000000;

   std::cout << "HURCHALLA_TARGET_BIT_WIDTH == " << HURCHALLA_TARGET_BIT_WIDTH << "\n";


   T max = hurchalla::ut_numeric_limits<T>::max() - 2;
//   T max = static_cast<T>(static_cast<T>(1) << 34);
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
