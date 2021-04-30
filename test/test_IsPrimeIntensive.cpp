// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "sample_primes_and_nonprimes.h"
#include "hurchalla/factoring/detail/is_prime_bruteforce.h"
#include "hurchalla/factoring/IsPrimeIntensive.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"

#include "gtest/gtest.h"
#include <cstdint>
#include <limits>
#include <cstddef>
#include <sstream>

namespace {


using namespace hurchalla;

#if 0
// I used this speed test to do a quick and dirty initial performance tuning of
// the default value for the macro HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE,
// used in ImplIsPrimeIntensive.h.  This test is not needed normally.
template <typename T>
void speed_test(T testlen)
{
    IsPrimeIntensive<T, false> isprime;
// choose between testing extremely large numbers, or extremely small
#  if 1
    T max = ut_numeric_limits<T>::max();    // extremely large
#  else
    if (testlen % 2 == 0)
       ++testlen;
    ASSERT_TRUE(testlen < ut_numeric_limits<T>::max() - 10000000);
    T max = testlen + 10000000;                 // extremely small
#  endif
    T dummy = 0;
    for (T x = max; x >= max - testlen; x = x-2) {
        bool b = isprime(x);
        // We need to prevent the compiler from completely removing
        // the is_prime calls due to b never being used.
        // So we'll add b to dummy just so it's used.
        dummy += b;
    }
    EXPECT_TRUE(dummy > 0);
}
TEST(HurchallaFactoringIsPrimeIntensive, speed_test64) {
    speed_test(static_cast<std::uint64_t>(10000000));
}
#  if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(HurchallaFactoringIsPrimeIntensive, speed_test128) {
    speed_test(static_cast<__uint128_t>(200000));
}
#  endif
#endif


TEST(HurchallaFactoringIsPrimeIntensive, exhaustive_uint16_t) {
    using T = std::uint16_t;
    {
        IsPrimeIntensive<T, true> isprime;
        for (T x = 0; x < ut_numeric_limits<T>::max(); ++x) {
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(isprime(x) == is_prime_bruteforce(x));
        }
        T x = ut_numeric_limits<T>::max();
        EXPECT_TRUE(isprime(x) == is_prime_bruteforce(x));
    }
    {
        IsPrimeIntensive<T, false> isprime;
        for (T x = 0; x < ut_numeric_limits<T>::max(); ++x) {
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(isprime(x) == is_prime_bruteforce(x));
        }
        T x = ut_numeric_limits<T>::max();
        EXPECT_TRUE(isprime(x) == is_prime_bruteforce(x));
    }
}


template <typename T, bool OPTIMIZE_PRIMES>
void test_sample_primes_and_nonprimes()
{
    IsPrimeIntensive<T, OPTIMIZE_PRIMES> isprime;

    constexpr int NUM_PRIMES =
                             sizeof(prime_numbers64)/sizeof(prime_numbers64[0]);
    for (int i=0; i<NUM_PRIMES; ++i) {
        if (prime_numbers64[i] <= ut_numeric_limits<T>::max()) {
            EXPECT_TRUE(isprime(static_cast<T>(prime_numbers64[i])));
        }
    }
    constexpr int NUM_NONPRIMES =
                       sizeof(nonprime_numbers64)/sizeof(nonprime_numbers64[0]);
    for (int i=0; i<NUM_NONPRIMES; ++i) {
        if (nonprime_numbers64[i] <= ut_numeric_limits<T>::max()) {
            EXPECT_FALSE(isprime(static_cast<T>(nonprime_numbers64[i])));
        }
    }
}

#if HURCHALLA_COMPILER_HAS_UINT128_T()
template <bool OPTIMIZE_PRIMES>
void test_sample_primes_and_nonprimes128()
{
    using T = __uint128_t;
    IsPrimeIntensive<T, OPTIMIZE_PRIMES> isprime;

    constexpr int NUM_PRIMES128 =
                           sizeof(prime_numbers128)/sizeof(prime_numbers128[0]);
    for (int i=0; i<NUM_PRIMES128; ++i) {
        EXPECT_TRUE(isprime(static_cast<T>(prime_numbers128[i])));
    }
    constexpr int NUM_NONPRIMES128 =
                     sizeof(nonprime_numbers128)/sizeof(nonprime_numbers128[0]);
    for (int i=0; i<NUM_NONPRIMES128; ++i) {
        EXPECT_FALSE(isprime(static_cast<T>(nonprime_numbers128[i])));
    }
}
#endif


TEST(HurchallaFactoringIsPrimeIntensive, basic_tests) {
    test_sample_primes_and_nonprimes<std::uint8_t, true>();
    test_sample_primes_and_nonprimes<std::uint16_t, true>();
#ifdef NDEBUG
    test_sample_primes_and_nonprimes<std::uint32_t, true>();
#endif
    test_sample_primes_and_nonprimes<std::uint64_t, true>();

    test_sample_primes_and_nonprimes<std::uint8_t, false>();
    test_sample_primes_and_nonprimes<std::uint16_t, false>();
#if 0
    // I verified this ran without error in both release and debug configs.
    // It's disabled by default here because it takes ~5 seconds to run in
    // release and 1.5 minutes in debug, and this is not worth the run-time
    // considering that the current and forseeable implementation of
    // IsPrimeIntensive<std::uint32_t, bool_value> with "bool_value" set to
    // false is exactly the same as for "bool_value" set to true.  We test
    // it with true above (when NDEBUG is defined), so essentially we'd be
    // repeating the exact same test if we tested it again here with false.
    test_sample_primes_and_nonprimes<std::uint32_t, false>();
#endif
    test_sample_primes_and_nonprimes<std::uint64_t, false>();

#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_sample_primes_and_nonprimes128<true>();
    test_sample_primes_and_nonprimes128<false>();
#endif
}


} // end namespace
