// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "sample_primes_and_nonprimes.h"
#include "hurchalla/factoring/detail/is_prime_bruteforce.h"
#include "hurchalla/factoring/IsPrimeIntensive.h"

#include "gtest/gtest.h"
#include <cstdint>
#include <limits>
#include <cstddef>
#include <sstream>

namespace {


TEST(HurchallaFactoringIsPrimeIntensive, exhaustive_uint16_t) {
    namespace hc = hurchalla;
    using T = std::uint16_t;
    {
        hc::IsPrimeIntensive<T, true> isprime;
        for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(isprime(x) == hc::is_prime_bruteforce(x));
        }
        T x = std::numeric_limits<T>::max();
        EXPECT_TRUE(isprime(x) == hc::is_prime_bruteforce(x));
    }
    {
        hc::IsPrimeIntensive<T, false> isprime;
        for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(isprime(x) == hc::is_prime_bruteforce(x));
        }
        T x = std::numeric_limits<T>::max();
        EXPECT_TRUE(isprime(x) == hc::is_prime_bruteforce(x));
    }
}


template <typename T, bool OPTIMIZE_PRIMES>
void test_sample_primes_and_nonprimes()
{
    namespace hc = hurchalla;
    hc::IsPrimeIntensive<T, OPTIMIZE_PRIMES> isprime;

    constexpr int NUM_PRIMES =
                             sizeof(prime_numbers64)/sizeof(prime_numbers64[0]);
    for (int i=0; i<NUM_PRIMES; ++i) {
        if (prime_numbers64[i] <= hc::ut_numeric_limits<T>::max())
            EXPECT_TRUE(isprime(static_cast<T>(prime_numbers64[i])));
    }

    constexpr int NUM_NONPRIMES =
                       sizeof(nonprime_numbers64)/sizeof(nonprime_numbers64[0]);
    for (int i=0; i<NUM_NONPRIMES; ++i) {
        if (nonprime_numbers64[i] <= hc::ut_numeric_limits<T>::max())
            EXPECT_FALSE(isprime(static_cast<T>(nonprime_numbers64[i])));
    }
}

#if HURCHALLA_COMPILER_HAS_UINT128_T()
template <bool OPTIMIZE_PRIMES>
void test_sample_primes_and_nonprimes128()
{
    namespace hc = hurchalla;
    using T = __uint128_t;
    hc::IsPrimeIntensive<T, OPTIMIZE_PRIMES> isprime;

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
    // IsPrimeIntensive with a template argument value of false is exactly the
    // same as for template argument true.  We test it with true above (at
    // least in release config), so essentially we'd be repeating the exact
    // same test if we tested it here with false.
    test_sample_primes_and_nonprimes<std::uint32_t, false>();
#endif
    test_sample_primes_and_nonprimes<std::uint64_t, false>();

#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_sample_primes_and_nonprimes128<true>();
    test_sample_primes_and_nonprimes128<false>();
#endif
}


} // end namespace
