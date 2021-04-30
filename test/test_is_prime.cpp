// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "sample_primes_and_nonprimes.h"
#include "hurchalla/factoring/detail/is_prime_bruteforce.h"
#include "hurchalla/factoring/is_prime.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"

#include "gtest/gtest.h"
#include <cstdint>
#include <limits>
#include <cstddef>
#include <sstream>

namespace {


TEST(HurchallaFactoringIsPrime, exhaustive_uint16_t) {
    namespace hc = hurchalla;
    using T = std::uint16_t;
    for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
        SCOPED_TRACE(testing::Message() << "x == " << x);
        EXPECT_TRUE(hc::is_prime(x) == hc::is_prime_bruteforce(x));
    }
    T x = std::numeric_limits<T>::max();
    EXPECT_TRUE(hc::is_prime(x) == hc::is_prime_bruteforce(x));
}

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 10/20/20.
TEST(HurchallaFactoringIsPrime, exhaustive_uint32_t) {
    namespace hc = hurchalla;
    using T = std::uint32_t;
    for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
        EXPECT_TRUE(hc::is_prime(x) == hc::is_prime_bruteforce(x));
    }
    T x = std::numeric_limits<T>::max();
    EXPECT_TRUE(hc::is_prime(x) == hc::is_prime_bruteforce(x));
}
#endif

#if 0
// I used this speed test to do a quick and dirty initial performance tuning of
// the default value for the macro HURCHALLA_ISPRIME_TRIALDIV_SIZE, used in
// impl_is_prime.h.  This test is not needed normally.
TEST(HurchallaFactoringIsPrime, speed_test32) {
    namespace hc = hurchalla;
    using T = std::uint32_t;
#if 1
    T max = std::numeric_limits<T>::max();
#else
    T max = 20065537;
#endif
    T dummy = 0;
    for (T x = max; x >= max - 20000000; x = x-2) {
        bool b = hc::is_prime(x);
        // We need to prevent the compiler from completely removing
        // the is_prime calls due to b never being used.
        // So we'll add b to dummy just so it's used.
        dummy += b;
    }
    EXPECT_TRUE(dummy > 0);
}
TEST(HurchallaFactoringIsPrime, speed_test64) {
    namespace hc = hurchalla;
    using T = std::uint64_t;
#if 1
    T max = std::numeric_limits<T>::max();
#else
    T max = 10065537;
#endif
    T dummy = 0;
    for (T x = max; x >= max - 10000000; x = x-2) {
        bool b = hc::is_prime(x);
        // We need to prevent the compiler from completely removing
        // the is_prime calls due to b never being used.
        // So we'll add b to dummy just so it's used.
        dummy += b;
    }
    EXPECT_TRUE(dummy > 0);
}
#endif


template <typename T>
void test_sample_primes_and_nonprimes()
{
    namespace hc = hurchalla;
    constexpr int NUM_PRIMES =
                             sizeof(prime_numbers64)/sizeof(prime_numbers64[0]);
    for (int i=0; i<NUM_PRIMES; ++i) {
        if (prime_numbers64[i] <= hc::ut_numeric_limits<T>::max()) {
            EXPECT_TRUE(hc::is_prime(static_cast<T>(prime_numbers64[i])));
        }
    }

    constexpr int NUM_NONPRIMES =
                       sizeof(nonprime_numbers64)/sizeof(nonprime_numbers64[0]);
    for (int i=0; i<NUM_NONPRIMES; ++i) {
        if (nonprime_numbers64[i] <= hc::ut_numeric_limits<T>::max()) {
            EXPECT_FALSE(hc::is_prime(static_cast<T>(nonprime_numbers64[i])));
        }
    }
}

TEST(HurchallaFactoringIsPrime, basic_tests) {
    namespace hc = hurchalla;

    test_sample_primes_and_nonprimes<std::uint8_t>();
    test_sample_primes_and_nonprimes<std::uint16_t>();
    test_sample_primes_and_nonprimes<std::uint32_t>();
    test_sample_primes_and_nonprimes<std::uint64_t>();

#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_sample_primes_and_nonprimes<__uint128_t>();
    using T = __uint128_t;
    constexpr int NUM_PRIMES128 =
                           sizeof(prime_numbers128)/sizeof(prime_numbers128[0]);
    for (int i=0; i<NUM_PRIMES128; ++i) {
        EXPECT_TRUE(hc::is_prime(static_cast<T>(prime_numbers128[i])));
    }
    constexpr int NUM_NONPRIMES128 =
                     sizeof(nonprime_numbers128)/sizeof(nonprime_numbers128[0]);
    for (int i=0; i<NUM_NONPRIMES128; ++i) {
        EXPECT_FALSE(hc::is_prime(static_cast<T>(nonprime_numbers128[i])));
    }
#endif
}


} // end namespace
