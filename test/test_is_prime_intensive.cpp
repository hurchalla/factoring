// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "sample_primes_and_nonprimes.h"
#include "hurchalla/factoring/is_prime.h"
#include "hurchalla/factoring/resource_intensive_api/is_prime_intensive.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"

#include "gtest/gtest.h"
#include <cstdint>
#include <limits>
#include <cstddef>
#include <sstream>

namespace {


using namespace hurchalla;
using namespace hurchalla::detail;


TEST(HurchallaFactoringIsPrimeIntense, exhaustive_uint8_t) {
    using T = std::uint8_t;
    for (T x = 0; x < ut_numeric_limits<T>::max(); ++x) {
        SCOPED_TRACE(testing::Message() << "x == " << x);
        EXPECT_TRUE(is_prime_intensive(x) == is_prime(x));
    }
    T x = ut_numeric_limits<T>::max();
    EXPECT_TRUE(is_prime_intensive(x) == is_prime(x));
}

TEST(HurchallaFactoringIsPrimeIntense, exhaustive_uint16_t) {
    using T = std::uint16_t;
    for (T x = 0; x < ut_numeric_limits<T>::max(); ++x) {
        SCOPED_TRACE(testing::Message() << "x == " << x);
        EXPECT_TRUE(is_prime_intensive(x) == is_prime(x));
    }
    T x = ut_numeric_limits<T>::max();
    EXPECT_TRUE(is_prime_intensive(x) == is_prime(x));
}

#if 0
    // Ordinarily you don't want to run this since it takes a while to
    // complete.  It passed when I tested it on 2/1/25.
    TEST(HurchallaFactoringIsPrimeIntense, exhaustive_uint32_t) {
        using T = std::uint32_t;
        for (T x = 0; x < ut_numeric_limits<T>::max(); ++x) {
            EXPECT_TRUE(is_prime_intensive(x) == is_prime(x));
        }
        T x = ut_numeric_limits<T>::max();
        EXPECT_TRUE(is_prime_intensive(x) == is_prime(x));
    }
#endif



template <typename T>
void test_sample_primes_and_nonprimes()
{
    using U = typename extensible_make_unsigned<T>::type;
    U tmax = static_cast<U>(ut_numeric_limits<T>::max());

    constexpr int NUM_PRIMES =
                             sizeof(prime_numbers64)/sizeof(prime_numbers64[0]);
    for (int i=0; i<NUM_PRIMES; ++i) {
        if (prime_numbers64[i] <= tmax) {
            EXPECT_TRUE(is_prime_intensive(static_cast<T>(prime_numbers64[i])));
        }
    }
    constexpr int NUM_NONPRIMES =
                       sizeof(nonprime_numbers64)/sizeof(nonprime_numbers64[0]);
    for (int i=0; i<NUM_NONPRIMES; ++i) {
      if (nonprime_numbers64[i] <= tmax) {
        EXPECT_FALSE(is_prime_intensive(static_cast<T>(nonprime_numbers64[i])));
      }
    }
}


TEST(HurchallaFactoringIsPrimeIntense, basic_tests) {
    test_sample_primes_and_nonprimes<std::uint8_t>();
    test_sample_primes_and_nonprimes<std::uint16_t>();
    test_sample_primes_and_nonprimes<std::uint32_t>();
    test_sample_primes_and_nonprimes<std::uint64_t>();

    test_sample_primes_and_nonprimes<std::int8_t>();
    test_sample_primes_and_nonprimes<std::int16_t>();

    test_sample_primes_and_nonprimes<std::int32_t>();
    test_sample_primes_and_nonprimes<std::int64_t>();

#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_sample_primes_and_nonprimes<__uint128_t>();
    using T = __uint128_t;
    constexpr int NUM_PRIMES128 =
                           sizeof(prime_numbers128)/sizeof(prime_numbers128[0]);
    for (int i=0; i<NUM_PRIMES128; ++i) {
        EXPECT_TRUE(is_prime_intensive(static_cast<T>(prime_numbers128[i])));
    }
    constexpr int NUM_NONPRIMES128 =
                     sizeof(nonprime_numbers128)/sizeof(nonprime_numbers128[0]);
    for (int i=0; i<NUM_NONPRIMES128; ++i) {
       EXPECT_FALSE(is_prime_intensive(static_cast<T>(nonprime_numbers128[i])));
    }
#endif

}


} // end namespace
