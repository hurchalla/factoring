// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


#include "sample_primes_and_nonprimes.h"
#include "hurchalla/factoring/detail/is_prime_bruteforce.h"
#include "hurchalla/factoring/resource_intensive_api/IsPrimeTable.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"

#include "gtest/gtest.h"
#include <cstdint>
//#include <limits>
//#include <cstddef>
//#include <sstream>

namespace {


using namespace hurchalla;


TEST(HurchallaFactoringIsPrimeTable, exhaustive_uint16_t) {
    using T = std::uint16_t;
    {
        IsPrimeTable<T> isprime;
        for (T x = 0; x < ut_numeric_limits<T>::max(); ++x) {
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(isprime(x) == is_prime_bruteforce::call(x));
        }
        T x = ut_numeric_limits<T>::max();
        EXPECT_TRUE(isprime(x) == is_prime_bruteforce::call(x));
    }
}


template <typename T>
void test_sample_primes_and_nonprimes()
{
    IsPrimeTable<T> isprime;

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


TEST(HurchallaFactoringIsPrimeTable, basic_tests) {
    test_sample_primes_and_nonprimes<std::uint8_t>();
    test_sample_primes_and_nonprimes<std::uint16_t>();
// this is a bit of a hack, but this test is extremely slow without optimization
// and so we try to skip it when building without any optimizations
#if (defined(_MSC_VER) && !defined(_DEBUG)) || defined(__OPTIMIZE__)
    test_sample_primes_and_nonprimes<std::uint32_t>();
#endif
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    // silence unused variable warnings
    static_cast<void>(prime_numbers128[0]);
    static_cast<void>(nonprime_numbers128[0]);
#endif
}


} // end namespace
