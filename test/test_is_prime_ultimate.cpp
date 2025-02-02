// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

// In debug mode, constructing the IsPrimeTable<std::uint32_t> objects takes too
// long to ordinarily be worth doing.
#if (defined(_MSC_VER) && !defined(_DEBUG)) || defined(__OPTIMIZE__)


#include "sample_primes_and_nonprimes.h"
#include "hurchalla/factoring/is_prime.h"
#include "hurchalla/factoring/resource_intensive_api/is_prime_ultimate.h"
#include "hurchalla/factoring/resource_intensive_api/IsPrimeTable.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"

#include "gtest/gtest.h"
#include <cstdint>
#include <limits>
#include <cstddef>
#include <sstream>

namespace {


using namespace hurchalla;
using namespace hurchalla::detail;


TEST(HurchallaFactoringIsPrimeUltimate, exhaustive_uint16_t) {
    IsPrimeTable<std::uint32_t> table;
    using T = std::uint16_t;
    for (T x = 0; x < ut_numeric_limits<T>::max(); ++x) {
        SCOPED_TRACE(testing::Message() << "x == " << x);
        EXPECT_TRUE(is_prime_ultimate(x, table) == is_prime(x));
    }
    T x = ut_numeric_limits<T>::max();
    EXPECT_TRUE(is_prime_ultimate(x, table) == is_prime(x));
}


template <typename T>
void test_sample_primes_and_nonprimes(IsPrimeTable<std::uint32_t>& table)
{
    using U = typename extensible_make_unsigned<T>::type;
    U tmax = static_cast<U>(ut_numeric_limits<T>::max());
    constexpr int NUM_PRIMES =
                             sizeof(prime_numbers64)/sizeof(prime_numbers64[0]);
    for (int i=0; i<NUM_PRIMES; ++i) {
        if (prime_numbers64[i] <= tmax) {
            EXPECT_TRUE(is_prime_ultimate(static_cast<T>(prime_numbers64[i]), table));
        }
    }
    constexpr int NUM_NONPRIMES =
                       sizeof(nonprime_numbers64)/sizeof(nonprime_numbers64[0]);
    for (int i=0; i<NUM_NONPRIMES; ++i) {
      if (nonprime_numbers64[i] <= tmax) {
        EXPECT_FALSE(is_prime_ultimate(static_cast<T>(nonprime_numbers64[i]), table));
      }
    }
}


#if HURCHALLA_COMPILER_HAS_UINT128_T()
template <typename T>
void test_sample_primes_and_nonprimes128(IsPrimeTable<std::uint32_t>& table)
{
    using U = typename extensible_make_unsigned<T>::type;
    U tmax = static_cast<U>(ut_numeric_limits<T>::max());
    constexpr int NUM_PRIMES128 =
                           sizeof(prime_numbers128)/sizeof(prime_numbers128[0]);
    for (int i=0; i<NUM_PRIMES128; ++i) {
        if (prime_numbers128[i] <= tmax) {
            EXPECT_TRUE(is_prime_ultimate(static_cast<T>(prime_numbers128[i]), table));
        }
    }
    constexpr int NUM_NONPRIMES128 =
                     sizeof(nonprime_numbers128)/sizeof(nonprime_numbers128[0]);
    for (int i=0; i<NUM_NONPRIMES128; ++i) {
        if (nonprime_numbers128[i] <= tmax) {
            EXPECT_FALSE(is_prime_ultimate(static_cast<T>(nonprime_numbers128[i]), table));
        }
    }
}
#endif


TEST(HurchallaFactoringIsPrimeUltimate, basic_tests) {
    IsPrimeTable<std::uint32_t> table;

    test_sample_primes_and_nonprimes<std::uint8_t>(table);
    test_sample_primes_and_nonprimes<std::uint16_t>(table);
    test_sample_primes_and_nonprimes<std::uint32_t>(table);
    test_sample_primes_and_nonprimes<std::uint64_t>(table);

    test_sample_primes_and_nonprimes<std::int8_t>(table);
    test_sample_primes_and_nonprimes<std::int16_t>(table);
    test_sample_primes_and_nonprimes<std::int32_t>(table);
    test_sample_primes_and_nonprimes<std::int64_t>(table);

# if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_sample_primes_and_nonprimes<__uint128_t>(table);
    test_sample_primes_and_nonprimes128<__uint128_t>(table);

    test_sample_primes_and_nonprimes<__int128_t>(table);
    test_sample_primes_and_nonprimes128<__int128_t>(table);
# endif
}


} // end namespace

#endif
