// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "is_prime_wheel210.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <cstddef>

namespace {


using namespace hurchalla;
using namespace hurchalla::detail;


TEST(HurchallaFactoringIsPrimeMillerRabin, super_simple_test) {
    std::uint64_t modulus = 53;
    MontgomeryForm<std::uint64_t> mf(modulus);
    EXPECT_TRUE(is_prime_miller_rabin::call_mont(mf));
}

TEST(HurchallaFactoringIsPrimeMillerRabin, integer_tests) {
    using is_prime_mr = is_prime_miller_rabin;
    std::int16_t primenum = 59;
    std::int16_t composite = 63;
    EXPECT_TRUE(is_prime_mr::call(primenum));
    EXPECT_FALSE(is_prime_mr::call(composite));
    EXPECT_TRUE(is_prime_mr::call(static_cast<std::uint16_t>(primenum)));
    EXPECT_FALSE(is_prime_mr::call(static_cast<std::uint16_t>(composite)));
    EXPECT_TRUE(is_prime_mr::call(static_cast<std::int32_t>(primenum)));
    EXPECT_FALSE(is_prime_mr::call(static_cast<std::int32_t>(composite)));
    EXPECT_TRUE(is_prime_mr::call(static_cast<std::uint32_t>(primenum)));
    EXPECT_FALSE(is_prime_mr::call(static_cast<std::uint32_t>(composite)));
    EXPECT_TRUE(is_prime_mr::call(static_cast<std::int64_t>(primenum)));
    EXPECT_FALSE(is_prime_mr::call(static_cast<std::int64_t>(composite)));
    EXPECT_TRUE(is_prime_mr::call(static_cast<std::uint64_t>(primenum)));
    EXPECT_FALSE(is_prime_mr::call(static_cast<std::uint64_t>(composite)));
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    __uint128_t prime_u127 = (static_cast<__uint128_t>(1) << 127) - 1;
    __uint128_t composite_u127 = prime_u127 - 2;
    EXPECT_TRUE(is_prime_mr::call(static_cast<__int128_t>(prime_u127)));
    EXPECT_FALSE(is_prime_mr::call(static_cast<__int128_t>(composite_u127)));
    EXPECT_TRUE(is_prime_mr::call(prime_u127));
    EXPECT_FALSE(is_prime_mr::call(composite_u127));
#endif
}

TEST(HurchallaFactoringIsPrimeMillerRabin, exhaustive_uint16_t) {
    using T = std::uint16_t;
    for (T m = ut_numeric_limits<T>::max(); m >= 3; m = static_cast<T>(m-2)) {
        SCOPED_TRACE(testing::Message() << "m == " << m);
        MontgomeryForm<T> mf(m);
        EXPECT_TRUE(is_prime_miller_rabin::call_mont(mf)==is_prime_wheel210(m));
    }
}

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 10/20/20.
    TEST(HurchallaFactoringIsPrimeMillerRabin, exhaustive_uint32_t) {
        using T = std::uint32_t;
        for (T m= ut_numeric_limits<T>::max(); m >= 3; m= static_cast<T>(m-2)) {
           MontgomeryForm<T> mf(m);
           EXPECT_TRUE(is_prime_miller_rabin::call_mont(mf) == is_prime_wheel210(m));
        }
    }
#endif


#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic warning "-Wstrict-overflow=2"
#endif
TEST(HurchallaFactoringIsPrimeMillerRabin, basic_test1) {
    using T = std::uint32_t;
    T modulus = 127;

    MontgomeryStandardMathWrapper<T> mWM(modulus);
    MontgomeryForm<T>    mFR(modulus);
    MontgomeryQuarter<T> mQR(modulus);

    EXPECT_TRUE(is_prime_miller_rabin::call_mont(mWM));
    EXPECT_TRUE(is_prime_miller_rabin::call_mont(mFR));
    EXPECT_TRUE(is_prime_miller_rabin::call_mont(mQR));
}
#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif


TEST(HurchallaFactoringIsPrimeMillerRabin, basic_test2) {
    using T = std::uint32_t;
    T modulus = 141;

    MontgomeryStandardMathWrapper<T> mWM(modulus);
    MontgomeryForm<T>    mFR(modulus);
    MontgomeryQuarter<T> mQR(modulus);

    EXPECT_FALSE(is_prime_miller_rabin::call_mont(mWM));
    EXPECT_FALSE(is_prime_miller_rabin::call_mont(mFR));
    EXPECT_FALSE(is_prime_miller_rabin::call_mont(mQR));
}

TEST(HurchallaFactoringIsPrimeMillerRabin, primes_close_to_twoPow64) {
    // Populate a vector of some of the largest primes less than (1<<64).
    // Primes obtained from  https://primes.utm.edu/lists/2small/0bit.html
    using std::uint64_t;
    uint64_t zero = 0;
    // rely on wrap-around when subtracting on next line:
    std::vector<uint64_t> primes = { zero-59, zero-83, zero-95, zero-179,
               zero-189, zero-257, zero-279, zero-323, zero-353, zero-363 };

    std::size_t prime_index = 0;
    for (auto i = zero - 1; i >= primes[9]; i-=2) {
        SCOPED_TRACE(testing::Message() << "i == " << i);
        using T = decltype(i);
        MontgomeryStandardMathWrapper<T> mWM(i);
        MontgomeryForm<T> mFR(i);
        bool is_primeSM = is_prime_miller_rabin::call_mont(mWM);
        bool is_primeFR = is_prime_miller_rabin::call_mont(mFR);
        if (i == primes[prime_index]) {
            EXPECT_TRUE(is_primeSM);
            EXPECT_TRUE(is_primeFR);
            ++prime_index;
        } else {
            EXPECT_FALSE(is_primeSM);
            EXPECT_FALSE(is_primeFR);
        }
    }
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    prime_index = 0;
    for (__uint128_t i = zero - 1; i >= primes[9]; i-=2) {
        using T = __uint128_t;
        MontgomeryStandardMathWrapper<T> mWM(i);
        MontgomeryForm<T>    mFR(i);
        MontgomeryQuarter<T> mQR(i);
        bool is_primeSM = is_prime_miller_rabin::call_mont(mWM);
        bool is_primeFR = is_prime_miller_rabin::call_mont(mFR);
        bool is_primeQR = is_prime_miller_rabin::call_mont(mQR);
        if (i == primes[prime_index]) {
            EXPECT_TRUE(is_primeSM);
            EXPECT_TRUE(is_primeFR);
            EXPECT_TRUE(is_primeQR);
            ++prime_index;
        } else {
            EXPECT_FALSE(is_primeSM);
            EXPECT_FALSE(is_primeFR);
            EXPECT_FALSE(is_primeQR);
        }
    }
#endif
}


#if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(HurchallaFactoringIsPrimeMillerRabin, primes_close_to_twoPow128) {
    // Populate a vector of some of the largest primes less than (1<<128).
    // Primes obtained from  https://primes.utm.edu/lists/2small/100bit.html
    __uint128_t zero = 0;
    // rely on wrap-around when subtracting on next line:
    std::vector<__uint128_t> primes = { zero-159, zero-173, zero-233,
                           zero-237, zero-275, zero-357, zero-675, zero-713,
                           zero-797, zero-1193 };

    std::size_t prime_index = 0;
    for (auto i = zero - 1; i >= primes[9]; i-=2) {
        using T = decltype(i);
        MontgomeryStandardMathWrapper<T> mWM(i);
        MontgomeryForm<T> mFR(i);
        bool is_primeSM = is_prime_miller_rabin::call_mont(mWM);
        bool is_primeFR = is_prime_miller_rabin::call_mont(mFR);
        if (i == primes[prime_index]) {
            EXPECT_TRUE(is_primeSM);
            EXPECT_TRUE(is_primeFR);
            ++prime_index;
        } else {
            EXPECT_FALSE(is_primeSM);
            EXPECT_FALSE(is_primeFR);
        }
    }
}
#endif


} // end namespace
