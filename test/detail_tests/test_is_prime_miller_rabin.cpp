// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

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
    EXPECT_TRUE(is_prime_miller_rabin(mf));
}

TEST(HurchallaFactoringIsPrimeMillerRabin, exhaustive_uint16_t) {
    using T = std::uint16_t;
    for (T m = ut_numeric_limits<T>::max(); m >= 3; m = static_cast<T>(m-2)) {
        SCOPED_TRACE(testing::Message() << "m == " << m);
        MontgomeryForm<T> mf(m);
        EXPECT_TRUE(is_prime_miller_rabin(mf) == is_prime_wheel210(m));
    }
}

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 10/20/20.
    TEST(HurchallaFactoringIsPrimeMillerRabin, exhaustive_uint32_t) {
        using T = std::uint32_t;
        for (T m= ut_numeric_limits<T>::max(); m >= 3; m= static_cast<T>(m-2)) {
            MontgomeryForm<T> mf(m);
            EXPECT_TRUE(is_prime_miller_rabin(mf) == is_prime_wheel210(m));
        }
    }
#endif

TEST(HurchallaFactoringIsPrimeMillerRabin, basic_test1) {
    using T = std::uint32_t;
    T modulus = 127;

    MontgomeryStandardMathWrapper<T> mWM(modulus);
    MontgomeryForm<T>    mFR(modulus);
    MontgomeryQuarter<T> mQR(modulus);

    EXPECT_TRUE(is_prime_miller_rabin(mWM));
    EXPECT_TRUE(is_prime_miller_rabin(mFR));
    EXPECT_TRUE(is_prime_miller_rabin(mQR));
}

TEST(HurchallaFactoringIsPrimeMillerRabin, basic_test2) {
    using T = std::uint32_t;
    T modulus = 141;

    MontgomeryStandardMathWrapper<T> mWM(modulus);
    MontgomeryForm<T>    mFR(modulus);
    MontgomeryQuarter<T> mQR(modulus);

    EXPECT_FALSE(is_prime_miller_rabin(mWM));
    EXPECT_FALSE(is_prime_miller_rabin(mFR));
    EXPECT_FALSE(is_prime_miller_rabin(mQR));
}

TEST(HurchallaFactoringIsPrimeMillerRabin, primes_close_to_twoPow64) {
    // Populate a vector of some of the largest primes less than 2^64.
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
        bool is_primeSM = is_prime_miller_rabin(mWM);
        bool is_primeFR = is_prime_miller_rabin(mFR);
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
        bool is_primeSM = is_prime_miller_rabin(mWM);
        bool is_primeFR = is_prime_miller_rabin(mFR);
        bool is_primeQR = is_prime_miller_rabin(mQR);
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
    // Populate a vector of some of the largest primes less than 2^128.
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
        bool is_primeSM = is_prime_miller_rabin(mWM);
        bool is_primeFR = is_prime_miller_rabin(mFR);
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
