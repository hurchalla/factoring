
#include "is_prime_wheel210.h"
#include "hurchalla/factoring/is_prime.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>

#include <iostream>

namespace {

    TEST(HurchallaFactoringIsPrime, exhaustive_uint16_t) {
        namespace hf = hurchalla::factoring;
        using T = uint16_t;
        for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
            EXPECT_TRUE(hf::is_prime(x) == hf::is_prime_wheel210(x));
        }
        T x = std::numeric_limits<T>::max();
        EXPECT_TRUE(hf::is_prime(x) == hf::is_prime_wheel210(x));
    }

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 6/7/20.
    TEST(HurchallaFactoringIsPrime, exhaustive_uint32_t) {
        namespace hf = hurchalla::factoring;
        using T = uint32_t;
        for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
            EXPECT_TRUE(hf::is_prime(x) == hf::is_prime_wheel210(x));
        }
        T x = std::numeric_limits<T>::max();
        EXPECT_TRUE(hf::is_prime(x) == hf::is_prime_wheel210(x));
    }
#endif


    TEST(HurchallaFactoringIsPrimeMillerRabin, super_simple_test) {
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;
        uint64_t modulus = 53;
        ma::MontgomeryForm<uint64_t> mont(modulus);
        EXPECT_TRUE(hf::is_prime_miller_rabin(mont));
    }

    TEST(HurchallaFactoringIsPrimeMillerRabin, exhaustive_uint16_t) {
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;
        using T = uint16_t;
        for (T m=std::numeric_limits<T>::max(); m >= 3; m=static_cast<T>(m-2)) {
            ma::MontgomeryForm<T> mont(m);
            EXPECT_TRUE(hf::is_prime_miller_rabin(mont) ==
                                                      hf::is_prime_wheel210(m));
        }
    }

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 6/1/20.
    TEST(HurchallaFactoringIsPrimeMillerRabin, exhaustive_uint32_t) {
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;
        using T = uint32_t;
        for (T m=std::numeric_limits<T>::max(); m >= 3; m=static_cast<T>(m-2)) {
            ma::MontgomeryForm<T> mont(m);
#  if 1
            EXPECT_TRUE(hf::is_prime_miller_rabin(mont) ==
                                                      hf::is_prime_wheel210(m));
#  else
            // As of 6/1/20, find_factor_monty() is not yet implemented, and so
            // this section is blocked out with #if.  But once it's available it
            // should be preferred over is_prime_wheel210(), since it should be
            // far more efficient when m is believed to be composite:
            bool is_pm = hf::is_prime_miller_rabin(mont);
            if (is_pm)
                EXPECT_TRUE(hf::is_prime_wheel210(m));
            else {
                int max_trials = 1000;
                // find_factor_monty() doesn't absolutely guarantee we will find
                // a factor, but as max_trials becomes large, the probability
                // becomes almost a certainty (if m is composite).
                T f = hf::find_factor_monty(mont, max_trials);
                // assume find_factor_monty()'s return is > 0 for composite m.
                EXPECT_TRUE(f != 0 && m % f == 0) << "f==" << f <<
                                                     ", m==" << m << "\n";
            }
#  endif
        }
    }
#endif

    TEST(HurchallaFactoringIsPrimeMillerRabin, basic_test1) {
        using T = uint32_t;
        T modulus = 127;

        namespace ma = hurchalla::montgomery_arithmetic;
        ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(modulus);
        ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(modulus);
        ma::MontgomeryForm<T, ma::MontyHalfRange<T>> mHR(modulus);
        ma::MontgomeryForm<T, ma::MontyQuarterRange<T>> mQR(modulus);
        ma::MontgomeryForm<T, ma::MontySqrtRange<T>> mSR(modulus);

        namespace hf = hurchalla::factoring;
        EXPECT_TRUE(hf::is_prime_miller_rabin(mWSM));
        EXPECT_TRUE(hf::is_prime_miller_rabin(mFR));
        EXPECT_TRUE(hf::is_prime_miller_rabin(mHR));
        EXPECT_TRUE(hf::is_prime_miller_rabin(mQR));
        EXPECT_TRUE(hf::is_prime_miller_rabin(mSR));
    }

    TEST(HurchallaFactoringIsPrimeMillerRabin, basic_test2) {
        using T = uint32_t;
        T modulus = 141;

        namespace ma = hurchalla::montgomery_arithmetic;
        ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(modulus);
        ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(modulus);
        ma::MontgomeryForm<T, ma::MontyHalfRange<T>> mHR(modulus);
        ma::MontgomeryForm<T, ma::MontyQuarterRange<T>> mQR(modulus);
        ma::MontgomeryForm<T, ma::MontySqrtRange<T>> mSR(modulus);

        namespace hf = hurchalla::factoring;
        EXPECT_FALSE(hf::is_prime_miller_rabin(mWSM));
        EXPECT_FALSE(hf::is_prime_miller_rabin(mFR));
        EXPECT_FALSE(hf::is_prime_miller_rabin(mHR));
        EXPECT_FALSE(hf::is_prime_miller_rabin(mQR));
        EXPECT_FALSE(hf::is_prime_miller_rabin(mSR));
    }

    TEST(HurchallaFactoringIsPrimeMillerRabin, primes_close_to_twoPow64) {
        // Populate a vector of some of the largest primes less than 2^64.
        // Primes obtained from  https://primes.utm.edu/lists/2small/0bit.html
        uint64_t zero = (uint64_t)0;
        // rely on wrap-around when subtracting on next line:
        std::vector<uint64_t> primes = { zero-59, zero-83, zero-95, zero-179,
                   zero-189, zero-257, zero-279, zero-323, zero-353, zero-363 };
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;

        size_t prime_index = 0;
        for (auto i = zero - 1; i >= primes[9]; i-=2) {
            using T = decltype(i);
            ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(i);
            ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(i);
            bool is_primeWSM = hf::is_prime_miller_rabin(mWSM);
            bool is_primeFR = hf::is_prime_miller_rabin(mFR);
            if (i == primes[prime_index]) {
                EXPECT_TRUE(is_primeWSM);
                EXPECT_TRUE(is_primeFR);
                ++prime_index;
            } else {
                EXPECT_FALSE(is_primeWSM);
                EXPECT_FALSE(is_primeFR);
            }
        }

#ifdef __SIZEOF_INT128__
        prime_index = 0;
        for (__uint128_t i = zero - 1; i >= primes[9]; i-=2) {
            using T = __uint128_t;
            ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(i);
            ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(i);
            ma::MontgomeryForm<T, ma::MontyHalfRange<T>> mHR(i);
            ma::MontgomeryForm<T, ma::MontyQuarterRange<T>> mQR(i);
            ma::MontgomeryForm<T, ma::MontySqrtRange<T>> mSR(i);
            bool is_primeWSM = hf::is_prime_miller_rabin(mWSM);
            bool is_primeFR = hf::is_prime_miller_rabin(mFR);
            bool is_primeHR = hf::is_prime_miller_rabin(mHR);
            bool is_primeQR = hf::is_prime_miller_rabin(mQR);
            bool is_primeSR = hf::is_prime_miller_rabin(mSR);
            if (i == primes[prime_index]) {
                EXPECT_TRUE(is_primeWSM);
                EXPECT_TRUE(is_primeFR);
                EXPECT_TRUE(is_primeHR);
                EXPECT_TRUE(is_primeQR);
                EXPECT_TRUE(is_primeSR);
                ++prime_index;
            } else {
                EXPECT_FALSE(is_primeWSM);
                EXPECT_FALSE(is_primeFR);
                EXPECT_FALSE(is_primeHR);
                EXPECT_FALSE(is_primeQR);
                EXPECT_FALSE(is_primeSR);
            }
        }
#endif
    }


#ifdef __SIZEOF_INT128__
    TEST(HurchallaFactoringIsPrimeMillerRabin, primes_close_to_twoPow128) {
        // Populate a vector of some of the largest primes less than 2^128.
        // Primes obtained from  https://primes.utm.edu/lists/2small/100bit.html
        __uint128_t zero = (__uint128_t)0;
        // rely on wrap-around when subtracting on next line: 
        std::vector<__uint128_t> primes = { zero-159, zero-173, zero-233,
                               zero-237, zero-275, zero-357, zero-675, zero-713,
                               zero-797, zero-1193 };
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;

        size_t prime_index = 0;
        for (auto i = zero - 1; i >= primes[9]; i-=2) {
            using T = decltype(i);
            ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(i);
            ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(i);
            bool is_primeWSM = hf::is_prime_miller_rabin(mWSM);
            bool is_primeFR = hf::is_prime_miller_rabin(mFR);
            if (i == primes[prime_index]) {
                EXPECT_TRUE(is_primeWSM);
                EXPECT_TRUE(is_primeFR);
                ++prime_index;
            } else {
                EXPECT_FALSE(is_primeWSM);
                EXPECT_FALSE(is_primeFR);
            }
        }
    }
#endif


} // end namespace
