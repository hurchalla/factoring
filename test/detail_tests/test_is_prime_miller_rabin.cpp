// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/factoring/detail/is_prime_wheel210.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <cstddef>

namespace {

    TEST(HurchallaFactoringIsPrimeMillerRabin, super_simple_test) {
        namespace hc = hurchalla;
        namespace hd = hurchalla::detail;
        std::uint64_t modulus = 53;
        hc::MontgomeryForm<std::uint64_t> mf(modulus);
        EXPECT_TRUE(hd::is_prime_miller_rabin(mf));
    }

    TEST(HurchallaFactoringIsPrimeMillerRabin, exhaustive_uint16_t) {
        namespace hc = hurchalla;
        namespace hd = hurchalla::detail;
        using T = std::uint16_t;
        for (T m=std::numeric_limits<T>::max(); m >= 3; m=static_cast<T>(m-2)) {
            SCOPED_TRACE(testing::Message() << "m == " << m);
            hc::MontgomeryForm<T> mf(m);
            EXPECT_TRUE(hd::is_prime_miller_rabin(mf) ==
                                                      hd::is_prime_wheel210(m));
        }
    }

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 10/20/20.
    TEST(HurchallaFactoringIsPrimeMillerRabin, exhaustive_uint32_t) {
        namespace hc = hurchalla;
        namespace hd = hurchalla::detail;
        using T = std::uint32_t;
        for (T m=std::numeric_limits<T>::max(); m >= 3; m=static_cast<T>(m-2)) {
            hc::MontgomeryForm<T> mf(m);
#  if 1
            EXPECT_TRUE(hd::is_prime_miller_rabin(mf) ==
                                                     hd::is_prime_wheel210(m));
#  else
            // As of 10/20/20, find_factor_pollard_rho() is not yet implemented,
            // so this section is blocked out with #if.  But once it's available
            // it should be preferred over is_prime_wheel210(), since it should
            // be far more efficient when m is believed to be composite:
            bool is_pm = hd::is_prime_miller_rabin(mf);
            if (is_pm)
                EXPECT_TRUE(hd::is_prime_wheel210(m));
            else {
                int max_trials = 100;
                // find_factor_pollard_rho() doesn't absolutely guarantee we
                // will find a factor, but as max_trials becomes large, the
                // probability becomes almost a certainty (if m is composite).
                T f = hd::find_factor_pollardrho(mf, max_trials);
                // assume the returned factor is > 1 for composite m.
                EXPECT_TRUE(f > 1 && m % f == 0) << "f==" << f <<
                                                     ", m==" << m << "\n";
            }
#  endif
        }
    }
#endif

    TEST(HurchallaFactoringIsPrimeMillerRabin, basic_test1) {
        using T = std::uint32_t;
        T modulus = 127;

        namespace hc = hurchalla;
        hc::MontgomeryStandardMathWrapper<T> mWM(modulus);
        hc::MontgomeryForm<T>    mFR(modulus);
        hc::MontgomeryQuarter<T> mQR(modulus);

        namespace hd = hurchalla::detail;
        EXPECT_TRUE(hd::is_prime_miller_rabin(mWM));
        EXPECT_TRUE(hd::is_prime_miller_rabin(mFR));
        EXPECT_TRUE(hd::is_prime_miller_rabin(mQR));
    }

    TEST(HurchallaFactoringIsPrimeMillerRabin, basic_test2) {
        using T = std::uint32_t;
        T modulus = 141;

        namespace hc = hurchalla;
        hc::MontgomeryStandardMathWrapper<T> mWM(modulus);
        hc::MontgomeryForm<T>    mFR(modulus);
        hc::MontgomeryQuarter<T> mQR(modulus);

        namespace hd = hurchalla::detail;
        EXPECT_FALSE(hd::is_prime_miller_rabin(mWM));
        EXPECT_FALSE(hd::is_prime_miller_rabin(mFR));
        EXPECT_FALSE(hd::is_prime_miller_rabin(mQR));
    }

    TEST(HurchallaFactoringIsPrimeMillerRabin, primes_close_to_twoPow64) {
        // Populate a vector of some of the largest primes less than 2^64.
        // Primes obtained from  https://primes.utm.edu/lists/2small/0bit.html
        using std::uint64_t;
        uint64_t zero = 0;
        // rely on wrap-around when subtracting on next line:
        std::vector<uint64_t> primes = { zero-59, zero-83, zero-95, zero-179,
                   zero-189, zero-257, zero-279, zero-323, zero-353, zero-363 };
        namespace hc = hurchalla;
        namespace hd = hurchalla::detail;

        std::size_t prime_index = 0;
        for (auto i = zero - 1; i >= primes[9]; i-=2) {
            SCOPED_TRACE(testing::Message() << "i == " << i);
            using T = decltype(i);
            hc::MontgomeryStandardMathWrapper<T> mWM(i);
            hc::MontgomeryForm<T> mFR(i);
            bool is_primeSM = hd::is_prime_miller_rabin(mWM);
            bool is_primeFR = hd::is_prime_miller_rabin(mFR);
            if (i == primes[prime_index]) {
                EXPECT_TRUE(is_primeSM);
                EXPECT_TRUE(is_primeFR);
                ++prime_index;
            } else {
                EXPECT_FALSE(is_primeSM);
                EXPECT_FALSE(is_primeFR);
            }
        }

#ifdef __SIZEOF_INT128__
        prime_index = 0;
        for (__uint128_t i = zero - 1; i >= primes[9]; i-=2) {
            using T = __uint128_t;
            hc::MontgomeryStandardMathWrapper<T> mWM(i);
            hc::MontgomeryForm<T>    mFR(i);
            hc::MontgomeryQuarter<T> mQR(i);
            bool is_primeSM = hd::is_prime_miller_rabin(mWM);
            bool is_primeFR = hd::is_prime_miller_rabin(mFR);
            bool is_primeQR = hd::is_prime_miller_rabin(mQR);
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


#ifdef __SIZEOF_INT128__
    TEST(HurchallaFactoringIsPrimeMillerRabin, primes_close_to_twoPow128) {
        // Populate a vector of some of the largest primes less than 2^128.
        // Primes obtained from  https://primes.utm.edu/lists/2small/100bit.html
        __uint128_t zero = 0;
        // rely on wrap-around when subtracting on next line: 
        std::vector<__uint128_t> primes = { zero-159, zero-173, zero-233,
                               zero-237, zero-275, zero-357, zero-675, zero-713,
                               zero-797, zero-1193 };
        namespace hc = hurchalla;
        namespace hd = hurchalla::detail;

        std::size_t prime_index = 0;
        for (auto i = zero - 1; i >= primes[9]; i-=2) {
            using T = decltype(i);
            hc::MontgomeryStandardMathWrapper<T> mWM(i);
            hc::MontgomeryForm<T> mFR(i);
            bool is_primeSM = hd::is_prime_miller_rabin(mWM);
            bool is_primeFR = hd::is_prime_miller_rabin(mFR);
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
