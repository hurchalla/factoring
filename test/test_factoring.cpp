
#include "is_prime_wheel210.h"
#include "hurchalla/factoring/is_prime_monty.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyWrappedStandardMath.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>


namespace {

    TEST(HurchallaFactoringIsPrimeMonty, super_simple_test) {
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;
        uint64_t modulus = 53;
        ma::MontgomeryForm mont(modulus);
        EXPECT_TRUE(hf::is_prime_monty(mont));
    }

    TEST(HurchallaFactoringIsPrimeMonty, exhaustive_uint16_t) {
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;
        using T = uint16_t;
        for (T m=std::numeric_limits<T>::max(); m >= 3; m=static_cast<T>(m-2)) {
            ma::MontgomeryForm mont(m);
            EXPECT_TRUE(hf::is_prime_monty(mont) == hf::is_prime_wheel210(m));
        }
    }

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 6/1/20.
    TEST(HurchallaFactoringIsPrimeMonty, exhaustive_uint32_t) {
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;
        using T = uint32_t;
        for (T m=std::numeric_limits<T>::max(); m >= 3; m=static_cast<T>(m-2)) {
            ma::MontgomeryForm mont(m);
#  if 1
            EXPECT_TRUE(hf::is_prime_monty(mont) == hf::is_prime_wheel210(m));
#  else
            // As of 6/1/20, find_factor_monty() is not yet implemented, and so
            // this section is blocked out with #if.  But once it's available it
            // should be preferred over is_prime_wheel210(), since it should be
            // far more efficient when m is believed to be composite:
            bool is_pm = hf::is_prime_monty(mont);
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

    TEST(HurchallaFactoringIsPrimeMonty, basic_test1) {
        using T = uint32_t;
        T modulus = 127;

        namespace ma = hurchalla::montgomery_arithmetic;
        ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(modulus);
        ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(modulus);
        ma::MontgomeryForm<T, ma::MontyHalfRange<T>> mHR(modulus);
        ma::MontgomeryForm<T, ma::MontyQuarterRange<T>> mQR(modulus);
        ma::MontgomeryForm<T, ma::MontySqrtRange<T>> mSR(modulus);

        namespace hf = hurchalla::factoring;
        EXPECT_TRUE(hf::is_prime_monty(mWSM));
        EXPECT_TRUE(hf::is_prime_monty(mFR));
        EXPECT_TRUE(hf::is_prime_monty(mHR));
        EXPECT_TRUE(hf::is_prime_monty(mQR));
        EXPECT_TRUE(hf::is_prime_monty(mSR));
    }

    TEST(HurchallaFactoringIsPrimeMonty, basic_test2) {
        using T = uint32_t;
        T modulus = 141;

        namespace ma = hurchalla::montgomery_arithmetic;
        ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(modulus);
        ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(modulus);
        ma::MontgomeryForm<T, ma::MontyHalfRange<T>> mHR(modulus);
        ma::MontgomeryForm<T, ma::MontyQuarterRange<T>> mQR(modulus);
        ma::MontgomeryForm<T, ma::MontySqrtRange<T>> mSR(modulus);

        namespace hf = hurchalla::factoring;
        EXPECT_FALSE(hf::is_prime_monty(mWSM));
        EXPECT_FALSE(hf::is_prime_monty(mFR));
        EXPECT_FALSE(hf::is_prime_monty(mHR));
        EXPECT_FALSE(hf::is_prime_monty(mQR));
        EXPECT_FALSE(hf::is_prime_monty(mSR));
    }

    TEST(HurchallaFactoringIsPrimeMonty, primes_close_to_twoPow64) {
        // Populate a vector of some of the largest primes less than 2^64.
        // Primes obtained from  https://primes.utm.edu/lists/2small/0bit.html
        uint64_t zero = (uint64_t)0;
        // rely on wrap-around when subtracting on next line:
        std::vector<uint64_t> primes = { zero-59, zero-83, zero-95, zero-179,
                   zero-189, zero-257, zero-279, zero-323, zero-353, zero-363 };
        namespace ma = hurchalla::montgomery_arithmetic;
        namespace hf = hurchalla::factoring;

        for (auto p : primes) {
            using T = decltype(p);
            ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(p);
            ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(p);
            EXPECT_TRUE(hf::is_prime_monty(mWSM));
            EXPECT_TRUE(hf::is_prime_monty(mFR));
        }
        // hack based on knowing all the primes in the vec differ by at least 10
        for (unsigned int i=2; i<10; i+=2) {
            for (auto p : primes) {
                using T = decltype(p);
                ma::MontgomeryForm<T,ma::MontyWrappedStandardMath<T>> mWSM(p+i);
                ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(p+i);
                EXPECT_FALSE(hf::is_prime_monty(mWSM));
                EXPECT_FALSE(hf::is_prime_monty(mFR));
            }
        }

#ifdef __SIZEOF_INT128__
        for (auto p : primes) {
            using T = __uint128_t;
            ma::MontgomeryForm<T, ma::MontyWrappedStandardMath<T>> mWSM(p);
            ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(p);
            ma::MontgomeryForm<T, ma::MontyHalfRange<T>> mHR(p);
            ma::MontgomeryForm<T, ma::MontyQuarterRange<T>> mQR(p);
            ma::MontgomeryForm<T, ma::MontySqrtRange<T>> mSR(p);
            EXPECT_TRUE(hf::is_prime_monty(mWSM));
            EXPECT_TRUE(hf::is_prime_monty(mFR));
            EXPECT_TRUE(hf::is_prime_monty(mHR));
            EXPECT_TRUE(hf::is_prime_monty(mQR));
            EXPECT_TRUE(hf::is_prime_monty(mSR));
        }
        // hack based on knowing all the primes in the vec differ by at least 10
        for (__uint128_t i=2; i<10; i+=2) {
            for (auto p : primes) {
                using T = __uint128_t;
                ma::MontgomeryForm<T,ma::MontyWrappedStandardMath<T>> mWSM(p+i);
                ma::MontgomeryForm<T, ma::MontyFullRange<T>> mFR(p+i);
                ma::MontgomeryForm<T, ma::MontyHalfRange<T>> mHR(p+i);
                ma::MontgomeryForm<T, ma::MontyQuarterRange<T>> mQR(p+i);
                ma::MontgomeryForm<T, ma::MontySqrtRange<T>> mSR(p+i);
                EXPECT_FALSE(hf::is_prime_monty(mWSM));
                EXPECT_FALSE(hf::is_prime_monty(mFR));
                EXPECT_FALSE(hf::is_prime_monty(mHR));
                EXPECT_FALSE(hf::is_prime_monty(mQR));
                EXPECT_FALSE(hf::is_prime_monty(mSR));
            }
        }
#endif
    }
}

