// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#include "factorize_bruteforce.h"
#include "hurchalla/factoring/factorize.h"
//#include "hurchalla/factoring/detail/wheel_factorization210.h"
//#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <vector>
//#include <array>
#include <algorithm>
#include <numeric>
#include <functional>

namespace {


    TEST(HurchallaFactoringFactorize, exhaustive_uint16_t) {
        namespace hf = hurchalla::factoring;
        using T = std::uint16_t;
        for (T x = std::numeric_limits<T>::max(); x >= 2; --x) {
            std::vector<T> answer = hf::factorize_bruteforce(x);
            std::sort(answer.begin(), answer.end());
            int num_factors;
            auto arr = hf::factorize(x, num_factors);
            std::sort(arr.begin(), arr.begin()+num_factors);
            EXPECT_TRUE(num_factors == static_cast<int>(answer.size()));
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(std::equal(arr.begin(), arr.begin()+num_factors,
                                                               answer.begin()));
        }
    }


// I used this speed test to do a quick and dirty initial performance tuning of
// the value for the macro HURCHALLA_POLLARD_RHO_TRIAL_GCD_THRESHOLD, used in
// PollardRhoTrial.  This test is not needed normally.
    TEST(HurchallaFactoringFactorize, speed_test32) {
        namespace hf = hurchalla::factoring;
        using T = std::uint32_t;
        T max = std::numeric_limits<T>::max()/2;
        for (T x = max; x >= max - 4000000; x = x-2) {
            int num_factors;
            auto arr = hf::factorize(x, num_factors);
            // We need to prevent the compiler from completely removing
            // the factorize calls due to arr never being used.
            // So we'll check arr[0] (which is never 0) just so it's used.
            EXPECT_TRUE(arr[0] != 0);
        }
    }
    TEST(HurchallaFactoringFactorize, speed_test64) {
        namespace hf = hurchalla::factoring;
        using T = std::uint64_t;
        T max = std::numeric_limits<T>::max();
        for (T x = max; x >= max - 200000; x = x-2) {
            int num_factors;
            auto arr = hf::factorize(x, num_factors);
            // We need to prevent the compiler from completely removing
            // the factorize calls due to arr never being used.
            // So we'll check arr[0] (which is never 0) just so it's used.
            EXPECT_TRUE(arr[0] != 0);
        }
    }


    template <typename T>
    T calculate_x(const std::vector<T>& answer)
    {
        return std::accumulate(answer.begin(), answer.end(), static_cast<T>(1),
                                                          std::multiplies<T>());
    }
    template <typename T>
    void test_factorize(const std::vector<T>& answer)
    {
        // multiply all the factors in answer to get the number to factorize.
        T x = calculate_x(answer);
        namespace hf = hurchalla::factoring;
        int num_factors;
        auto arr = hf::factorize(x, num_factors);
        EXPECT_TRUE(num_factors == static_cast<int>(answer.size()));
        // at this time, I haven't made a guarantee for factorize()
        // that the destination range will be sorted, so we'll sort it here.
        std::sort(arr.begin(), arr.begin()+num_factors);
        EXPECT_TRUE(std::equal(arr.begin(), arr.begin()+num_factors,
                                                               answer.begin()));
    }

    TEST(HurchallaFactoringFactorize, hard_semi_primes) {
        using U = std::uint64_t;
        U twoPow32 = static_cast<U>(1) << 32;
        // use largest primes < 2^32:
        // 2^32 minus { 5, 17, 65, 99, 107, 135, 153, 185, 209, 267 }
        std::vector<U> answer = { twoPow32 - 17, twoPow32 - 5 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
        test_factorize(answer);
    }

#ifdef __SIZEOF_INT128__
    TEST(HurchallaFactoringFactorize, hard_semi_primes128_32) {
        using U = __uint128_t;
        U twoPow32 = static_cast<U>(1) << 32;
        // use largest primes < 2^32:
        // 2^32 minus { 5, 17, 65, 99, 107, 135, 153, 185, 209, 267 }
        std::vector<U> answer =
                          { twoPow32-99, twoPow32-65, twoPow32-17, twoPow32-5 };
        test_factorize<U>(answer);
    }

    TEST(HurchallaFactoringFactorize, hard_semi_primes128_42) {
        using U = __uint128_t;
        U twoPow42 = static_cast<U>(1) << 42;
        // use largest primes < 2^42:
        // 2^42 minus { 11, 17, 33, 53, 65, 143, 161, 165, 215, 227 }
        std::vector<U> answer = { twoPow42-33, twoPow42-17, twoPow42-11 };
        test_factorize<U>(answer);
    }

/*
// These tests are fine, but they take multiple seconds or minutes to run.
// Hence we typically will want to skip them.
//
    TEST(HurchallaFactoringFactorize, hard_semi_primes128_52) {
        using U = __uint128_t;
        U twoPow52 = static_cast<U>(1) << 52;
        // use largest primes < 2^52:
        // 2^52 minus { 47, 143, 173, 183, 197, 209, 269, 285, 335, 395 }
        std::vector<U> answer = { twoPow52-143, twoPow52-47 };
        test_factorize<U>(answer);
    }

    TEST(HurchallaFactoringFactorize, hard_semi_primes128_64) {
        using U = __uint128_t;
        U twoPow64 = static_cast<U>(1) << 64;
        // use largest primes < 2^64:
        // 2^64 minus { 59, 83, 95, 179, 189, 257, 279, 323, 353, 363  }
//        std::vector<U> answer = { twoPow64-83, twoPow64-59 };
//        std::vector<U> answer = { twoPow64-179, twoPow64-95 };
        std::vector<U> answer = { twoPow64-257, twoPow64-189 };
        test_factorize<U>(answer);
    }
*/
#endif

    TEST(HurchallaFactoringFactorize, basic_tests) {
        using U = std::uint64_t;
        std::vector<U> answer = { 2, 3, 5, 13, 17 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
        test_factorize<U>(answer);
    }

#ifdef __SIZEOF_INT128__
    TEST(HurchallaFactoringFactorize, basic_tests_128bit) {
        using U = __uint128_t;
        std::vector<U> answer = { 2, 3, 5, 13, 17 };
        test_factorize<U>(answer);
    }
#endif


} // end namespace
