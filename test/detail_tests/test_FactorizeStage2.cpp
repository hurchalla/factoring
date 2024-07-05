// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// FactorizeStage2.h by default uses PollardRhoBrentSwitchingTrial for its
// Pollard-Rho trials.
// If for some reason you instead want to test a different Pollard-Rho trial
// (typically from the detail/experimental folder), then uncomment the #define
// below and specify the name of the Pollard-Rho trial that you want to use.
// As an example, the commented-out #define below specifies the name
// PollardRhoTrial.  Ordinarily there's no reason for you to test any of the
// experimental trials.

//#define HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME PollardRhoTrial


#include "../factorize_bruteforce.h"
#include "hurchalla/factoring/detail/FactorizeStage2.h"
#include "hurchalla/factoring/detail/impl_factorize.h"
#include "hurchalla/factoring/detail/IsPrimeFactor.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>

namespace {


using namespace hurchalla;
using namespace hurchalla::detail;


TEST(HurchallaFactoringFactorizePollardRho, exhaustive_uint16_t) {
    using T = std::uint16_t;
    T max = ut_numeric_limits<T>::max();
    if (max % 2 == 0)
        --max;
    for (T x = max; x >= 3; x = static_cast<T>(x-2)) {
        std::vector<T> answer = factorize_bruteforce(x);
        std::sort(answer.begin(), answer.end());
        std::vector<T> vec;

#ifndef HURCHALLA_FACTORING_ECM_THRESHOLD_BITS
# error "HURCHALLA_FACTORING_ECM_THRESHOLD_BITS must be defined"
#endif
        constexpr int EcmMinBits = HURCHALLA_FACTORING_ECM_THRESHOLD_BITS;
        constexpr int MaxBitsX = ut_numeric_limits<T>::digits;
        auto is_prime_functor = IsPrimeFactor();
        T always_prime_limit = 0;
        bool expect_arbitrary_size_factors = true;
        FactorizeStage2<EcmMinBits, MaxBitsX, T>
                                            factorize_stage2(always_prime_limit,
                                                 expect_arbitrary_size_factors);
        factorize_stage2(std::back_inserter(vec), is_prime_functor, x);

        std::sort(vec.begin(), vec.end());
        SCOPED_TRACE(testing::Message() << "x == " << x);
        EXPECT_TRUE(vec.size() == answer.size());
        EXPECT_TRUE(std::equal(vec.begin(), vec.end(), answer.begin()));
    }
}


template <typename T>
T calculate_x(const std::vector<T>& answer)
{
    return std::accumulate(answer.begin(), answer.end(), static_cast<T>(1),
                                                      std::multiplies<T>());
}
template <typename T>
void test_factorize(const std::vector<T>& answer,
                    bool expect_arbitrary_size_factors)
{
    // multiply all the factors in answer to get the number to factorize.
    T x = calculate_x(answer);
    EXPECT_TRUE(x > 0);  // the test is bad if this fails
    std::vector<T> vec;

//    constexpr bool allowHalfRange = ut_numeric_limits<T>::is_signed;
    using U = typename extensible_make_unsigned<T>::type;

    //convenient adapter to correctly cast between (possible)signed and unsigned
    struct FactorVectorAdapter {
        using value_type = U;
        explicit FactorVectorAdapter(std::vector<T>& a) : v(a) {}
        void push_back(const value_type& val)
        {
            HPBC_ASSERT2(val <= static_cast<U>(ut_numeric_limits<T>::max()));
            v.push_back(static_cast<T>(val));
        }
    private:
        std::vector<T>& v;
    };
    FactorVectorAdapter fva(vec);

    auto iter = std::back_inserter(fva);
    while (x % 2 == 0) {
        *iter++ = 2;
        x = x/2;
    }

#ifndef HURCHALLA_FACTORING_ECM_THRESHOLD_BITS
# error "HURCHALLA_FACTORING_ECM_THRESHOLD_BITS must be defined"
#endif
    constexpr int EcmMinBits = HURCHALLA_FACTORING_ECM_THRESHOLD_BITS;
    constexpr int MaxBitsX = ut_numeric_limits<T>::digits;
    auto is_prime_functor = IsPrimeFactor();
    U always_prime_limit = 0;
    FactorizeStage2<EcmMinBits, MaxBitsX, U>
                    factorize_stage2(always_prime_limit,
                                                 expect_arbitrary_size_factors);
    factorize_stage2(iter, is_prime_functor, static_cast<U>(x));


    EXPECT_TRUE(vec.size() == answer.size());
    // at this time, I haven't made a guarantee for factorize()
    // that the destination range will be sorted, so we'll sort it here.
    std::sort(vec.begin(), vec.end());
    EXPECT_TRUE(std::equal(vec.begin(), vec.end(), answer.begin()));
}

TEST(HurchallaFactoringFactorizePollardRho, hard_semi_primes) {
    using U = std::uint64_t;
    U twoPow32 = static_cast<U>(1) << 32;
    // use largest primes < (1<<32):
    // (1<<32) minus { 5, 17, 65, 99, 107, 135, 153, 185, 209, 267 }
    std::vector<U> answer = { twoPow32 - 99, twoPow32 - 65 };
    SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
    test_factorize(answer, true);
    test_factorize(answer, false);
}

TEST(HurchallaFactoringFactorizePollardRho, signed_hard_semi_primes32) {
    using T = std::int32_t;
    T twoPow15 = static_cast<T>(1) << 15;
    // use largest primes < (1<<15):
    // (1<<15) minus { 19, 49, 51, 55, 61, 75, 81, 115, 121, 135 }
    std::vector<T> answer = { twoPow15 - 49, twoPow15 - 19 };
    SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
    test_factorize(answer, true);
    test_factorize(answer, false);
}
TEST(HurchallaFactoringFactorizePollardRho, signed_hard_semi_primes64) {
    using T = std::int64_t;
    T twoPow31 = static_cast<T>(1) << 31;
    // use largest primes < (1<<31):
    // (1<<31) minus { 1, 19, 61, 69, 85, 99, 105, 151, 159, 171 }
    std::vector<T> answer = { twoPow31 - 19, twoPow31 - 1 };
    SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
    test_factorize(answer, true);
    test_factorize(answer, false);
}
#if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(HurchallaFactoringFactorizePollardRho, signed_hard_semi_primes128) {
    using T = __int128_t;
    T twoPow33 = static_cast<T>(1) << 33;
    // use largest primes < (1<<33):
    // (1<<33) minus { 9, 25, 49, 79, 105, 285, 301, 303, 321, 355 }
    std::vector<T> answer = { twoPow33 - 25, twoPow33 - 9 };
    test_factorize(answer, true);
}
#endif

#if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(HurchallaFactoringFactorizePollardRho, hard_semi_primes128_32) {
    using U = __uint128_t;
    U twoPow32 = static_cast<U>(1) << 32;
    // use largest primes < (1<<32):
    // (1<<32) minus { 5, 17, 65, 99, 107, 135, 153, 185, 209, 267 }
    std::vector<U> answer =
                     { twoPow32-185, twoPow32-153, twoPow32-135, twoPow32-107 };
    test_factorize<U>(answer, true);
}
#endif


TEST(HurchallaFactoringFactorizePollardRho, basic_tests) {
    {
        using U = std::uint64_t;
        std::vector<U> answer = { 3, 5, 19, 23, 59, 127 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
        test_factorize<U>(answer, true);
        test_factorize<U>(answer, false);
    }
    {
        using U = std::uint32_t;
        std::vector<U> answer = { 2, 2, 2, 43, 59, 59, 113 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
        test_factorize<U>(answer, true);
        test_factorize<U>(answer, false);
    }
    {
        using U = std::uint32_t;
        std::vector<U> answer = { 32771, 32771 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
        test_factorize<U>(answer, true);
        test_factorize<U>(answer, false);
    }
}


#if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(HurchallaFactoringFactorizePollardRho, basic_tests_128bit) {
    using U = __uint128_t;
    std::vector<U> answer = { 2, 3, 5, 13, 17 };
    test_factorize<U>(answer, true);
    test_factorize<U>(answer, false);
}
#endif


} // end namespace
