// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


// factorize_pollard_rho.h by default uses PollardRhoBrentTrial.
// If for some reason you instead want to test one of the experimental Pollard-
// Rho trials in the detail/experimental folder, then uncomment the #define
// below and specify the name of the experimental trial that you want to use.
// As an example, the commented-out #define below specifies the name
// PollardRhoTrial.  Ordinarily there's no reason for you to test any of the
// experimental trials.

//#define HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME PollardRhoTrial


#include "../factorize_bruteforce.h"
#include "hurchalla/factoring/detail/factorize_pollard_rho.h"
#include "hurchalla/factoring/detail/PollardRhoIsPrime.h"
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
        factorize_pollard_rho::call(std::back_inserter(vec), x,
                                                           PollardRhoIsPrime());
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
void test_factorize(const std::vector<T>& answer)
{
    // multiply all the factors in answer to get the number to factorize.
    T x = calculate_x(answer);
    EXPECT_TRUE(x > 0);  // the test is bad if this fails
    std::vector<T> vec;
    auto iter = std::back_inserter(vec);
    while (x % 2 == 0) {
        *iter++ = 2;
        x = x/2;
    }
    factorize_pollard_rho::call(iter, x, PollardRhoIsPrime());
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
    test_factorize(answer);
}


#if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(HurchallaFactoringFactorizePollardRho, hard_semi_primes128_32) {
    using U = __uint128_t;
    U twoPow32 = static_cast<U>(1) << 32;
    // use largest primes < (1<<32):
    // (1<<32) minus { 5, 17, 65, 99, 107, 135, 153, 185, 209, 267 }
    std::vector<U> answer =
                     { twoPow32-185, twoPow32-153, twoPow32-135, twoPow32-107 };
    test_factorize<U>(answer);
}
#endif


TEST(HurchallaFactoringFactorizePollardRho, basic_tests) {
    {
        using U = std::uint64_t;
        std::vector<U> answer = { 3, 5, 19, 23, 59, 127 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
        test_factorize<U>(answer);
    }
    {
        using U = std::uint32_t;
        std::vector<U> answer = { 2, 2, 2, 43, 59, 59, 113 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
        test_factorize<U>(answer);
    }
    {
        using U = std::uint32_t;
        std::vector<U> answer = { 32771, 32771 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer));
        test_factorize<U>(answer);
    }
}


#if HURCHALLA_COMPILER_HAS_UINT128_T()
TEST(HurchallaFactoringFactorizePollardRho, basic_tests_128bit) {
    using U = __uint128_t;
    std::vector<U> answer = { 2, 3, 5, 13, 17 };
    test_factorize<U>(answer);
}
#endif


} // end namespace
