// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/factorize_trialdivision.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <vector>
#include <numeric>
#include <functional>
#include <algorithm>


#include "gtest/gtest.h"

namespace {


using namespace hurchalla;
using namespace hurchalla::detail;


template <typename T>
T calculate_x(const std::vector<T>& answer)
{
    return std::accumulate(answer.begin(), answer.end(), static_cast<T>(1),
                                                      std::multiplies<T>());
}


template <template<class,int> class TTD, int SIZE, typename T>
void ftd_test(std::vector<T>& sorted_answer)
{
    // multiply all the factors in sorted_answer to get the number to factorize.
    T x = calculate_x(sorted_answer);

    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    EXPECT_TRUE(x >= 2);  // 0 and 1 do not have prime factorizations.

    T q, nextprime;
    std::vector<T> vec;
    factorize_trialdivision::call<TTD,SIZE>(std::back_inserter(vec), q,
                                                                  nextprime, x);
    std::sort(vec.begin(), vec.end());
    for (size_t i=0; i<sorted_answer.size(); ++i) {
        auto n = sorted_answer[i];
        if (n < nextprime) {
            EXPECT_TRUE(i < vec.size());
            EXPECT_TRUE(vec[i] == n);
        }
        else
            break;
    }
}


template <template<class,int> class TTD, int SIZE, typename T>
void ftd_sized_tests()
{
    static_assert(1 < SIZE);
    static_assert(!ut_numeric_limits<T>::is_signed);
    {
        std::vector<T> answer = { 2, 3 };
        ftd_test<TTD, SIZE>(answer);
    }
    {
        std::vector<T> answer = { 5, 5, 7 };
        ftd_test<TTD, SIZE>(answer);
    }
    if constexpr (ut_numeric_limits<T>::digits >= 16)
    {
        {
            std::vector<T> answer = { 7, 7, 11, 13 };
            ftd_test<TTD, SIZE>(answer);
        }
        {
            std::vector<T> answer = { 7, 31, 31 };
            ftd_test<TTD, SIZE>(answer);
        }
        {
            std::vector<T> answer = { 251 };
            ftd_test<TTD, SIZE>(answer);
        }
    }
    if constexpr (ut_numeric_limits<T>::digits >= 32)
    {
        {
            std::vector<T> answer = { 31, 257, 257 };
            ftd_test<TTD, SIZE>(answer);
        }
        {
            std::vector<T> answer = { 17, 65537 };
            ftd_test<TTD, SIZE>(answer);
        }
        {
            std::vector<T> answer = { 2, 3, 5, 7, 11, 13, 17, 19, 29 };
            ftd_test<TTD, SIZE>(answer);
        }
    }
    if constexpr (ut_numeric_limits<T>::digits >= 64)
    {
        {
            std::vector<T> answer = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 65537 };
            ftd_test<TTD, SIZE>(answer);
        }
        {
            std::vector<T> answer = { 2, 3, 5, 5, 65537, 65537 };
            ftd_test<TTD, SIZE>(answer);
        }
    }
}


template <template<class,int> class TTD, typename T>
void ftd_typed_tests()
{
    ftd_sized_tests<TTD, 2, T>();
    ftd_sized_tests<TTD, 10, T>();
    ftd_sized_tests<TTD, 54, T>();
    if constexpr (ut_numeric_limits<T>::digits >= 16) {
        ftd_sized_tests<TTD, 55, T>();
        ftd_sized_tests<TTD, 198, T>();
        ftd_sized_tests<TTD, 1000, T>();
    }
    // Using SIZE as large as 2501 (and even 1000 above) faces some risk of
    // getting a compiler failure of
    // "constexpr evaluation hit maximum step limit; possible infinite loop?"
    // MSVC 2017 doesn't even compile it, giving an error that seems to have
    // nothing to do with the real problem of exceeding a size limit.
    // It's not a bug in any of the code under test or the test itself, it's
    // just a limit of the compiler with constexpr initialization.  There's
    // probably a compiler flag you can set to increase the limit, but if
    // you run into this error I'd simply suggest commenting out the test 
    // here that caused it - even the test with SIZE of 1000 is far larger
    // than I expect we'd ever use in practice for trial division.  I'd only
    // try to fix such an error if it happens with SIZE < 200, since that's
    // a more practical range (though ~200 is still fairly large).
#if !defined(_MSC_VER) || (_MSC_VER >= 1927)
    if constexpr (ut_numeric_limits<T>::digits >= 32) {
        ftd_sized_tests<TTD, 2501, T>();
    }
#endif
}


template <template<class,int> class TTD>
void ftd_ttd_tests()
{
    ftd_typed_tests<TTD, std::uint8_t>();
    ftd_typed_tests<TTD, std::uint16_t>();
    ftd_typed_tests<TTD, std::uint32_t>();
    ftd_typed_tests<TTD, std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    ftd_typed_tests<TTD, __uint128_t>();
#endif
}


TEST(HurchallaFactoringFactorizeTrialDivision, basic_tests) {
    ftd_ttd_tests<PrimeTrialDivisionWarren>();
    ftd_ttd_tests<PrimeTrialDivisionMayer>();
}


} // end namespace
