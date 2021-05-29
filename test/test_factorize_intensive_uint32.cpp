// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "hurchalla/factoring/resource_intensive_api/factorize_intensive_uint32.h"
#include "hurchalla/factoring/resource_intensive_api/IsPrimeIntensive.h"
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

#if 0
// I used this speed test to do a quick and dirty initial performance tuning of
// PollardRhoTrial and PollardRhoBrentTrial.  This test is not needed normally.
    const IsPrimeIntensive<std::uint32_t,true>& get_ipi32()
    {
        static IsPrimeIntensive<std::uint32_t,true> ipi;
        return ipi;
    }
    TEST(HurchallaFactoringFactorizeIntensiveUint32, speed_test_primer) {
        get_ipi32();
    }
    TEST(HurchallaFactoringFactorizeIntensiveUint32, speed_test32_intensive) {
        const auto& ipi = get_ipi32();
        using T = std::uint32_t;
        T max = ut_numeric_limits<T>::max()/2;
        for (T x = max; x >= max - 4000000; x = x-2) {
            int num_factors;
            auto arr = factorize_intensive_uint32(x, num_factors, ipi);
            // We need to prevent the compiler from completely removing
            // the factorize calls due to arr never being used.
            // So we'll check arr[0] (which is never 0) just so it's used.
            EXPECT_TRUE(arr[0] != 0);
        }
    }
#endif


// this is a bit of a hack, but the test is extremely slow without optimization
// and so we try to skip it when not building with any optimizations
#if (defined(_MSC_VER) && !defined(_DEBUG)) || defined(__OPTIMIZE__)

template <typename T>
T calculate_x(const std::vector<T>& answer)
{
    return std::accumulate(answer.begin(), answer.end(), static_cast<T>(1),
                                                      std::multiplies<T>());
}

template <typename T>
void test_factorize(const std::vector<T>& answer,
                    const IsPrimeIntensive<std::uint32_t,true>& ipi)
{
    // multiply all the factors in answer to get the number to factorize.
    T x = calculate_x(answer);
    int num_factors;
    auto arr = factorize_intensive_uint32(x, num_factors, ipi);
    EXPECT_TRUE(num_factors == static_cast<int>(answer.size()));
    // at this time, I haven't made a guarantee for factorize_intensive_uint32()
    // that the destination range will be sorted, so we'll sort it here.
    std::sort(arr.begin(), arr.begin()+num_factors);
    EXPECT_TRUE(std::equal(arr.begin(), arr.begin()+num_factors,
                                                           answer.begin()));
}

TEST(HurchallaFactoringFactorizeIntensiveUint32, basic_tests_and_hard_semiprimes) {
    using U = std::uint32_t;
    IsPrimeIntensive<std::uint32_t,true> ipi;

    // basic test
    std::vector<U> answer1 = { 2, 3, 5, 13, 17 };
    SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer1));
    test_factorize<U>(answer1, ipi);

    // hard semiprimes
    U twoPow16 = static_cast<U>(1) << 16;
    // use largest primes < 2^16:
    // 2^16 minus { 15, 17, 39, 57, 87, 89, 99, 113, 117, 123 }
    std::vector<U> answer2 = { twoPow16 - 17, twoPow16 - 15 };
    SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer2));
    test_factorize(answer2, ipi);
}

#endif

} // end namespace
