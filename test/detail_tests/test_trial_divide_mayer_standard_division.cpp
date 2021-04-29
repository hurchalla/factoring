// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "hurchalla/factoring/detail/trial_divide_mayer.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>

#include "gtest/gtest.h"

namespace {


#undef HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE
#define HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE 1


using namespace hurchalla::detail;

template <typename T>
void trial_divide_test(T x, T divisor)
{
    EXPECT_TRUE(divisor % 2 == 1);   // required to use trial_divide_mayer

    T div_result;
    bool isdivisible = trial_divide_mayer(div_result, x, divisor);
    EXPECT_TRUE(isdivisible == (x % divisor == 0));
    if (isdivisible) {
        EXPECT_TRUE(div_result == x/divisor);
    }
}

template <typename T>
void trial_divide_typed_tests()
{
    static_assert(!std::numeric_limits<T>::is_signed);
    T midpoint = (static_cast<T>(0) - 1)/2;
    if (midpoint % 2 == 0)
        ++midpoint;

    T midpoint_minus10 = static_cast<T>(midpoint - 10);

    T max = static_cast<T>(static_cast<T>(0) - 1);
    if (max % 2 == 0)
        --max;

    for (T x = 0; x < 20; ++x) {
        for (T n = 1; n < 20; n = static_cast<T>(n+2))
            trial_divide_test(x, n);
        for (T n = max; n > max - 20; n = static_cast<T>(n-2))
            trial_divide_test(x, n);
        for (T n = midpoint_minus10; n < midpoint + 10; n = static_cast<T>(n+2))
            trial_divide_test(x, n);
    }
    for (T x = max; x > max - 20; --x) {
        for (T n = 1; n < 20; n = static_cast<T>(n+2))
            trial_divide_test(x, n);
        for (T n = max; n > max - 20; n = static_cast<T>(n-2))
            trial_divide_test(x, n);
        for (T n = midpoint_minus10; n < midpoint + 10; n = static_cast<T>(n+2))
            trial_divide_test(x, n);
    }
    for (T x = midpoint_minus10; x < midpoint + 10; ++x) {
        for (T n = 1; n < 20; n = static_cast<T>(n+2))
            trial_divide_test(x, n);
        for (T n = max; n > max - 20; n = static_cast<T>(n-2))
            trial_divide_test(x, n);
        for (T n = midpoint_minus10; n < midpoint + 10; n = static_cast<T>(n+2))
            trial_divide_test(x, n);
    }
}

TEST(HurchallaFactoringTrialDivideMayer, mayer_tests_standard_division) {
    trial_divide_typed_tests<std::uint8_t>();
    trial_divide_typed_tests<std::uint16_t>();
    trial_divide_typed_tests<std::uint32_t>();
    trial_divide_typed_tests<std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    trial_divide_typed_tests<__uint128_t>();
#endif
}


} // end namespace
