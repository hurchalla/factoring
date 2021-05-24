// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla


#undef HURCHALLA_PREFER_EUCLIDEAN_GCD
#undef HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE

// test_greatest_common_divisor_euclid.cpp defines the following macro,
// and that file has nothing else in it except a #include of this file
#ifdef TESTING_GREATEST_COMMON_DIVISOR_EUCLID
   // in order to test the Euclidean gcd, rather than the normal Binary gcd,
   // these two macros need to be defined:
#  define HURCHALLA_PREFER_EUCLIDEAN_GCD 1
#  define HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE 1
#endif


#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>

#include "gtest/gtest.h"

namespace {


using namespace hurchalla;

template <typename T>
void test_gcd()
{
    T a = 6; T b = 8;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(2));
    a = 110; b = 121;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(11));
    a = 210; b = 150;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(30));
    a = 231; b = 189;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(21));
    a = 1; b = 17;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(1));
    a = 19; b = 1;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(1));
    a = 0; b = 17;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(17));
    a = 19; b = 0;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(19));
    a = 19; b = 17;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(1));
    a = 17; b = 19;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(1));
    a = 255; b = 255;
    EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(255));
    if constexpr (ut_numeric_limits<T>::digits >= 16) {
        a = 21945; b = 63525;
        EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(1155));
        a = 40755; b = 7623;
        EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(33));
    }
    if constexpr (ut_numeric_limits<T>::digits >= 32) {
        a = UINT32_C(2908157904); b = UINT32_C(1141161890);
        EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(65042));
    }
    if constexpr (ut_numeric_limits<T>::digits >= 64) {
        a = UINT64_C(434276666443008); b = UINT64_C(3846826911345880);
        EXPECT_TRUE(greatest_common_divisor(a, b) ==
                    static_cast<T>(UINT64_C(1677313784)));
        a = UINT64_C(278020828800); b = UINT64_C(513269738478);
        EXPECT_TRUE(greatest_common_divisor(a, b) == static_cast<T>(342));
    }
}

#ifdef TESTING_GREATEST_COMMON_DIVISOR_EUCLID
TEST(HurchallaFactoringGreatestCommonDivisor, gcd_euclid) {
#else
TEST(HurchallaFactoringGreatestCommonDivisor, gcd_binary) {
#endif
    test_gcd<uint8_t>();
    test_gcd<uint16_t>();
    test_gcd<uint32_t>();
    test_gcd<uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    test_gcd<__uint128_t>();
#endif
}


} // end namespace
