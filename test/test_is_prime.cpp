// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#include "is_prime_bruteforce.h"
#include "hurchalla/factoring/is_prime.h"

#include "gtest/gtest.h"
#include <cstdint>
#include <limits>
#include <cstddef>
#include <sstream>

namespace {


    TEST(HurchallaFactoringIsPrime, exhaustive_uint16_t) {
        namespace hc = hurchalla;
        using T = std::uint16_t;
        for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(hc::is_prime(x) == hc::is_prime_bruteforce(x));
        }
        T x = std::numeric_limits<T>::max();
        EXPECT_TRUE(hc::is_prime(x) == hc::is_prime_bruteforce(x));
    }

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 10/20/20.
    TEST(HurchallaFactoringIsPrime, exhaustive_uint32_t) {
        namespace hc = hurchalla;
        using T = std::uint32_t;
        for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
            EXPECT_TRUE(hc::is_prime(x) == hc::is_prime_bruteforce(x));
        }
        T x = std::numeric_limits<T>::max();
        EXPECT_TRUE(hc::is_prime(x) == hc::is_prime_bruteforce(x));
    }
#endif

/*
// I used this speed test to do a quick and dirty initial performance tuning of
// the value for the macro HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR, used in
// impl_is_prime.h.  This test is not needed normally.
    TEST(HurchallaFactoringIsPrime, speed_test32) {
        namespace hc = hurchalla;
        using T = std::uint32_t;
        T max = std::numeric_limits<T>::max()/2;
        T dummy = 0;
        for (T x = max; x >= max - 20000000; x = x-2) {
            bool b = hc::is_prime(x);
            // We need to prevent the compiler from completely removing
            // the is_prime calls due to b never being used.
            // So we'll add be to dummy just so it's used.
            dummy += b;
        }
        EXPECT_TRUE(dummy > 0);
    }
    TEST(HurchallaFactoringIsPrime, speed_test64) {
        namespace hc = hurchalla;
        using T = std::uint64_t;
        T max = std::numeric_limits<T>::max();
        T dummy = 0;
        for (T x = max; x >= max - 5000000; x = x-2) {
            bool b = hc::is_prime(x);
            // We need to prevent the compiler from completely removing
            // the is_prime calls due to b never being used.
            // So we'll add be to dummy just so it's used.
            dummy += b;
        }
        EXPECT_TRUE(dummy > 0);
    }
*/


    TEST(HurchallaFactoringIsPrime, basic_tests) {
        namespace hc = hurchalla;

        EXPECT_TRUE(hc::is_prime(UINT32_C(53)));
        EXPECT_TRUE(hc::is_prime(UINT32_C(127)));
        EXPECT_FALSE(hc::is_prime(UINT32_C(141)));
        EXPECT_FALSE(hc::is_prime(UINT32_C(140)));
        EXPECT_FALSE(hc::is_prime(UINT32_C(256)));
        EXPECT_TRUE(hc::is_prime(UINT32_C(67967)));
        EXPECT_TRUE(hc::is_prime(UINT32_C(67979)));
        EXPECT_FALSE(hc::is_prime(UINT32_C(67981)));
        EXPECT_FALSE(hc::is_prime(UINT32_C(67982)));

        EXPECT_TRUE(hc::is_prime(UINT64_C(53)));
        EXPECT_TRUE(hc::is_prime(UINT64_C(127)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(141)));
        EXPECT_TRUE(hc::is_prime(UINT64_C(67967)));
        EXPECT_TRUE(hc::is_prime(UINT64_C(67979)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(67981)));

        EXPECT_TRUE(hc::is_prime(UINT64_C(40000000003)));
        EXPECT_TRUE(hc::is_prime(UINT64_C(40000000031)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(40000000029)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(40000000027)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(40000000025)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(40000000024)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(8589934592)));  // 2^33

        EXPECT_TRUE(hc::is_prime(UINT64_C(18446744073709551557)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551558)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551555)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551554)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551553)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551551)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551549)));
        EXPECT_TRUE(hc::is_prime(UINT64_C(18446744073709551533)));
        EXPECT_TRUE(hc::is_prime(UINT64_C(18446744073709551521)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551523)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551525)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551527)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551529)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551530)));
        EXPECT_FALSE(hc::is_prime(UINT64_C(18446744073709551531)));
    }

#ifdef __SIZEOF_INT128__
    TEST(HurchallaFactoringIsPrime, basic_tests_128bit) {
        namespace hc = hurchalla;
        __uint128_t zero = 0;
        // rely on 2^128 wrap-around in the subtractions below
        EXPECT_TRUE(hc::is_prime(zero-159));
        EXPECT_TRUE(hc::is_prime(zero-173));
        EXPECT_TRUE(hc::is_prime(zero-233));
        EXPECT_FALSE(hc::is_prime(zero-160));
        EXPECT_FALSE(hc::is_prime(zero-161));
        EXPECT_FALSE(hc::is_prime(zero-163));
        EXPECT_FALSE(hc::is_prime(zero-165));
        EXPECT_FALSE(hc::is_prime(zero-167));

        EXPECT_TRUE(hc::is_prime(static_cast<__uint128_t>(40000000003)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(40000000001)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(40000000002)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(40000000005)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(40000000007)));

        EXPECT_TRUE(hc::is_prime(static_cast<__uint128_t>(53)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(54)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(55)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(49)));

        EXPECT_TRUE(hc::is_prime(static_cast<__uint128_t>(67967)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(67968)));
        EXPECT_FALSE(hc::is_prime(static_cast<__uint128_t>(67969)));
    }
#endif


} // end namespace
