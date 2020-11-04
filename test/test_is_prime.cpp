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
        namespace hf = hurchalla::factoring;
        using T = std::uint16_t;
        for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(hf::is_prime(x) == hf::is_prime_bruteforce(x));
        }
        T x = std::numeric_limits<T>::max();
        EXPECT_TRUE(hf::is_prime(x) == hf::is_prime_bruteforce(x));
    }

#if 0
// Ordinarily you don't want to run this since it takes ~2.5 hours to complete.
// It passed when I tested it on 10/20/20.
    TEST(HurchallaFactoringIsPrime, exhaustive_uint32_t) {
        namespace hf = hurchalla::factoring;
        using T = std::uint32_t;
        for (T x = 0; x < std::numeric_limits<T>::max(); ++x) {
            EXPECT_TRUE(hf::is_prime(x) == hf::is_prime_bruteforce(x));
        }
        T x = std::numeric_limits<T>::max();
        EXPECT_TRUE(hf::is_prime(x) == hf::is_prime_bruteforce(x));
    }
#endif

    TEST(HurchallaFactoringIsPrime, basic_tests) {
        namespace hf = hurchalla::factoring;

        EXPECT_TRUE(hf::is_prime(UINT32_C(53)));
        EXPECT_TRUE(hf::is_prime(UINT32_C(127)));
        EXPECT_FALSE(hf::is_prime(UINT32_C(141)));
        EXPECT_FALSE(hf::is_prime(UINT32_C(140)));
        EXPECT_FALSE(hf::is_prime(UINT32_C(256)));
        EXPECT_TRUE(hf::is_prime(UINT32_C(67967)));
        EXPECT_TRUE(hf::is_prime(UINT32_C(67979)));
        EXPECT_FALSE(hf::is_prime(UINT32_C(67981)));
        EXPECT_FALSE(hf::is_prime(UINT32_C(67982)));

        EXPECT_TRUE(hf::is_prime(UINT64_C(53)));
        EXPECT_TRUE(hf::is_prime(UINT64_C(127)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(141)));
        EXPECT_TRUE(hf::is_prime(UINT64_C(67967)));
        EXPECT_TRUE(hf::is_prime(UINT64_C(67979)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(67981)));

        EXPECT_TRUE(hf::is_prime(UINT64_C(40000000003)));
        EXPECT_TRUE(hf::is_prime(UINT64_C(40000000031)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(40000000029)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(40000000027)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(40000000025)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(40000000024)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(8589934592)));  // 2^33

        EXPECT_TRUE(hf::is_prime(UINT64_C(18446744073709551557)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551558)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551555)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551554)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551553)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551551)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551549)));
        EXPECT_TRUE(hf::is_prime(UINT64_C(18446744073709551533)));
        EXPECT_TRUE(hf::is_prime(UINT64_C(18446744073709551521)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551523)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551525)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551527)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551529)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551530)));
        EXPECT_FALSE(hf::is_prime(UINT64_C(18446744073709551531)));
    }

#ifdef __SIZEOF_INT128__
    TEST(HurchallaFactoringIsPrime, basic_tests_128bit) {
        namespace hf = hurchalla::factoring;
        __uint128_t zero = 0;
        // rely on 2^128 wrap-around in the subtractions below
        EXPECT_TRUE(hf::is_prime(zero-159));
        EXPECT_TRUE(hf::is_prime(zero-173));
        EXPECT_TRUE(hf::is_prime(zero-233));
        EXPECT_FALSE(hf::is_prime(zero-160));
        EXPECT_FALSE(hf::is_prime(zero-161));
        EXPECT_FALSE(hf::is_prime(zero-163));
        EXPECT_FALSE(hf::is_prime(zero-165));
        EXPECT_FALSE(hf::is_prime(zero-167));

        EXPECT_TRUE(hf::is_prime(static_cast<__uint128_t>(40000000003)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(40000000001)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(40000000002)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(40000000005)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(40000000007)));

        EXPECT_TRUE(hf::is_prime(static_cast<__uint128_t>(53)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(54)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(55)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(49)));

        EXPECT_TRUE(hf::is_prime(static_cast<__uint128_t>(67967)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(67968)));
        EXPECT_FALSE(hf::is_prime(static_cast<__uint128_t>(67969)));
    }
#endif


} // end namespace
