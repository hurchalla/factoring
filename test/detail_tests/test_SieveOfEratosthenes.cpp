// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "hurchalla/factoring/detail/SieveOfEratosthenes.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include <cstdint>

#include "gtest/gtest.h"

namespace {


using namespace hurchalla::detail;


void exhaustive_sieve_test(const SieveOfEratosthenes& sieve)
{
    EXPECT_TRUE(sieve.size() > 0);
    EXPECT_FALSE(sieve[0]);
    EXPECT_FALSE(sieve[1]);
    EXPECT_TRUE(sieve[2]);
    for (uint64_t i=3; i<sieve.size()-1; i+=2) {
        EXPECT_TRUE(sieve[i] == is_prime_miller_rabin_integral(i));
        EXPECT_FALSE(sieve[i+1]);
    }
    uint64_t i = sieve.size() - 1;
    if (i % 2 == 0)
        EXPECT_FALSE(sieve[i]);
    else
        EXPECT_TRUE(sieve[i] == is_prime_miller_rabin_integral(i));
}

TEST(HurchallaFactoringSieve, sieve_uint8_t) {
    SieveOfEratosthenes sieve(static_cast<std::uint64_t>(1) << 8);
    exhaustive_sieve_test(sieve);
}
TEST(HurchallaFactoringSieve, sieve_uint16_t) {
    SieveOfEratosthenes sieve(static_cast<std::uint64_t>(1) << 16);
    exhaustive_sieve_test(sieve);
}
#if 0
// Ordinarily you don't want to run this since it takes ~5 minutes to complete.
TEST(HurchallaFactoringSieve, sieve_uint32_t) {
    SieveOfEratosthenes sieve(static_cast<std::uint64_t>(1) << 32);
    exhaustive_sieve_test(sieve);
}
#endif


} // end namespace
