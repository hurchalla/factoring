// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/SieveOfEratosthenes.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/detail/is_prime_trialdivision.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <vector>

#include "gtest/gtest.h"

namespace {


using namespace hurchalla;
using namespace hurchalla::detail;


template <typename T>
bool get_primality(T x)
{
    static_assert(!ut_numeric_limits<T>::is_signed);
    if (x < 2)
        return false;
    if (x % 2 == 0)
        return x == 2;
    return is_prime_miller_rabin_integral(x);
}


template <int SIZE, typename T>
void iptd_test(T x, const std::vector<std::uint64_t>& primevec)
{
    static_assert(!ut_numeric_limits<T>::is_signed);
    using size_type = std::vector<std::uint64_t>::size_type;
    // a failure here would mean this test is buggy
    static_assert(ut_numeric_limits<int>::digits
                  <= ut_numeric_limits<size_type>::digits);
    // a failure here would mean this test is buggy
    EXPECT_TRUE(primevec.size() == SIZE + 1);
    std::uint64_t nextprime = primevec[static_cast<size_type>(SIZE)];
    // These tests need for nextprime*nextprime to still fit in a uint64_t.
    // The tests would be at fault if this is false, but mostly just for
    // choosing a SIZE that's way larger than needed.
    EXPECT_TRUE(nextprime < static_cast<std::uint64_t>(1)<<32);

    bool success, isprime;
    isprime = is_prime_trialdivision<PrimeTrialDivisionMayer, SIZE>(x, success);
    if (x < nextprime*nextprime) {
        EXPECT_TRUE(success);
    }
    if (success) {
        EXPECT_TRUE(isprime == get_primality(x));
    }

    isprime = is_prime_trialdivision<PrimeTrialDivisionWarren,SIZE>(x, success);
    if (x < nextprime*nextprime) {
        EXPECT_TRUE(success);
    }
    if (success) {
        EXPECT_TRUE(isprime == get_primality(x));
    }
}


template <int SIZE, typename T>
void iptd_sized_tests(const SieveOfEratosthenes& sieve)
{
    static_assert(1 < SIZE);
    static_assert(!ut_numeric_limits<T>::is_signed);

    std::vector<std::uint64_t> primevec;
    primevec.push_back(2);
    for (std::uint64_t i=3; ; i+=2) {
        EXPECT_TRUE(i < sieve.size());  // false would mean this test is buggy
        if (sieve[i]) {
            primevec.push_back(i);
            if (primevec.size() > SIZE)
                break;
        }
    }
    // a failure here would mean this test is buggy
    EXPECT_TRUE(primevec.size() == SIZE + 1);

    using size_type = std::vector<std::uint64_t>::size_type;
    // a failure here would mean this test is buggy
    static_assert(ut_numeric_limits<int>::digits
                  <= ut_numeric_limits<size_type>::digits);

    T max = static_cast<T>(static_cast<T>(0) - 1);
    T midpoint = max/2;
    T midpoint_minus50 = static_cast<T>(midpoint - 50);
    for (T x = 0; x < 255; ++x)
        iptd_test<SIZE>(x, primevec);
    for (T x = max; x > max - 100; --x)
        iptd_test<SIZE>(x, primevec);
    for (T x = midpoint_minus50; x < midpoint + 50; ++x)
        iptd_test<SIZE>(x, primevec);

    std::vector<int> indices = { 0, 1, 2, SIZE, SIZE-1, SIZE/2, SIZE/2 + 1 };
    if (5 < primevec.size())
        indices.push_back(5);
    for (auto index : indices) {
        // a failure here would mean this test is buggy
        EXPECT_TRUE(static_cast<unsigned int>(index) < primevec.size());
        std::uint64_t prime = primevec[static_cast<size_type>(index)];
        // it's fine if we under/overflow below- it will wrap around harmlessly.
        iptd_test<SIZE>(static_cast<T>(prime-2), primevec);
        iptd_test<SIZE>(static_cast<T>(prime-1), primevec);
        iptd_test<SIZE>(static_cast<T>(prime+0), primevec);
        iptd_test<SIZE>(static_cast<T>(prime+1), primevec);
        iptd_test<SIZE>(static_cast<T>(prime+2), primevec);
    }
}


template <typename T>
void iptd_typed_tests(const SieveOfEratosthenes& sieve)
{
    iptd_sized_tests<2,T>(sieve);
    iptd_sized_tests<10,T>(sieve);
    iptd_sized_tests<54,T>(sieve);
    if constexpr (ut_numeric_limits<T>::digits >= 16) {
        iptd_sized_tests<55,T>(sieve);
        iptd_sized_tests<198,T>(sieve);
        iptd_sized_tests<1000,T>(sieve);
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
        iptd_sized_tests<2501,T>(sieve);
    }
#endif
}


TEST(HurchallaFactoringIsPrimeTrialDivision, iptd_tests) {
    SieveOfEratosthenes sieve(static_cast<std::uint64_t>(1)<<24);

    iptd_typed_tests<std::uint8_t>(sieve);
    iptd_typed_tests<std::uint16_t>(sieve);
    iptd_typed_tests<std::uint32_t>(sieve);
    iptd_typed_tests<std::uint64_t>(sieve);
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    iptd_typed_tests<__uint128_t>(sieve);
#endif
}


} // end namespace
