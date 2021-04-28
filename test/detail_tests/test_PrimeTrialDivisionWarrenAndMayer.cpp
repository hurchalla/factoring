// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include <cstdint>
#include <limits>

#include "gtest/gtest.h"

namespace {


using namespace hurchalla;
using namespace hurchalla::detail;

//PTD can be PrimeTrialDivisionMayer<T,SIZE> or PrimeTrialDivisionWarren<T,SIZE>
template <typename PTD>
void isDivisible_test(typename PTD::value_type x, int index)
{
    using T = typename PTD::value_type;
    constexpr int SIZE = PTD::size();
    // the caller of this test is faulty if this assertion is false
    EXPECT_TRUE(0 <= index && index < SIZE);

    T div_result;
    bool is_divisible = PTD::isDivisible(div_result, x, index);
    T divisor = PTD::oddPrime(index);
    EXPECT_TRUE(is_divisible == (x % divisor == 0));
    if (is_divisible)
        EXPECT_TRUE(div_result == x/divisor);

    // We don't actually care about squaring the divisor -
    // it's just a convenient way to test oddPrimeSquared
    auto divisor_squared = PTD::oddPrimeSquared(index);
    // We certainly expect all primes in PTD are under 4 billion (more precisely
    // under UINT32_MAX), considering that they are created at compile-time in a
    // static array
    EXPECT_TRUE(divisor < (static_cast<std::uint64_t>(1) << 32));
    // Thus divisor squared should fit in a uint64_t, and we can test that
    // oddPrimeSquared() gave us the correct result
    std::uint64_t divisor64 = static_cast<std::uint64_t>(divisor);
    EXPECT_TRUE(divisor_squared == divisor64 * divisor64);
}


//PTD can be PrimeTrialDivisionMayer<T,SIZE> or PrimeTrialDivisionWarren<T,SIZE>
template <typename PTD>
void PTD_tests()
{
    using T = typename PTD::value_type;
    constexpr int SIZE = PTD::size();
    static_assert(SIZE > 0);
    static_assert(PTD::oddPrime(0) == 3);

    // get the first prime larger than the last prime used by PTD
    constexpr auto next_prime = PTD::nextPrimePastEnd();
    // get the square of PTD::nextPrimePastEnd() without overflow.
    constexpr auto next_prime_squared = PTD::nextPrimePastEndSquared();

    // We don't expect next_prime to have a type larger than uint32_t,
    // and thus next_prime squared should fit in a uint64_t.
    static_assert(std::numeric_limits<decltype(next_prime)>::is_integer, "");
    static_assert(!std::numeric_limits<decltype(next_prime)>::is_signed, "");
    static_assert(std::numeric_limits<decltype(next_prime)>::digits <= 32, "");
    static_assert(next_prime_squared == static_cast<std::uint64_t>(next_prime) *
                                        static_cast<std::uint64_t>(next_prime));

    static_assert(!std::numeric_limits<T>::is_signed);
    T midpoint = (static_cast<T>(0) - 1)/2;
    if (midpoint % 2 == 0)
        ++midpoint;

    T max = static_cast<T>(static_cast<T>(0) - 1);
    if (max % 2 == 0)
        --max;

    constexpr T range = 30;
    constexpr int indexrange = (SIZE > 40) ? 10 : SIZE/4;

    for (T x = 0; x < range; ++x) {
        for (int i = 0; i < indexrange; ++i)
            isDivisible_test<PTD>(x, i);
        for (int i = SIZE-1; i > SIZE - indexrange; --i)
            isDivisible_test<PTD>(x, i);
        for (int i = SIZE/2 - indexrange/2; i < SIZE/2 + indexrange/2; ++i)
            isDivisible_test<PTD>(x, i);
    }
    for (T x = max; x > max - range; --x) {
        for (int i = 0; i < indexrange; ++i)
            isDivisible_test<PTD>(x, i);
        for (int i = SIZE-1; i > SIZE - indexrange; --i)
            isDivisible_test<PTD>(x, i);
        for (int i = SIZE/2 - indexrange/2; i < SIZE/2 + indexrange/2; ++i)
            isDivisible_test<PTD>(x, i);
    }
    for (T x = midpoint - range/2; x < midpoint + range/2; ++x) {
        for (int i = 0; i < indexrange; ++i)
            isDivisible_test<PTD>(x, i);
        for (int i = SIZE-1; i > SIZE - indexrange; --i)
            isDivisible_test<PTD>(x, i);
        for (int i = SIZE/2 - indexrange/2; i < SIZE/2 + indexrange/2; ++i)
            isDivisible_test<PTD>(x, i);
    }
}


template <template <class, int> class TPTD>
void run_PTD_tests_uint8()
{
    // there are 53 odd primes that are less than 256
    using PTD8_53 = TPTD<std::uint8_t, 53>;
    static_assert(PTD8_53::nextPrimePastEnd() == 257);
    PTD_tests<PTD8_53>();

#if 0
    // If we use 54 or above we should get a compilation error.  I enabled this
    // section once to test it, and as desired, it failed to compile.
    using PTD8_54 = TPTD<std::uint8_t, 54>;
    PTD8_54::oddPrimes(0);  // this should cause instantiation, and thus failure
#endif

    using PTD8_52 = TPTD<std::uint8_t, 52>;
    static_assert(PTD8_52::nextPrimePastEnd() == 251);
    PTD_tests<PTD8_52>();

    using PTD8_1 = TPTD<std::uint8_t, 1>;
    static_assert(PTD8_1::nextPrimePastEnd() == 5);
    PTD_tests<PTD8_1>();
}

template <template <class, int> class TPTD, typename T>
void run_PTD_tests_larger()
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits >= 16, "");

    // there are 53 primes that are less than 256
    using PTD_53 = TPTD<T, 53>;
    static_assert(PTD_53::nextPrimePastEnd() == 257);
    PTD_tests<PTD_53>();

    // 54 and above should work fine for uint16_t and larger types
    using PTD_54 = TPTD<T, 54>;
    static_assert(PTD_54::nextPrimePastEnd() == 263);
    PTD_tests<PTD_54>();

    using PTD_52 = TPTD<T, 52>;
    static_assert(PTD_52::nextPrimePastEnd() == 251);
    PTD_tests<PTD_52>();

    using PTD_1 = TPTD<T, 1>;
    static_assert(PTD_1::nextPrimePastEnd() == 5);
    PTD_tests<PTD_1>();

    using PTD_2000 = TPTD<T, 800>;
    static_assert(PTD_2000::nextPrimePastEnd() == 6151);
    PTD_tests<PTD_2000>();
}


TEST(HurchallaFactoringPrimeTrialDivisionMayer, PTDMayer_tests) {
    run_PTD_tests_uint8<PrimeTrialDivisionMayer>();
    run_PTD_tests_larger<PrimeTrialDivisionMayer, std::uint16_t>();
    run_PTD_tests_larger<PrimeTrialDivisionMayer, std::uint32_t>();
    run_PTD_tests_larger<PrimeTrialDivisionMayer, std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    run_PTD_tests_larger<PrimeTrialDivisionMayer, __uint128_t>();
#endif
}

TEST(HurchallaFactoringPrimeTrialDivisionWarren, PTDWarren_tests) {
    run_PTD_tests_uint8<PrimeTrialDivisionWarren>();
    run_PTD_tests_larger<PrimeTrialDivisionWarren, std::uint16_t>();
    run_PTD_tests_larger<PrimeTrialDivisionWarren, std::uint32_t>();
    run_PTD_tests_larger<PrimeTrialDivisionWarren, std::uint64_t>();
#if HURCHALLA_COMPILER_HAS_UINT128_T()
    run_PTD_tests_larger<PrimeTrialDivisionWarren, __uint128_t>();
#endif
}


} // end namespace
