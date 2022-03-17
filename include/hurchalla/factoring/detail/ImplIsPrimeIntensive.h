// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IMPL_IS_PRIME_INTENSIVE_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_IS_PRIME_INTENSIVE_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/detail/is_prime_trialdivision.h"
#include "hurchalla/factoring/detail/SieveOfEratosthenes.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace detail {


#ifndef HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE
// Some short perf tests on Haswell suggest 75 would be a decent value for
// HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE.  Tested with IsPrimeIntensive
// having OPTIMIZE_PRIMES==false and using uint64_t and __uint128_t.
// (FYI, size 54 would trial all prime factors < 256)
#  define HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE (75)
#endif

#ifndef HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_TYPE
#  define HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_TYPE PrimeTrialDivisionWarren
#endif


// Primary template declaration
template <typename T, bool OPTIMIZE_PRIMES = false, typename DUMMY=void>
struct ImplIsPrimeIntensive;


// helper parent class that utilizes Sieve of Eratosthenes.
template <typename T, typename DUMMY>
struct SieveImplIsPrimeIntensive {
private:
    static_assert(std::is_same<DUMMY, void>::value, "");
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits <= 32, "");
    const SieveOfEratosthenes sieve;
public:
    SieveImplIsPrimeIntensive() :
         sieve(static_cast<std::uint64_t>(1)<<(ut_numeric_limits<T>::digits)) {}
    bool operator()(T x) const
    {
        return sieve.isPrime(x);
    }
};
// 8 bit, uses Sieve of Eratosthenes (which is equally fast for testing of
// prime and composite numbers).
template <bool OPTIMIZE_PRIMES, typename DUMMY>
struct ImplIsPrimeIntensive<std::uint8_t, OPTIMIZE_PRIMES, DUMMY> :
                             SieveImplIsPrimeIntensive<std::uint8_t, DUMMY> {};
// 16 bit, uses Sieve of Eratosthenes (which is equally fast for testing of
// prime and composite numbers).
template <bool OPTIMIZE_PRIMES, typename DUMMY>
struct ImplIsPrimeIntensive<std::uint16_t, OPTIMIZE_PRIMES, DUMMY> :
                             SieveImplIsPrimeIntensive<std::uint16_t, DUMMY> {};
// 32 bit, uses Sieve of Eratosthenes (which is equally fast for testing of
// prime and composite numbers).
template <bool OPTIMIZE_PRIMES, typename DUMMY>
struct ImplIsPrimeIntensive<std::uint32_t, OPTIMIZE_PRIMES, DUMMY> :
                             SieveImplIsPrimeIntensive<std::uint32_t, DUMMY> {};


// 64bit, specialized version for when we expect the numbers we test will likely
// be prime
template <typename DUMMY>
struct ImplIsPrimeIntensive<std::uint64_t, true, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
    using T = std::uint64_t;
    bool operator()(T x) const
    {
        if (x % 2 == 0)
            return (x == 2);
        if (x < 2)
            return false;

// Visual Studio 2017 gets an internal compiler error (compiler bug) when
// using TOTAL_BASES = 2.  VS2017 works if we set TOTAL_BASES = 3.
#if !defined(_MSC_VER) || (_MSC_VER >= 1927)
        constexpr std::size_t TOTAL_BASES = 2;
#else
        // 64 bit primality testing with 3 bases uses a much smaller hash table
        // than testing with 2 bases, and thus uses less CPU cache - see
        // is_prime_miller_rabin.h for info.
        // For some CPUs the large table might not fit in cache (though I'd
        // still suspect 2 bases would be faster even so).
        constexpr std::size_t TOTAL_BASES = 3;
#endif
        // Since we expect x is usually prime (we're optimizing for primes here)
        // we use the largest possible trial size for miller-rabin testing to
        // maximize instruction level parallelism, knowing that we would very
        // likely need to go through all the bases even if we used the smallest
        // possible trial size of one single base - miller rabin always tests
        // all its bases whenever it gets a prime number.
        constexpr std::size_t TRIAL_SIZE = TOTAL_BASES;

#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  error "HURCHALLA_TARGET_BIT_WIDTH must be defined"
#endif
#if (HURCHALLA_TARGET_BIT_WIDTH <= 32)
        if (x <= ut_numeric_limits<std::uint32_t>::max()) {
            using U = std::uint32_t;
            if (x < (static_cast<U>(1) << 30)) {
                auto mf = MontgomeryQuarter<U>(static_cast<U>(x));
                return MillerRabinMontgomery<decltype(mf), 32, TRIAL_SIZE,
                                                     TOTAL_BASES>::is_prime(mf);
            }
            else {
                auto mf = MontgomeryForm<U>(static_cast<U>(x));
                return MillerRabinMontgomery<decltype(mf), 32, TRIAL_SIZE,
                                                     TOTAL_BASES>::is_prime(mf);
            }
        }
#endif
        if (x < (static_cast<T>(1) << 62)) {
            auto mf = MontgomeryQuarter<T>(x);
            return MillerRabinMontgomery<decltype(mf), 64, TRIAL_SIZE,
                                                 TOTAL_BASES>::is_prime(mf);
        }
        else {
            auto mf = MontgomeryForm<T>(x);
            return MillerRabinMontgomery<decltype(mf), 64, TRIAL_SIZE,
                                                 TOTAL_BASES>::is_prime(mf);
        }
    }
};


// 64bit, specialized version for when we don't particularly expect the numbers
// we test will be prime
template <typename DUMMY>
struct ImplIsPrimeIntensive<std::uint64_t, false, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
    using T = std::uint64_t;
    bool operator()(T x) const
    {
        bool success;
        // Using trial division on average boosts our performance (so long as x
        // is not especially likely to be prime), because it avoids miller-rabin
        // for composites that have a small enough factor.
        bool isPrime = is_prime_trialdivision::call<
                         HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_TYPE,
                         HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE>(x, success);
        if (success)
            return isPrime;
        // is_prime_trialdivision should have successfully handled any even x,
        // and any x < 2.
        HPBC_ASSERT2(x % 2 != 0);
        HPBC_ASSERT2(x >= 2);

// Visual Studio 2017 gets an internal compiler error (compiler bug) when
// using TOTAL_BASES = 2.  VS2017 works if we set TOTAL_BASES = 3.
#if !defined(_MSC_VER) || (_MSC_VER >= 1927)
        constexpr std::size_t TOTAL_BASES = 2;
        constexpr std::size_t TRIAL_SIZE = 1;
#else
        constexpr std::size_t TOTAL_BASES = 3;
        constexpr std::size_t TRIAL_SIZE = 2;
#endif

#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  error "HURCHALLA_TARGET_BIT_WIDTH must be defined"
#endif
#if (HURCHALLA_TARGET_BIT_WIDTH <= 32)
        if (x <= ut_numeric_limits<std::uint32_t>::max()) {
            using U = std::uint32_t;
            if (x < (static_cast<U>(1) << 30)) {
                auto mf = MontgomeryQuarter<U>(static_cast<U>(x));
                return MillerRabinMontgomery<decltype(mf), 32, TRIAL_SIZE,
                                                     TOTAL_BASES>::is_prime(mf);
            }
            else {
                auto mf = MontgomeryForm<U>(static_cast<U>(x));
                return MillerRabinMontgomery<decltype(mf), 32, TRIAL_SIZE,
                                                     TOTAL_BASES>::is_prime(mf);
            }
        }
#endif
        if (x < (static_cast<T>(1) << 62)) {
            auto mf = MontgomeryQuarter<T>(x);
            return MillerRabinMontgomery<decltype(mf), 64, TRIAL_SIZE,
                                                 TOTAL_BASES>::is_prime(mf);
        }
        else {
            auto mf = MontgomeryForm<T>(x);
            return MillerRabinMontgomery<decltype(mf), 64, TRIAL_SIZE,
                                                 TOTAL_BASES>::is_prime(mf);
        }
    }
};


// primary template implementation.  Requires T is a 128 bit type
template <typename T, bool OPTIMIZE_PRIMES, typename DUMMY>
struct ImplIsPrimeIntensive {
private:
    // primary template handles any 128 bit unsigned integer type
    static_assert(std::is_same<DUMMY, void>::value, "");
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits == 128, "");
#if (HURCHALLA_TARGET_BIT_WIDTH <= 64)
    const ImplIsPrimeIntensive<std::uint64_t, false, DUMMY> ipi64;
public:
    ImplIsPrimeIntensive() : ipi64() {}
#endif
public:
    bool operator()(T x) const
    {
#if (HURCHALLA_TARGET_BIT_WIDTH <= 64)
        if (x <= ut_numeric_limits<uint64_t>::max())
            return ipi64(static_cast<std::uint64_t>(x));
#endif
        // For now, we won't do anything special for when we expect x is prime
        // or not prime (we ignore OPTIMIZE_PRIMES).  Currently the primality
        // testing of 128 bit numbers is not well optimized because it depends
        // on a total of 128 bases with probabilistic miller-rabin testing.
        // Such a large number of bases is inherently slow for primes, so the
        // fact that this function does some preliminary trial division that's
        // pointless for a prime probably won't matter much in comparison.

        // Using trial division on average boosts our performance (so long as x
        // is not especially likely to be prime), because it avoids miller-rabin
        // for composites that have a small enough factor.
        bool success;
        bool isPrime = is_prime_trialdivision::call<
                         HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_TYPE,
                         HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE>(x, success);
        if (success)
            return isPrime;
        // is_prime_trialdivision should have successfully handled any even x,
        // and any x < 2.
        HPBC_ASSERT2(x % 2 != 0);
        HPBC_ASSERT2(x >= 2);
        return is_prime_miller_rabin::call(x);
    }
};


}}  // end namespace

#endif
