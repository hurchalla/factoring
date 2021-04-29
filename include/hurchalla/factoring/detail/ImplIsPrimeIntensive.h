// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IMPL_IS_PRIME_INTENSIVE_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_IS_PRIME_INTENSIVE_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/detail/is_prime_trialdivision.h"
#include "hurchalla/factoring/detail/SieveOfEratosthenes.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace detail {


// TODO: empirically find good default(s) for
// HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE

#ifndef HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE
// Some short perf testing on Haswell suggest 100 would be a good value for
// intensive use of PrimeTrialDivisionWarren with uint64_t.  Presumably this
// would be more or less ok with __uint128_t but I haven't tried it.
// FYI, size 54 would trial all prime factors < 256
#  define HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE (100)
#endif


// primary template
template <typename T, bool OPTIMIZE_PRIMES = false, typename DUMMY=void>
struct ImplIsPrimeIntensive {
    // primary template handles any 128 bit unsigned integer type
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits == 128, "");
    static_assert(std::is_same<DUMMY, void>::value, "");
    T operator()(T x) const
    {
        // For now, we won't do anything special for when we expect x is prime
        // or not prime (we ignore OPTIMIZE_PRIMES).  Currently the primality
        // testing of 128 bit numbers is not well optimized because it depends
        // on a total of 128 bases with probabilistic miller-rabin testing.
        // Such a large number of bases is inherently slow for primes, so the
        // fact that this function does some preliminary trial division that's
        // pointless for a prime won't matter much in comparison.

        // try all prime factors < 256
        bool success;
        bool res = is_prime_trialdivision<PrimeTrialDivisionWarren,
                         HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE>(x, success);
        if (success)
            return res;
        // is_prime_trialdivision() should have successfully handled any even x,
        // and any x < 2.
        HPBC_ASSERT2(x % 2 != 0);
        HPBC_ASSERT2(x >= 2);

        return detail::is_prime_miller_rabin_integral(x);
    }
};


// helper parent class that utilizes Sieve of Eratosthenes.
template <typename T, typename DUMMY>
struct SieveImplIsPrimeIntensive {
private:
    static_assert(std::is_same<DUMMY, void>::value, "");
    static_assert(ut_numeric_limits<T>::digits < 64, "");
    const SieveOfEratosthenes sieve;
public:
    SieveImplIsPrimeIntensive() :
         sieve(static_cast<std::uint64_t>(1)<<(ut_numeric_limits<T>::digits)) {}
    T operator()(T x) const
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
    T operator()(T x) const
    {
        if (x % 2 == 0)
            return (x == 2);
        if (x < 2)
            return false;
#if 1
        constexpr std::size_t TOTAL_BASES = 2;
#else
        // 64 bit primality testing with 3 bases uses a much smaller hash table
        // than 2 testing with 2 bases, and thus uses less CPU cache - see
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
    T operator()(T x) const
    {
        bool success;
        // try all prime factors < 256
        bool res = is_prime_trialdivision<PrimeTrialDivisionWarren,
                         HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE>(x, success);
        if (success)
            return res;
        // is_prime_trialdivision() should have successfully handled any even x,
        // and any x < 2.
        HPBC_ASSERT2(x % 2 != 0);
        HPBC_ASSERT2(x >= 2);
        return detail::is_prime_miller_rabin_integral(x);
    }
};


}}  // end namespace

#endif
