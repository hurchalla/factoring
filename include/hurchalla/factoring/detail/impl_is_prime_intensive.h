// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FUNCTION_IMPL_IS_PRIME_INTENSIVE_H_INCLUDED
#define HURCHALLA_FACTORING_FUNCTION_IMPL_IS_PRIME_INTENSIVE_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_trialdivision.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstddef>
#include <type_traits>

namespace hurchalla { namespace detail {


// Note: we use a struct and static functions in order to disallow ADL
struct impl_is_prime_intensive {

  template <typename MontType>
  static bool mont_miller_rabin(const MontType& mf)
  {
    using T = typename MontType::IntegerType;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    using U = typename extensible_make_unsigned<T>::type;

    constexpr int digitsT = ut_numeric_limits<T>::digits;
    // MillerRabinBases produces bases that are (usually) type uint16_t, and
    // is_miller_rabin.h needs the bases to fit in type T.  Thus, T with < 16
    // digits won't work.
    static_assert(digitsT >= 16);

    T x = mf.getModulus();
    HPBC_PRECONDITION2(x > 1);
    U xu = static_cast<U>(x);

    if constexpr(digitsT <= 32) {
        constexpr std::size_t TOTAL_BASES = 1;
        constexpr std::size_t TRIAL_SIZE = 1;
        const auto bases = MillerRabinBases<digitsT, TOTAL_BASES>::get(xu);
        return IPMR_internal::miller_rabin_trials<TRIAL_SIZE>(mf, bases);
    }
    else if constexpr(32 < digitsT && digitsT <= 64) {
        constexpr std::size_t TOTAL_BASES = 2;
        // Specifying TRIAL_SIZE 2 causes both tests to run at the same time.
        constexpr std::size_t TRIAL_SIZE = 2;
        const auto bases = MillerRabinBases<digitsT, TOTAL_BASES>::get(xu);
        return IPMR_internal::miller_rabin_trials<TRIAL_SIZE>(mf, bases);
        //
        // Note: we *always* run *two* miller-rabin tests above.  Normally an
        // implementation would do a single MR test against the first base
        // (which can almost always detect a composite), and then it would do
        // the remainder of the tests if still needed.  [FYI with hashing, we'll
        // only have a max of two tests here.]  We choose to *always* run both
        // at the same time to take advantage of instruction level parallelism
        // of the CPU (pipelined and super scalar execution), which would not
        // be otherwise utilized due to the long dependency chain present in a
        // Miller-Rabin test and in the montgomery multiplications.  We can run
        // 2 tests at the same time at an only slightly greater cost than a
        // single test.  Since we have already done a trial division step before
        // calling this function, an arbitrary 64 bit number that we get here
        // (this 'if' clause handles 64 bit numbers) would (with the default
        // TRIAL_DIVISION_SIZE) have very roughly about a 0.25 chance of being
        // prime and 0.75 chance of being composite.
        // So, let's assume we hit this clause four times.  We'd expect we'd get
        // on average 3 composites and 1 prime.  The normal implementation
        // aproach would need 1 Miller-Rabin test for each composite, and 2 M-R
        // tests for the prime, which is a total cost of 5 MR tests.  If we say
        // for the sake of argument that running 2 tests at the same time costs
        // only about 1.25 times as much as running one normal MR test, then if
        // we *always* run two MR tests at the same time, testing each composite
        // would cost about 1.25 (normal) MR tests, and then the prime would
        // also cost 1.25 (normal) MR tests.  This would be a total cost of 5
        // (normal) MR tests, which is the same as we estimated for the normal
        // implementation.  But in practice, the performance I've measured for 
        // two simultaneous MR tests has cost less than 1.25 (normal) MR tests.
        // Furthermore, when we test for example 48 bit numbers rather than full
        // 64 bit numbers, the odds of the number(after surviving trial division
        // and reaching here) being prime are greater than 0.25 because it is
        // only a 48 bit number. [The odds go up as the number becomes smaller.]
        // This hurts the performance of the "normal" implementation approach
        // because using the "normal" approach, the cost of evaluating a prime
        // is significantly more than the cost of evaluating a composite.  In
        // contrast, the approach we use here does not lose any performance as
        // the proportion of primes increases, because it always runs exactly
        // two tests regardless of whether a number is prime or composite.
        // Finally, with the traditional approach, the worst case performance
        // for determining primality of any given number, is much worse than its
        // average performance.  For the approach we use above, the worst case
        // performance and the average case performance are the same.
        //
        // The explanation above is theory...  In the actual performance
        // measurements I've done, the approach above was faster than the
        // traditional approach.  This approach (it uses two bases from hashing)
        // also measured faster than my experiments with Baillie PSW (BPSW).
        // https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test
        // BPSW provides a deterministic primality test for up to 64 bit numbers
        // just like the approach here.  Yet if BPSW measures no faster, then
        // there is not much reason to use BPSW for 64 bit primality testing,
        // unless we are running in an environment in which it is harmful or not
        // possible to store in memory (or access) the 448 KB hash table that
        // is needed by our 2 base Miller-Rabin hashing.  Also as a practical
        // matter, BPSW is comprised of an MR test and a Lucas test, and the
        // Lucas portion would add extra code over the MR testing we have
        // implemented here, which would increase the potential of coding error.
    }
    else if constexpr(64 < digitsT && digitsT <= 128) {
        // For the rationale behind the following static_assert, see the struct
        // MillerRabinMontgomery<MontType, 128, TRIAL_SIZE, 128>
        // inside is_prime_miller_rabin.h
        static_assert(ut_numeric_limits<T>::digits == 128 ||
                           (ut_numeric_limits<T>::is_signed &&
                            ut_numeric_limits<T>::digits == 127));

        // Use a 13 base test if x is small enough.  It's a lot faster than the
        // 128 base test below.
        // Note: 3317044064679887385961981 == (179817<<64) + 5885577656943027709
        constexpr U limit13 =
                 (static_cast<U>(179817) << 64) + UINT64_C(5885577656943027709);
        if (xu < limit13) {
            constexpr std::size_t TRIAL_SIZE = 3;
            return is_prime_miller_rabin_special::
                          case_3317044064679887385961981_128_13<TRIAL_SIZE>(mf);
        }

        // 128 bit miller-rabin with 128 bases is going to be slow no matter
        // what, but a trial size of 3 will usually improve performance over
        // trial size 1, due to more efficient use of the CPU's pipelined and/or
        // superscalar execution units.
        // Strictly from a performance standpoint we could probably do (much?)
        // better by adding one or more Lucas tests similar to BPSW
        // https://en.wikipedia.org/wiki/Baillie%E2%80%93PSW_primality_test
        // perhaps thus greatly reducing the need for large numbers of
        // Miller-Rabin tests, like the admittedly overly large 128 tests here.
        // I've only used BPSW (i.e. Lucas) for experimental purposes so far-
        // hurchalla/factoring/detail/experimental/bpsw/
        // It would take too much time for me to fully understand Lucas tests,
        // and BPSW, and write my own fully trustworthy code for the Lucas test.
        // Correctness is by far the most important goal for this project.
        // [This project doesn't include algorithms or code that I don't fully
        // understand, with the sole exception of ECM and that's only because
        // the ecm code trivially checks its results for correctness and falls
        // back to a much slower method (pollard-rho) if unexpectedly it's in
        // error.  There's no way to trivially check any primality test, which
        // includes BPSW.]
        constexpr std::size_t TRIAL_SIZE = 3;
        const auto& bases = MillerRabinProbabilisticBases128<>::bases;
#if 0
        return IPMR_internal::miller_rabin_trials<TRIAL_SIZE>(mf, bases);
#else
        // using miller_rabin_trials128 improves perf by tuning to the modulus size
        return IPMR_internal::miller_rabin_trials128<TRIAL_SIZE>(mf, bases);
#endif
    }
    else {
        // C++ treats static_assert in constexpr-if as ill-formed if it is
        // always false and does not depend on a template param.  So we use
        // sizeof(T)==0 here instead of plain 'false'.
        static_assert(sizeof(T) == 0, "digitsT <= 128 required");
    }
  }



  template <typename T>
  constexpr static unsigned int defaultTrialDivisionSize()
  {
    return (ut_numeric_limits<T>::digits > 32) ? 150u : 80u;
  }


  template <unsigned int TRIAL_DIVISION_SIZE, typename T>
  static bool call(T x)
  {
    static_assert(ut_numeric_limits<T>::is_integer);
    HPBC_PRECONDITION2(x >= 0);

    using U = typename extensible_make_unsigned<T>::type;
    constexpr int digitsU = ut_numeric_limits<U>::digits;

    if constexpr (digitsU > 32) {
        static_assert(digitsU % 2 == 0);
        // if possible, recurse to a call that will complete faster.
        if (x < (static_cast<T>(1) << (digitsU/2))) {
            using uint_halfbits_U = typename sized_uint<digitsU/2>::type;
            return call<TRIAL_DIVISION_SIZE>(static_cast<uint_halfbits_U>(x));
        }
    }
    constexpr int digitsT = ut_numeric_limits<T>::digits;
    static_assert(digitsT <= 128, "We disallow > 128 bit types because we have"
           "no primality checking algorithm implemented for numbers > 128bit.");

    // First try small trial divisions to find easy factors.
    // If primality still unknown, use miller-rabin to prove prime or not prime.

    if constexpr (digitsU <= 8) {
        // note: the 6th prime is 13, which is the last prime under 2^4 - thus
        // trial division using all primes up to and including the 6th prime
        // will be sufficient to determine primality of any number under 2^8.
        constexpr int SIZE = 6;
        bool success;
        bool isprime = is_prime_trialdivision::call<PrimeTrialDivisionWarren,
                                              SIZE>(static_cast<U>(x), success);
        // trial division with the first 6 primes should always succeed here
        HPBC_ASSERT2(success);
        return isprime;
    } else if constexpr (digitsU <= 16) {
        // note: the 54th prime is 251, which is the last prime under 2^8 - thus
        // trial division using all primes up to and including the 54th prime
        // will be sufficient to determine primality of any number under 2^16.
        constexpr int SIZE = 54;
        bool success;
        bool isprime = is_prime_trialdivision::call<PrimeTrialDivisionWarren,
                                              SIZE>(static_cast<U>(x), success);
        // trial division with the first 54 primes should always succeed here
        HPBC_ASSERT2(success);
        return isprime;
    } else {
        if constexpr (TRIAL_DIVISION_SIZE > 0)
        {
            bool success;
            bool isprime= is_prime_trialdivision::call<PrimeTrialDivisionWarren,
                               TRIAL_DIVISION_SIZE>(static_cast<U>(x), success);
            if (success)
                return isprime;
            // is_prime_trialdivision::call should have successfully handled any
            // even x, and any x < 2.
            HPBC_ASSERT2(x % 2 != 0);
            HPBC_ASSERT2(x >= 2);
        } else {
            // Montgomery arithmetic requires odd numbers, and checking for
            // divisibility by two is trivial anyway.
            if (x % 2 == 0)
                return (x == 2);
            if (x < 2)
                return false;
        }
        HPBC_ASSERT2(x % 2 != 0);
        HPBC_ASSERT2(x >= 2);

        // At this point, we couldn't detect whether x is prime.  So we fall
        // back to determining primality via miller-rabin.

        if constexpr (std::is_same<typename MontgomeryForm<T>::MontyTag,
                                   TagMontyQuarterrange>::value) {
            return mont_miller_rabin(MontgomeryForm<T>(x));
        } else {
            static_assert(digitsU >= 2);
            constexpr U Rdiv4 = static_cast<U>(1) << (digitsU - 2);
            if (static_cast<U>(x) < Rdiv4)
                return mont_miller_rabin(MontgomeryQuarter<T>(x));
            else
                return mont_miller_rabin(MontgomeryForm<T>(x));
        }
    }
  }

}; // end struct impl_is_prime_intensive


}}  // end namespace

#endif
