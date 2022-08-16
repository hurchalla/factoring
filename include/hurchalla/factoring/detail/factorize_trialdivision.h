// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FACTORIZE_TRIALDIVISION_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_TRIALDIVISION_H_INCLUDED


#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace detail {


// Do *NOT* change any values of the macros in the code immediately below.  If
// you wish to use different values for one or more of the macros (which is
// completely okay), please predefine the macro when compiling.
//
//
// HURCHALLA_TRIAL_DIVISION_TEMPLATE can be either PrimeTrialDivisionWarren or
// PrimeTrialDivisionMayer.  PrimeTrialDivisionWarren is usually faster, but it
// uses around 5-10x more memory (e.g. for T = uint64_t, 2.5KB vs. 0.25KB).
#ifndef HURCHALLA_TRIAL_DIVISION_TEMPLATE
#  define HURCHALLA_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionWarren
#endif
// FYI there are 54 primes <= 256.  Thus this function with SIZE 54 could be
// used to guarantee complete factoring of any value x < 257*257, which includes
// all values of type uint16_t.
// On benchmarks on Haswell, SIZE 139 worked well for PrimeTrialDivisionWarren
#ifndef HURCHALLA_TRIAL_DIVISION_CROSSOVER_BITS
#  define HURCHALLA_TRIAL_DIVISION_CROSSOVER_BITS 45
#endif
#ifndef HURCHALLA_TRIAL_DIVISION_SIZE_SMALL
#  define HURCHALLA_TRIAL_DIVISION_SIZE_SMALL 109
#endif
#ifndef HURCHALLA_TRIAL_DIVISION_SIZE_LARGE
#  define HURCHALLA_TRIAL_DIVISION_SIZE_LARGE 139
#endif


// Let LIMIT = (x < (1<<HURCHALLA_TRIAL_DIVISION_CROSSOVER_BITS))
//            ? HURCHALLA_TRIAL_DIVISION_SIZE_SMALL
//            : HURCHALLA_TRIAL_DIVISION_SIZE_LARGE
// factorize_trialdivision() partially (and sometimes completely)
// factorizes the number x using the first LIMIT prime numbers as potential
// divisors for divisibility testing.
//
// Preconditions for factorize_trialdivision():
//   Requires x >= 2.
//
// Postconditions:
// 1) The return value is an output iterator to the position one past the last
//   factor that the function wrote to the destination range (iterated by the
//   function's parameter 'iter').  The destination range consists of all the
//   prime factors that the function was able to find for x.  This range will
//   be empty if the function could not find any factors for x and could not
//   determine whether x was prime.  The range will consist of the single
//   element x if the function determined that x was prime.
// 2) q will be set to the quotient of x divided by all the elements written to
//   the destination range (iterated by the function's parameter 'iter').  If
//   nothing was written to the range (if it's empty), then q will be set to x.
//   There are specific details that naturally follow from these facts and from
//   Postcondition 1: if q gets set to 1, then this indicates the function was
//   able to completely factor x and the destination range consists of all the
//   factors.  If q > 1, then this indicates the function was not able to
//   completely factorize x, and q represents the value still remaining to be
//   factored.  q will never be set to zero (or a value < 0).

// The function will trial a total of LIMIT potential prime factors, unless it
// finishes early by completely factoring 'x'.
//
// The function will set next_prime to the next prime larger than the last prime
// factor that it will potentially trial.  The value of next_prime is the same
// regardless of whether or not the function successfully finishes and returns
// before trialing all its potential factors.


// Note: we use a struct with static functions in order to disallow ADL
struct factorize_trialdivision {
  // overload for uint8_t (not a partial specialization)
  template <class OutputIt>
  static OutputIt call(OutputIt iter, std::uint8_t& HURCHALLA_RESTRICT q,
           std::uint8_t& HURCHALLA_RESTRICT next_prime, std::uint8_t x, int = 0)
  {
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    using std::uint8_t;

    next_prime = 17;

    q = x;
    // we only need to try primes < 16, since 16*16==256 is greater than any
    // possible uint8_t value.
    while (q % 2 == 0) {
        *iter++ = 2;
        q = static_cast<uint8_t>(q / 2);
        if (q == 1) return iter;  // we completely factored x
    }
    while (q % 3 == 0) {
        *iter++ = 3;
        q = static_cast<uint8_t>(q / 3);
        if (q == 1) return iter;  // we completely factored x
    }
    while (q % 5 == 0) {
        *iter++ = 5;
        q = static_cast<uint8_t>(q / 5);
        if (q == 1) return iter;  // we completely factored x
    }
    while (q % 7 == 0) {
        *iter++ = 7;
        q = static_cast<uint8_t>(q / 7);
        if (q == 1) return iter;  // we completely factored x
    }
    while (q % 11 == 0) {
        *iter++ = 11;
        q = static_cast<uint8_t>(q / 11);
        if (q == 1) return iter;  // we completely factored x
    }
    while (q % 13 == 0) {
        *iter++ = 13;
        q = static_cast<uint8_t>(q / 13);
        if (q == 1) return iter;  // we completely factored x
    }
    HPBC_ASSERT2(q > 1);
    // q must be prime, because x < 256 (due to uint8_t) and thus q < 256, and
    // at this point q is not divisible by any primes < 16 == sqrt(256).
    *iter++ = q;
    q = 1;   // we completely factored x
    return iter;
  }


  template <class OutputIt, typename T>
  static OutputIt call(OutputIt iter, T& HURCHALLA_RESTRICT q,
                       T& HURCHALLA_RESTRICT next_prime,
                       T x, int prime_index = 0)
  {
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    constexpr int SIZE = HURCHALLA_TRIAL_DIVISION_SIZE_LARGE;
    static_assert(SIZE > 1);
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    int limit;
    {
        static_assert(HURCHALLA_TRIAL_DIVISION_SIZE_SMALL <= SIZE);
        // Since we're using uint64_t during the left shift below, we must shift
        // less than 64 bits (change to a bigger type if you need >= 64).
        static_assert(HURCHALLA_TRIAL_DIVISION_CROSSOVER_BITS < 64);
        constexpr uint64_t crossover =
            static_cast<uint64_t>(1) << HURCHALLA_TRIAL_DIVISION_CROSSOVER_BITS;
        limit = (x < crossover) ? HURCHALLA_TRIAL_DIVISION_SIZE_SMALL : SIZE;
        HPBC_ASSERT2(limit <= SIZE);
    }

    // HURCHALLA_TRIAL_DIVISION_TEMPLATE can be PrimeTrialDivisionWarren or
    // PrimeTrialDivisionMayer.  These classes include only odd primes -
    // hence they don't use the prime 2, and so we set their TD_SIZE to SIZE-1.
    constexpr int TD_SIZE = SIZE-1;
    using TD = HURCHALLA_TRIAL_DIVISION_TEMPLATE<T, TD_SIZE>;
    static_assert(TD::oddPrime(0) == 3);

    int td_limit = limit - 1;

    constexpr auto tmp = TD::nextPrimePastEnd();
    // assert the next prime fits in type T
    static_assert(0 <= tmp && tmp <= ut_numeric_limits<T>::max());
    // (paranoia) use std::integral_constant to guarantee compile-time values.
    // Constexpr guarantees only that a variable *can* be constant evaluated (at
    // compile-time), but technically, doesn't guarantee that it *will* be, when
    // there's nothing to prevent run-time eval.
    next_prime = static_cast<T>(
                             std::integral_constant<decltype(tmp), tmp>::value);
    constexpr auto tmp2 = TD::nextPrimePastEndSquared();
    // (paranoia) use std::integral_constant to guarantee a compile-time value.
    constexpr auto next_prime_squared =
                            std::integral_constant<decltype(tmp2), tmp2>::value;

    // try the only even prime, 2, as a special case potential factor
    q = x;
    while (q % 2 == 0) {
        *iter++ = 2;
        q = static_cast<T>(q / 2);
        if (q == 1) return iter;  // we completely factored x
    }
    HPBC_ASSERT2(q > 1);

    if constexpr (ut_numeric_limits<T>::digits > HURCHALLA_TARGET_BIT_WIDTH) {
        if ((q >> HURCHALLA_TARGET_BIT_WIDTH) == 0) {
            using T2 = sized_uint<HURCHALLA_TARGET_BIT_WIDTH>::type;
            T2 q2 = static_cast<T2>(q);
            T2 next_prime2;
            iter = call<SIZE>(iter, q2, next_prime2, q2, 0);
            q = q2;
            next_prime = next_prime2;
            return iter;
        }
    }

    for (int i=prime_index; i<td_limit; ++i) {
        HPBC_ASSERT2(q > 1);
        if (TD::oddPrimeSquared(i) > q) {
            // Since no primes <= sqrt(q) are factors of q, q must be a prime
            *iter++ = q;
            q = 1;  // we completely factored x
            return iter;
        }
        T div_result;
        // next line is roughly equivalent to  while (q % TD::oddPrime(i) == 0)
        while (TD::isDivisible(div_result, q, i)) {
            const T prime = TD::oddPrime(i);
            // Since we reached this point, we know (q % prime == 0).
            *iter++ = prime;
            // Since the prime divides q, div_result was set to q/prime.
            q = div_result;
            if (q == 1) return iter;  // we completely factored x

            if constexpr (ut_numeric_limits<T>::digits
                          > HURCHALLA_TARGET_BIT_WIDTH) {
                if ((q >> HURCHALLA_TARGET_BIT_WIDTH) == 0) {
                    using T2 = sized_uint<HURCHALLA_TARGET_BIT_WIDTH>::type;
                    T2 q2 = static_cast<T2>(q);
                    T2 next_prime2;
                    iter = call<SIZE>(iter, q2, next_prime2, q2, i);
                    q = q2;
                    next_prime = next_prime2;
                    return iter;
                }
            }
        }
        HPBC_ASSERT2(q > 1);
    }
    HPBC_ASSERT2(q > 1);

    if (q < next_prime_squared) {
        // we know q is prime, because we tested all possible factors for q.
        *iter++ = q;
        q = 1;  // we completely factored x
    }
    return iter;
  }
}; // end struct factorize_trialdivision


}}  // end namespace

#endif
