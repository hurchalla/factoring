// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FACTORIZE_TRIALDIVISION_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_TRIALDIVISION_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// factorize_trialdivision() partially (and sometimes completely)
// factorizes the number x using the first SIZE prime numbers as potential
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

// The function will trial a total of SIZE potential prime factors.  SIZE is a
// compile-time parameter for fine tuning of performance, that helps us to
// achieve the best possible overall performance of a factorization method in
// which this function is called.  Ideally any/all compile-time choices for SIZE
// would be guided by performance measurement tests.
//
// The function will set next_prime to the next prime larger than the last prime
// factor that it will potentially trial.  The value of next_prime is the same
// regardless of whether or not the function successfully finishes and returns
// before trialing all its potential factors.
//
// The default SIZE 54 tries all primes < 256 as potential divisors.  Thus this
// function with SIZE 54 could be used to guarantee complete factoring of any
// value x < 257*257, which includes all values of type uint16_t.

// The template-template param TTD should be either PrimeTrialDivisionWarren or
// PrimeTrialDivisionMayer.  PrimeTrialDivisionWarren is usually faster, but it
// uses around 5-10x more memory (PrimeTrialDivisionMayer typically uses around
// 2*SIZE bytes of memory).


// overload for uint8_t (not a partial specialization, which is impossible)
template <template<class,int> class TTD, int SIZE=54, class OutputIt>
OutputIt factorize_trialdivision(OutputIt iter,
                                 std::uint8_t& HURCHALLA_RESTRICT q,
                                 std::uint8_t& HURCHALLA_RESTRICT next_prime,
                                 std::uint8_t x)
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


template <template<class,int>class TTD, int SIZE=54, class OutputIt, typename T>
OutputIt factorize_trialdivision(OutputIt iter,
                                 T& HURCHALLA_RESTRICT q,
                                 T& HURCHALLA_RESTRICT next_prime,
                                 T x)
{
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    static_assert(SIZE > 1);
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations

    // PrimeTrialDivisionWarren and PrimeTrialDivisionMayer include only odd
    // primes - hence they don't use the prime 2, and so we set their TD_SIZE to
    // SIZE-1
    constexpr int TD_SIZE = SIZE-1;
    using TD = TTD<T, TD_SIZE>;
    static_assert(TD::oddPrime(0) == 3);

    constexpr auto tmp = TD::nextPrimePastEnd();
    // assert the next prime fits in type T
    static_assert(0 <= tmp && tmp <= ut_numeric_limits<T>::max());
    next_prime = static_cast<T>(tmp);
    constexpr auto next_prime_squared = TD::nextPrimePastEndSquared();

    // try the only even prime, 2, as a special case potential factor
    q = x;
    while (q % 2 == 0) {
        *iter++ = 2;
        q = static_cast<T>(q / 2);
        if (q == 1) return iter;  // we completely factored x
    }
    HPBC_ASSERT2(q > 1);

    for (int i=0; i<TD_SIZE; ++i) {
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


}}  // end namespace

#endif
