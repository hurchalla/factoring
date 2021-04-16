// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_FACTORIZE_TRIALDIVISION_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_TRIALDIVISION_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// overload for uint8_t (not a partial specialization, which is impossible)
template <template<class,int> class TTD, int SIZE=54, class OutputIt>
OutputIt factorize_trialdivision(OutputIt iter, std::uint8_t& q,
    std::uint8_t& next_prime, std::uint8_t x, int)
{
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    using std::uint8_t;

    next_prime = 0;  // we'll always completely factor x, so any value would do

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




// factorize_trialdivision() partially (and sometimes completely)
// factorizes the number x using the first SIZE or size_limit (whichever is
// smaller) prime numbers as potential divisors for divisibility testing.
//
// Preconditions for factorize_trialdivision():
//   Requires x >= 2.
//
// Postconditions:
// (note these are the same as for factorize_wheel210())
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

// The function will trial a total of either SIZE or size_limit potential prime
// factors, whichever is smaller.  SIZE is a compile-time parameter for fine
// tuning of performance, that helps us to achieve the best possible overall
// performance of a factorization method in which this function is called.
// Likewise, size_limit is a run-time parameter intended for fine tuning, given
// a particular number to factor.  The idea behind size_limit is that when
// factoring smaller numbers it *might* improve performance to use size_limit
// to stop trial division before trialing all SIZE primes.  Ideally any/all
// compile-time and run-time choices for SIZE and size_limit would be guided by
// performance measurement tests.  
//
// If the function does not complete factorize x, it will copy into next_prime
// the next prime larger than the last prime factor that it trialed.  For some
// factorization methods you may want/need to know the next prime that wasn't
// yet trialed in order to continue factoring, while for others methods (e.g
// pollard-rho) it's irrelevant.  The value of next_prime is unspecified if the
// function completely factors x.
//
// The default SIZE 54 tries all primes < 256 as potential divisors (unless
// size_limit < SIZE).  Thus this function with SIZE 54 could be used to
// guarantee complete factoring of any value x < 257*257, which includes all
// values of type uint16_t.

// the template-template param TTD should be either TrialDivisionWarren or
// TrialDivisionMayer.
template <template<class,int> class TTD,
          int SIZE=54, class OutputIt, typename T>
OutputIt factorize_trialdivision(OutputIt iter, T& q, T& next_prime, T x,
                                                                 int size_limit)
{
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    static_assert(SIZE > 1);
    HPBC_PRECONDITION2(size_limit > 0);
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations

    // TrialDivisionWarren and TrialDivisionMayer include only odd primes -
    // hence they don't use the prime 2, and so we set their TD_SIZE to SIZE-1
    constexpr int TD_SIZE = SIZE-1;  
    int td_size_limit = size_limit - 1;  // same as above

    using TD = TTD<T, TD_SIZE>;  
    static_assert(TD::oddPrime(0) == 3);

    // try the only even prime, 2, as a special case potential factor
    q = x;
    while (q % 2 == 0) {
        *iter++ = 2;
        q = static_cast<T>(q / 2);
        if (q == 1) return iter;  // we completely factored x
    }
    HPBC_ASSERT2(q > 1);

    using U = decltype(TD::nextPrimePastEndSquared());
    U next_prime_squared;
    if (td_size_limit >= TD_SIZE) {
        td_size_limit = TD_SIZE;
        constexpr auto tmp = TD::nextPrimePastEnd();
        static_assert(ut_numeric_limits<decltype(tmp)>::is_integer);
        // assert the next prime fits in type T
        static_assert(0 <= tmp && tmp <= ut_numeric_limits<T>::max());
        next_prime = static_cast<T>(tmp);
        next_prime_squared = TD::nextPrimePastEndSquared();
    }
    else {
        next_prime = TD::oddPrime(td_size_limit);
        auto tmp = TD::oddPrimeSquared(td_size_limit);
        static_assert(ut_numeric_limits<decltype(tmp)>::is_integer);
        static_assert(!ut_numeric_limits<decltype(tmp)>::is_signed);
        static_assert(ut_numeric_limits<U>::max() >=
                      ut_numeric_limits<decltype(tmp)>::max());
        next_prime_squared = tmp;
    }
    HPBC_ASSERT2(td_size_limit <= TD_SIZE);

    for (int i=0; i<td_size_limit; ++i) {
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
