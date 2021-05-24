// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IS_PRIME_TRIALDIVISION_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_TRIALDIVISION_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace detail {


// is_prime_trialdivision(): attempts to determine primality for the number x,
// by trial dividing x by the first SIZE prime numbers.
//
// Preconditions for is_prime_trialdivision():
//   Requires SIZE >= 2
//
// Postconditions:
//   The function returns false if it finds any factor for x (and sets
// is_successful=true), or returns true if it determines that x is prime (and
// sets is_successful=true).  If the function is unable to determine primality
// for the input x then it sets is_successful=false and returns an undefined
// boolean value (either true or false).
//
// Template and Function Parameters:
//   The function will trial a total of SIZE potential prime factors.  SIZE is
// a compile-time parameter for fine tuning of performance, which helps you to
// achieve the best possible overall performance when combining this function
// with additional primality tests.  Ideally any/all choices for SIZE should be
// guided by performance measurements.
//   The default SIZE 54 tries all primes < 256 as potential divisors.  Thus
// this function with SIZE 54 could be used to guarantee success at primality
// testing any value x < 257*257, which includes all values of type uint16_t.
//   Template-template parameter TTD should be either PrimeTrialDivisionWarren
// or PrimeTrialDivisionMayer.  PrimeTrialDivisionWarren is usually faster, but
// it uses around 5-10x more memory (PrimeTrialDivisionMayer typically uses
// around 2*SIZE bytes of memory).


// (The implementations in ths file are adapted from factorize_trialdivision())

// overload for uint8_t (not a partial specialization, which is impossible)
template <template<class,int> class TTD, int SIZE=54>
bool is_prime_trialdivision(std::uint8_t x, bool& is_successful)
{
    using std::uint8_t;
    is_successful = true;  // this function overload will always succeed
    if (x < 2)
        return false;
    // we only need to try primes < 16, since 16*16==256 is greater than any
    // possible uint8_t value.
    if (x % 2 == 0)
        return (x == 2);
    if (x % 3 == 0)
        return (x == 3);
    if (x % 5 == 0)
        return (x == 5);
    if (x % 7 == 0)
        return (x == 7);
    if (x % 11 == 0)
        return (x == 11);
    if (x % 13 == 0)
        return (x == 13);
    // x must be prime, because we know x < 256 (due to uint8_t) and x was
    // not divisible by any primes < 16 == sqrt(256)
    return true;
}

template <template<class,int> class TTD, int SIZE=54, typename T>
bool is_prime_trialdivision(T x, bool& is_successful)
{
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    static_assert(SIZE > 1);

    // PrimeTrialDivisionWarren and PrimeTrialDivisionMayer include only odd
    // primes - hence they don't use the prime 2, and so we set their TD_SIZE to
    // SIZE-1
    constexpr int TD_SIZE = SIZE-1;
    using TD = TTD<T, TD_SIZE>;  
    static_assert(TD::oddPrime(0) == 3);

    is_successful = true;
    if (x < 2)
        return false;
    // try the only even prime, 2, as a special case potential factor
    if (x % 2 == 0)
        return (x == 2);

    constexpr auto next_prime_squared = TD::nextPrimePastEndSquared();
    for (int i=0; i<TD_SIZE; ++i) {
        if (TD::oddPrimeSquared(i) > x) {
            // Since no primes <= sqrt(x) were factors of x, x must be a prime
            return true;
        }
        T div_result;
        // next line is roughly equivalent to  if (x % TD::oddPrime(i) == 0)
        if (TD::isDivisible(div_result, x, i))
            return false;
    }

    if (x < next_prime_squared) {
        // we know x is prime, because we tested all possible factors for x.
        return true;
    }
    else {
        // We weren't able to determine if x is prime.
        is_successful = false;
        return false;  // it doesn't matter what bool value we return.
    }
}


}}  // end namespace

#endif
