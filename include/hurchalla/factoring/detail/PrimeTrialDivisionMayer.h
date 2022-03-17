// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_PRIME_TRIAL_DIVISION_MAYER_H_INCLUDED
#define HURCHALLA_FACTORING_PRIME_TRIAL_DIVISION_MAYER_H_INCLUDED


#include "hurchalla/factoring/detail/trial_divide_mayer.h"
#include "hurchalla/factoring/detail/OddPrimes.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <array>
#include <type_traits>

namespace hurchalla { namespace detail {


// See trial_divide_mayer.h for details of the algorithm used for isDivisible()
// and for a proof of its correctness.

template <typename T, int SIZE>
class PrimeTrialDivisionMayer {
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    static_assert(SIZE > 0);
public:
    using value_type = T;
    static constexpr int size() { return SIZE; }

    static constexpr T oddPrime(int index)
    {
        HPBC_CONSTEXPR_PRECONDITION(0 <= index && index < SIZE);
        return oddprimes[static_cast<std::size_t>(index)];
    }

    // Returns oddPrime() squared without overflow.
    static constexpr auto oddPrimeSquared(int index)
    {
        HPBC_CONSTEXPR_PRECONDITION(0 <= index && index < SIZE);
        U prime = oddprimes[static_cast<std::size_t>(index)];
        // get the smallest type that can always square an element of oddprimes
        // without overflow.
        using U2 =
              decltype(OddPrimes::get_constant_squared<U, oddprimes[SIZE-1]>());
        using P2 = typename safely_promote_unsigned<U2>::type;
        return static_cast<U2>(static_cast<P2>(prime) * static_cast<P2>(prime));
    }

    // Returns the first prime larger than the last prime used by this class
    static constexpr auto nextPrimePastEnd()
    {
        constexpr auto next = OddPrimes::get_next_prime<U, oddprimes[SIZE-1]>();
        return next;
    }

    // Returns nextPrimePastEnd() squared without overflow.
    static constexpr auto nextPrimePastEndSquared()
    {
        constexpr auto prime = nextPrimePastEnd();  // compile-time init
        // square the constant without overflow.
        constexpr auto square =
                      OddPrimes::get_constant_squared<decltype(prime), prime>();
        return square;
    }

    // Returns true if oddprimes[divisor_index] divides dividend and otherwise
    // returns false.
    // If oddprimes[divisor_index] divides dividend, then the quotient is placed
    // in quotient.  If oddprimes[divisor_index] does not divide dividend, then
    // the value of quotient is unspecified.
    static bool isDivisible(T& quotient, T dividend, int divisor_index)
    {
        HPBC_PRECONDITION2(0 <= divisor_index && divisor_index < SIZE);
        T divisor = oddprimes[static_cast<std::size_t>(divisor_index)];
        return trial_divide_mayer::call(quotient, dividend, divisor);
    }

private:
    // get the first N=SIZE odd primes
    static constexpr auto oddprimes = OddPrimes::get_array<SIZE>();
    using U = typename decltype(oddprimes)::value_type;
    static_assert(ut_numeric_limits<U>::is_integer);
    static_assert(!ut_numeric_limits<U>::is_signed);
    static_assert(std::is_same_v<
      const std::array<U,static_cast<std::size_t>(SIZE)>, decltype(oddprimes)>);
    // assert any element of the oddprimes array fits in type T
    static_assert(ut_numeric_limits<T>::digits >= ut_numeric_limits<U>::digits);
};


}}  // end namespace

#endif
