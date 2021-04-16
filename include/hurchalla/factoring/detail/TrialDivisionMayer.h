// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_TRIAL_DIVISION_MAYER_H_INCLUDED
#define HURCHALLA_FACTORING_TRIAL_DIVISION_MAYER_H_INCLUDED


#include "hurchalla/factoring/detail/trial_divide_mayer.h"
#include "hurchalla/factoring/detail/odd_primes.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <array>
#include <type_traits>

namespace hurchalla { namespace detail {


// See trial_divide_mayer.h for details of the algorithm used here and for a
// proof of its correctness.

template <typename T, int SIZE>
class TrialDivisionMayer {
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    static_assert(SIZE > 0);
public:
    using value_type = T;
    static constexpr int size() { return SIZE; }

    static constexpr T oddPrime(int divisor_index)
    {
        HPBC_CONSTEXPR_PRECONDITION(0 <= divisor_index && divisor_index < SIZE);
        return oddprimes[static_cast<std::size_t>(divisor_index)];
    }

    // Returns oddPrime() squared without overflow.
    static constexpr auto oddPrimeSquared(int divisor_index)
    {
        HPBC_CONSTEXPR_PRECONDITION(0 <= divisor_index && divisor_index < SIZE);
        U prime = oddprimes[static_cast<std::size_t>(divisor_index)];
        // get the smallest type that can always square an element of oddprimes
        // without overflow.
        using U2 = decltype(get_constant_squared<U, oddprimes[SIZE-1]>());
        using P2 = typename safely_promote_unsigned<U2>::type;
        return static_cast<U2>(static_cast<P2>(prime) * static_cast<P2>(prime));
    }

    // Returns the first prime larger than the last prime used by this class
    static constexpr auto nextPrimePastEnd()
    {
        constexpr auto next = get_next_prime<U, oddprimes[SIZE-1]>();
        return next;
    }

    // Returns nextPrimePastEnd() squared without overflow.
    static constexpr auto nextPrimePastEndSquared()
    {
        constexpr auto prime = nextPrimePastEnd();  // compile-time init
        // square the constant without overflow.
        constexpr auto square = get_constant_squared<decltype(prime), prime>();
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
        return trial_divide_mayer(quotient, dividend, divisor);
    }

private:
    // get the first N=SIZE odd primes
    static constexpr auto oddprimes = get_odd_primes<SIZE>();
    using U = typename decltype(oddprimes)::value_type;
    static_assert(std::is_same_v<
      const std::array<U,static_cast<std::size_t>(SIZE)>, decltype(oddprimes)>);
    // assert any element of the oddprimes array fits in type T
    static_assert(ut_numeric_limits<T>::max() >= ut_numeric_limits<U>::max());
};


}}  // end namespace

#endif
