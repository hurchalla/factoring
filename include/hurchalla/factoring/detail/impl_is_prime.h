// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_wheel210.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// TODO: empirically find a decent value for max_trial_factor.
#ifndef HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR
#  define HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR (16 + 1*210)
#endif

template <typename T>
bool impl_is_prime(T x)
{
    HPBC_PRECONDITION2(x >= 0);
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    // First try small trial divisions to find easy factors or easy primality.
    // If primality still unknown, use miller-rabin to prove prime or composite.

    constexpr T max_trial_factor =
             (HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR <= ut_numeric_limits<T>::max())
             ? HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR : ut_numeric_limits<T>::max();
    bool isSuccessful;
    bool isPrime = is_prime_wheel210(x, &isSuccessful, max_trial_factor);
    if (isSuccessful)
        return isPrime;

    // At this point, we didn't find any factors and we couldn't detect whether
    // x is prime.  We'll fall back to determining primality via miller-rabin.
    return is_prime_miller_rabin_integral(x);
}


}}  // end namespace


#endif
