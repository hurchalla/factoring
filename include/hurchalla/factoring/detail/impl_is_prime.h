// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_trialdivision.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"

namespace hurchalla { namespace detail {


#ifndef HURCHALLA_ISPRIME_TRIALDIV_SIZE
// Some short perf testing on Haswell suggests 15 would be a good value for
// PrimeTrialDivisionMayer, and 54 a good value for PrimeTrialDivisionWarren.
// We use Mayer instead of Warren, since it's lower overhead on static memory
// and it works with the macro HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE, unlike
// Warren.  You can change the implementation to use PrimeTrialDivisionWarren
// if you want, by simply subsistuting that name for PrimeTrialDivisionMayer in
// the template function call below and including its header file.
// Note that this impl_is_prime() function (and the associated is_prime) is
// intended to be relatively lightweight, while in contrast IsPrimeIntensive is
// intended to be a more heavyweight option for repeated intensive primality
// testing.  That is why I chose for this lightweight function to use the
// lightweight PrimeTrialDivisionMayer rather than PrimeTrialDivisionWarren.

// FYI, size 54 would trial all prime factors < 256
#  define HURCHALLA_ISPRIME_TRIALDIV_SIZE (15)
#endif


// Note: we use a struct with static functions in order to disallow ADL
struct impl_is_prime {
  template <typename T>
  static bool call(T x)
  {
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    HPBC_PRECONDITION2(x >= 0);

    if constexpr (HURCHALLA_TARGET_BIT_WIDTH < ut_numeric_limits<T>::digits) {
        using U = typename sized_uint<HURCHALLA_TARGET_BIT_WIDTH>::type;
        static_assert(ut_numeric_limits<U>::is_integer);
        if (x <= ut_numeric_limits<U>::max())
            return impl_is_prime::call(static_cast<U>(x));
    }

    // First try small trial divisions to find easy factors.
    // If primality still unknown, use miller-rabin to prove prime or not prime.
    bool success;
    bool isPrime = is_prime_trialdivision::call<PrimeTrialDivisionMayer,
                                   HURCHALLA_ISPRIME_TRIALDIV_SIZE>(x, success);
    if (success)
        return isPrime;

    // At this point, we couldn't detect whether x is prime.  We'll fall back to
    // determining primality via miller-rabin.
    // is_prime_trialdivision should have handled evens and 0 and 1.
    HPBC_ASSERT2(x % 2 == 1);
    HPBC_ASSERT2(x > 1);
    return is_prime_miller_rabin::call(x);
  }
};


}}  // end namespace

#endif
