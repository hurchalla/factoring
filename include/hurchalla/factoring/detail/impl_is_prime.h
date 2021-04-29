// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_trialdivision.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// TODO: empirically find good default(s) for HURCHALLA_ISPRIME_TRIALDIV_SIZE

#ifndef HURCHALLA_ISPRIME_TRIALDIV_SIZE
// Some short perf testing on Haswell suggest 15 would be a good value for
// PrimeTrialDivisionMayer, and 54 a good value for PrimeTrialDivisionWarren.
// We use Mayer instead of Warren, since it's lower overhead on static memory
// and it works with the macro HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE, unlike
// Warren.
// FYI, size 54 would trial all prime factors < 256
#  define HURCHALLA_ISPRIME_TRIALDIV_SIZE (15)
#endif


template <typename T>
bool impl_is_prime(T x)
{
    HPBC_PRECONDITION2(x >= 0);
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    // First try small trial divisions to find easy factors.
    // If primality still unknown, use miller-rabin to prove prime or not prime.
    bool success;
    bool isPrime = is_prime_trialdivision<PrimeTrialDivisionMayer,
                                   HURCHALLA_ISPRIME_TRIALDIV_SIZE>(x, success);
    if (success)
        return isPrime;

    // At this point, we couldn't detect whether x is prime.  We'll fall back to
    // determining primality via miller-rabin.
    return is_prime_miller_rabin_integral(x);
}


}}  // end namespace

#endif
