// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_FACTORIZE_DISPATCH_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_DISPATCH_H_INCLUDED


#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/factorize_trialdivision.h"
#include "hurchalla/factoring/detail/factorize_pollard_rho.h"
#include "hurchalla/factoring/detail/factorize_wheel210.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// TODO determine good max_factor
#ifndef HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR
// 256 + 0*210  results in skipping factorize_wheel210(), and only using
// small_division256().  We want to use 256 + some multiple of 210, since that
// makes an optimal max trial factor for factorize_wheel210.
#  define HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR (256 + 0*210)
#endif


#define HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionWarren
//#define HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionMayer


// Initial tests so far suggest it might be fastest to always define
// HURCHALLA_USE_PR_TRIAL_DIVISION.  If it is not defined, wheel factorization
// will be used.  But it's possible wheel factorization might even be slower
// than doing almost nothing prior to Pollard Rho?
#define HURCHALLA_USE_PR_TRIAL_DIVISION

// FYI there are 54 primes below 256
#define HURCHALLA_PR_TRIAL_DIVISION_SIZE 135

// I'll probably want to get rid of this macro entirely, and also change
// the function not to take an index limit argument
#define HURCHALLA_PR_TRIAL_DIVISION_INDEX_LIMIT HURCHALLA_PR_TRIAL_DIVISION_SIZE



template <class OutputIt, typename T, class PrimalityFunctor>
OutputIt factorize_dispatch(OutputIt iter, T x, const PrimalityFunctor& is_prime_mf)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits % 2 == 0, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    constexpr T sqrtR = static_cast<T>(1)<<(ut_numeric_limits<T>::digits/2);

    T q;
#ifdef HURCHALLA_USE_PR_TRIAL_DIVISION
    T next_prime;
    int index_limit = HURCHALLA_PR_TRIAL_DIVISION_INDEX_LIMIT;

    iter = factorize_trialdivision<HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE,
                                   HURCHALLA_PR_TRIAL_DIVISION_SIZE>
                                          (iter, q, next_prime, x, index_limit);
    HPBC_ASSERT2(q >= 1);  // factorize_trialdivision() guarantees this
    if (q == 1)   // if factorize_trialdivision() completely factored x
        return iter;
    // factorize_trialdivision() guarantees that any factor of x that is less
    // than next_prime*next_prime must be prime.
    T threshold_always_prime = (next_prime < sqrtR) ?
          static_cast<T>(next_prime * next_prime) : ut_numeric_limits<T>::max();
#else
    // Use Wheel Factorization--
    // Ensure max_trial_factor * max_trial_factor never overflows.
    // Note that no factors ever exist >= sqrtR, so limiting to sqrtR-1 is fine.
    constexpr T max_trial_factor =
           (HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR < sqrtR)
           ? HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR : static_cast<T>(sqrtR - 1);
    iter = factorize_wheel210(iter, q, x, max_trial_factor);
    HPBC_ASSERT2(q >= 1);  // factorize_wheel210 guarantees this
    if (q == 1)   // if factorize_wheel210 completely factored x
        return iter;
    constexpr T threshold_always_prime =
                              static_cast<T>(max_trial_factor*max_trial_factor);
#endif

    iter = factorize_pollard_rho(iter, q, is_prime_mf, threshold_always_prime);
    return iter;
}


}} // end namespace

#endif
