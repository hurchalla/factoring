// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_FACTORIZE_DISPATCH_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_DISPATCH_H_INCLUDED


#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/factorize_trialdivision.h"
#include "hurchalla/factoring/detail/factorize_pollard_rho.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


#define HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionWarren
//#define HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionMayer

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

    T q, next_prime;
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

    iter = factorize_pollard_rho(iter, q, is_prime_mf, threshold_always_prime);
    return iter;
}


}} // end namespace

#endif
