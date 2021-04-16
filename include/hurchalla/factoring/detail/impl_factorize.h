// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_IMPL_FACTORIZE_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_FACTORIZE_H_INCLUDED


#include "hurchalla/factoring/detail/TrialDivisionWarren.h"
#include "hurchalla/factoring/detail/TrialDivisionMayer.h"

#include "hurchalla/factoring/detail/factorize_trialdivision.h"
#include "hurchalla/factoring/detail/pollard_rho_factorize.h"
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


#define TRIAL_DIVISION_TEMPLATE TrialDivisionWarren
//#define TRIAL_DIVISION_TEMPLATE TrialDivisionMayer


// Initial tests so far suggest it might be fastest to entirely skip wheel
// factorization prior to pollard rho.
// A client should enable a wheel factorization stage (if desired) by
// predefining the following macro, though we can define it here for quick and
// dirty testing.
//#define HURCHALLA_USE_WHEEL_FACTORIZATION


#define HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048

// estimate min at 318
// 1165 2541
// estimate min at 320
// 1147 2540

#define HURCHALLA_NUM_PRIMES_UNDER_2048 309

//#define HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048_INDEX_LIMIT 310
// 1138 2538
//#define HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048_INDEX_LIMIT 200
// 1157 2537
//#define HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048_INDEX_LIMIT 250
// 1143 2537
#define HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048_INDEX_LIMIT 280
// 1139 2537
// 1142 2540

// 1098 2536

// 1374 2548

// 1101 2537


template <template<class> class Functor, class OutputIt, typename T>
T impl_factorize(OutputIt iter, T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits % 2 == 0, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    constexpr T sqrtR = static_cast<T>(1)<<(ut_numeric_limits<T>::digits/2);

    T q;
#ifdef HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048
    T next_prime;
    int index_limit = HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048_INDEX_LIMIT;

    iter = factorize_trialdivision<TRIAL_DIVISION_TEMPLATE,
                                   HURCHALLA_NUM_PRIMES_UNDER_2048>
                                          (iter, q, next_prime, x, index_limit);
    HPBC_ASSERT2(q >= 1);  // small_trial_division2048 guarantees this
    if (q == 1)   // if small_trial_division2048 completely factored x
        return 0;
    // small_trial_division2048() guarantees that any factor of x that is less
    // than next_prime*next_prime must be prime.
    T threshold_always_prime = (next_prime < sqrtR) ?
          static_cast<T>(next_prime * next_prime) : ut_numeric_limits<T>::max();
#else
#ifdef HURCHALLA_USE_WHEEL_FACTORIZATION
    // Ensure max_trial_factor * max_trial_factor never overflows.
    // Note that no factors ever exist >= sqrtR, so limiting to sqrtR-1 is fine.
    constexpr T max_trial_factor =
           (HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR < sqrtR)
           ? HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR : static_cast<T>(sqrtR - 1);
    iter = factorize_wheel210(iter, q, x, max_trial_factor);
    HPBC_ASSERT2(q >= 1);  // factorize_wheel210 guarantees this
    if (q == 1)   // if factorize_wheel210 completely factored x
        return 0;
    constexpr T threshold_always_prime =
                              static_cast<T>(max_trial_factor*max_trial_factor);
#else
    // factor out all primes < 256
    iter = small_trial_division256(iter, q, x);
    HPBC_ASSERT2(q >= 1);  // small_trial_division256 guarantees this
    if (q == 1)   // if small_trial_division256 completely factored x
        return 0;
    // small_trial_division256() guarantees that any factor of x that is less
    // than 257*257 must be prime.
    constexpr T threshold_always_prime = (257<sqrtR) ? static_cast<T>(257*257) :
                                                    ut_numeric_limits<T>::max();
#endif
#endif

    T iterations_performed;
    iter = pollard_rho_factorize<TRIAL_DIVISION_TEMPLATE, Functor>(iter, q,
              threshold_always_prime, static_cast<T>(1), &iterations_performed);
    return iterations_performed;
}


}} // end namespace

#endif
