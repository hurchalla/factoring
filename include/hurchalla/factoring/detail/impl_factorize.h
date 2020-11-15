// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_IMPL_FACTORIZE_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_FACTORIZE_H_INCLUDED


#include "hurchalla/factoring/detail/pollard_rho_factorize.h"
#include "hurchalla/factoring/detail/wheel_factorization210.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace factoring {


// TODO determine good max_factor
#ifndef HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR
// 256 + 0*210  results in skipping wheel_factorization210(), and only using
// small_division256().  We want to use 256 + some multiple of 210, since that
// makes an optimal max trial factor for wheel_factorization210.
#  define HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR (256 + 0*210)
#endif


// Initial tests so far suggest it might be fastest to entirely skip wheel
// factorization prior to pollard rho.
// A client should enable a wheel factorization stage (if desired) by
// predefining the following macro, though we can define it here for quick and
// dirty testing.
//#define HURCHALLA_USE_WHEEL_FACTORIZATION



template <template<class> class Functor, class OutputIt, typename T>
T impl_factorize(OutputIt iter, T x)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");
    static_assert(ma::ma_numeric_limits<T>::digits % 2 == 0, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations

    T q;
#ifdef HURCHALLA_USE_WHEEL_FACTORIZATION
    constexpr T sqrtR = static_cast<T>(1)<<(ma::ma_numeric_limits<T>::digits/2);
    // Ensure max_trial_factor * max_trial_factor never overflows.
    // Note that no factors ever exist >= sqrtR, so limiting to sqrtR-1 is fine.
    constexpr T max_trial_factor =
          (HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR < sqrtR)
          ? HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR : static_cast<T>(sqrtR - 1);
    iter = wheel_factorization210(iter, q, x, max_trial_factor);
    HPBC_ASSERT2(q >= 1);  // wheel_factorization210 guarantees this
    if (q == 1)   // if wheel_factorization210 completely factored x
        return 0;
    T threshold_always_prime= static_cast<T>(max_trial_factor*max_trial_factor);
#else
    // factor out all primes < 256
    iter = small_trial_division256(iter, q, x);
    HPBC_ASSERT2(q >= 1);  // small_trial_division256 guarantees this
    if (q == 1)   // if small_trial_division256 completely factored x
        return 0;
    // small_trial_division256() guarantees that any factor of x that is less
    // than 257*257 must be prime.
    T threshold_always_prime = static_cast<T>(257 * 257);
#endif

    T iterations_performed;
    iter = pollard_rho_factorize<Functor>(iter, q, threshold_always_prime,
                                      static_cast<T>(1), &iterations_performed);
    return iterations_performed;
}


}} // end namespace

#endif
