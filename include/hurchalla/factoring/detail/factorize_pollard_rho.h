// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FACTORIZE_POLLARD_RHO_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_POLLARD_RHO_H_INCLUDED


#include "hurchalla/factoring/detail/PollardRhoBrentSwitchingTrial.h"
#ifndef HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME
#  define HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME PollardRhoBrentSwitchingTrial
#else
#  include "hurchalla/factoring/detail/PollardRhoBrentTrial.h"
#  include "hurchalla/factoring/detail/PollardRhoBrentTrialParallel.h"
#  include "hurchalla/factoring/detail/experimental/PollardRhoTrial.h"
#  include "hurchalla/factoring/detail/experimental/PollardRhoBrentMontgomeryTrial.h"
#  include "hurchalla/factoring/detail/experimental/PollardRhoBrentMontgomeryTrialParallel.h"
#endif
#include "hurchalla/factoring/detail/factorize_wheel210.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unreachable.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// defined lower in this file
template <int log2ModulusLimit, class MF, class PrimalityFunctor,
          class OutputIt, typename T>
OutputIt hurchalla_factorize_pr_internal(OutputIt iter, T x,
                  const PrimalityFunctor& is_prime_pr, T threshold_always_prime,
                  T base_c, T expected_iterations);


#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#  pragma warning(disable : 4702)
#endif

// Dispatch function to the fastest  hurchalla_factorize_pr_internal  template
// function instantiation available for x (parameterized on the fastest
// MontgomeryForm for x).

// Note: we use a struct with static functions in order to disallow ADL
struct factorize_pollard_rho {
  template <int log2ModulusLimit, class PrimalityFunctor,
            class OutputIt, typename T>
  static OutputIt call(OutputIt iter, T x, const PrimalityFunctor& is_prime_pr,
          T threshold_always_prime = 0, T base_c = 1, T expected_iterations = 0)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    HPBC_PRECONDITION2(x % 2 == 1);  // x odd is required for montgomery form

    using S = sized_uint<HURCHALLA_TARGET_BIT_WIDTH>::type;

    // factor using a native integer type or smaller, if possible
    if (ut_numeric_limits<T>::digits <= ut_numeric_limits<S>::digits ||
                                             x <= ut_numeric_limits<S>::max()) {
        using U = typename std::conditional<
                  (ut_numeric_limits<T>::digits < ut_numeric_limits<S>::digits),
                  T, S>::type;
#if defined(HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH)
        using MF = MontgomeryStandardMathWrapper<U>;
        return hurchalla_factorize_pr_internal<log2ModulusLimit, MF>(
                                   iter, x, is_prime_pr, threshold_always_prime,
                                   base_c, expected_iterations);
#else
        static_assert(ut_numeric_limits<U>::digits >= 2);
        constexpr U Udiv4 = static_cast<U>(static_cast<U>(1) <<
                                            (ut_numeric_limits<U>::digits - 2));
        constexpr U Udiv2 = static_cast<U>(static_cast<U>(1) <<
                                            (ut_numeric_limits<U>::digits - 1));
        if (x < Udiv4) {
            using MF = MontgomeryQuarter<U>;
            return hurchalla_factorize_pr_internal<log2ModulusLimit, MF>(
                                   iter, x, is_prime_pr, threshold_always_prime,
                                   base_c, expected_iterations);
        } else if (x < Udiv2) {
            using MF = MontgomeryHalf<U>;
            return hurchalla_factorize_pr_internal<log2ModulusLimit, MF>(
                                   iter, x, is_prime_pr, threshold_always_prime,
                                   base_c, expected_iterations);
        } else {
            constexpr bool allowHalfRange =
                       (log2ModulusLimit == ut_numeric_limits<T>::digits - 1);
            // allowHalfRange applies to T, so we can only use MontgomeryHalf on
            // U if T and U are the same size.  Note that if U is smaller than
            // the native bit width, then MontgomeryForm<U> can sometimes map to
            // a more efficient MontyType than MontgomeryHalf<U>.  And if U is
            // larger than the native bit width (though it won't be), we'd again
            // prefer MontgomeryForm (it tends to work better at huge sizes).
            // Thus we only use MontgomeryHalf if U is the same size as the
            // native bit width, and also the same size as T, and allowHalfRange
            // is true.
            using MF = typename std::conditional<(allowHalfRange &&
                 ut_numeric_limits<U>::digits == ut_numeric_limits<T>::digits &&
                 ut_numeric_limits<U>::digits == HURCHALLA_TARGET_BIT_WIDTH),
                 MontgomeryHalf<U>, MontgomeryForm<U>>::type;
            return hurchalla_factorize_pr_internal<log2ModulusLimit, MF>(
                                   iter, x, is_prime_pr, threshold_always_prime,
                                   base_c, expected_iterations);
        }
#endif
    }
    // If we reach this point, the following clause will be true.  We explicitly
    // use an 'if' anyway, so that we can be sure the compiler will not generate
    // any code for it when T's bit width <= HURCHALLA_TARGET_BIT_WIDTH.
    if constexpr (ut_numeric_limits<T>::digits > HURCHALLA_TARGET_BIT_WIDTH) {
#if defined(HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH)
        using MF = MontgomeryStandardMathWrapper<T>;
        return hurchalla_factorize_pr_internal<log2ModulusLimit, MF>(
                                   iter, x, is_prime_pr, threshold_always_prime,
                                   base_c, expected_iterations);
#else
        static_assert(ut_numeric_limits<T>::digits >= 2);
        T Rdiv4 = static_cast<T>(
                       static_cast<T>(1) << (ut_numeric_limits<T>::digits - 2));
        if (x < Rdiv4) {
            using MF = MontgomeryQuarter<T>;
            return hurchalla_factorize_pr_internal<log2ModulusLimit, MF>(
                                   iter, x, is_prime_pr, threshold_always_prime,
                                   base_c, expected_iterations);
        } else {
            using MF = MontgomeryForm<T>;
            return hurchalla_factorize_pr_internal<log2ModulusLimit, MF>(
                                   iter, x, is_prime_pr, threshold_always_prime,
                                   base_c, expected_iterations);
        }
#endif
    } else {
        ::hurchalla::unreachable();
        HPBC_ASSERT2(false);  // it ought to be impossible to reach this code.
        return iter;
    }
#if defined(__INTEL_COMPILER)
    // for some reason icc doesn't see that all paths return prior to this point
    // so we'll uselessly return iter to quiet icc's warning/error
    return iter;
#endif
  }
};

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif



template <int log2ModulusLimit, class MF, class PrimalityFunctor,
          class OutputIt, typename T>
OutputIt hurchalla_factorize_pr_internal(OutputIt iter, T x,
                  const PrimalityFunctor& is_prime_pr, T threshold_always_prime,
                  T base_c, T expected_iterations)
{
    using S = typename MF::IntegerType;
    using C = typename MF::CanonicalValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<S>::is_integer, "");
    static_assert(!ut_numeric_limits<S>::is_signed, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    HPBC_PRECONDITION2(x % 2 == 1);  // odd x is required for montgomery form
    HPBC_PRECONDITION2(x <= ut_numeric_limits<S>::max());

    if (x < threshold_always_prime) {
        *iter++ = x;  // x is prime
        return iter;
    }
    MF mf(static_cast<S>(x));

    if (is_prime_pr(mf)) {
        *iter++ = x;  // x is prime
        return iter;
    }
    using P = typename safely_promote_unsigned<T>::type;

    HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME<MF> pr_trial_func;

    // we don't want to use a sequence  x[i+1] = x[i]*x[i] + c  where c == 0 or
    // c == -2.  See JM Pollard "A Monte Carlo method for factorization".  The
    // following only guarantees we avoid those sequences on the first iteration
    // of the loop, but unless x is very small, it's extremely unlikely next_c
    // would ever grow large enough during loop iterations to reach x-2.  If we
    // do end up with one of those sequences, it's fairly harmless - they just
    // have low success rate, making them inefficient since we may need another
    // loop iteration after it.
    if (base_c >= x-2 || base_c == 0)
        base_c = 1;
    C unity = mf.getUnityValue();

    C cc = mf.getCanonicalValue(mf.convertIn(static_cast<S>(base_c)));
    for (T i = 0; i < x; ++i) {
        HPBC_ASSERT2(expected_iterations <= ut_numeric_limits<S>::max());
        S tmp_expected = static_cast<S>(expected_iterations);
        T tmp_factor = static_cast<T>(pr_trial_func(mf, tmp_expected, cc));
        HPBC_ASSERT2(tmp_expected <= ut_numeric_limits<T>::max());
        expected_iterations = static_cast<T>(tmp_expected);
        if (tmp_factor >= 2) {    // we found a good factor (maybe prime).
            // Next_c could overflow, but that's okay.  We'd prefer for
            // efficiency that it didn't, but any T value would be valid.
            T next_c = static_cast<T>(base_c + static_cast<P>(i) + 1);
            // Try to factor the factor (it may or may not be prime)
            iter = factorize_pollard_rho::call<log2ModulusLimit>(
                          iter, tmp_factor, is_prime_pr, threshold_always_prime,
                          next_c, expected_iterations);
            // Next try to factor the quotient.
            // since 1 < tmp_factor < x, we know 1 < (x/tmp_factor) < x
            T quotient = static_cast<T>(x/tmp_factor);  
            iter = factorize_pollard_rho::call<log2ModulusLimit>(
                          iter, quotient, is_prime_pr, threshold_always_prime,
                          next_c, expected_iterations);
            return iter;
        }
        else {
            // tmp_factor < 2 indicates pr_trial_func failed to find a factor.
            // This is a low probability, but expected possibility.  We want
            // to retry the trial_func in the next loop iteration using an
            // incremented value of 'i', and keep doing so until we get a good
            // factor, or until 'i' gets to be as large as x (which is so
            // unlikely to happen that it's nearly impossible).
            // Meanwhile, we don't need to do anything here.
        }
        cc = mf.add(cc, unity);
    }
    // We went through every allowed value of i, but didn't find a factor.  This
    // is so unlikely that we could assert it never happens, since it's more
    // likely that a coding error would result in reaching this point than that
    // we could ever legitimately reach this point.
    // Technically... we shouldn't assert that we never reach this point, but it
    // can be useful for testing- if it happens, we might investigate whether it
    // was legit that the algorithm should fail for the current value of x, or
    // whether a coding error was responsible, such as trying to use pollard rho
    // trials on a prime number.
    HPBC_ASSERT2(false);

    // Since we didn't find a factor, fall back to slow trial division.
    // factorize_wheel210 in principle will always completely factor x, though
    // for values of x > (1<<64), it may be too slow to be usable in practice.
    return factorize_wheel210::call(iter, x);
}


}}  // end namespace

#endif
