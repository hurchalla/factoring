// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_FACTORIZE_POLLARD_RHO_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_POLLARD_RHO_H_INCLUDED


#include "hurchalla/factoring/detail/PollardRhoTrial.h"
#include "hurchalla/factoring/detail/PollardRhoBrentTrial.h"
#include "hurchalla/factoring/detail/factorize_wheel210.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unreachable.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


#ifndef HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME
#if 1
#  define HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME PollardRhoBrentTrial
#else
#  define HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME PollardRhoTrial
#endif
#endif


// defined lower in this file
template <class PrimalityFunctor, class OutputIt, typename T>
OutputIt factorize_pollard_rho(OutputIt iter, T x,
                               const PrimalityFunctor& is_prime_pr,
                               T threshold_always_prime = 0, T base_c = 1);


namespace prf_detail {

template <class MF, class PrimalityFunctor, class OutputIt, typename T>
OutputIt factorize_pr(OutputIt iter, T x, const PrimalityFunctor& is_prime_pr,
                                             T threshold_always_prime, T base_c)
{
    using S = typename MF::T_type;
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

    HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME<MF> pr_trial_func;
    if (base_c >= x)
        base_c = 1;
    C unity = mf.getUnityValue();
    C cc = mf.getCanonicalValue(mf.convertIn(static_cast<S>(base_c)));
    for (T c = base_c; c < x; ++c,cc = mf.getCanonicalValue(mf.add(cc,unity))) {
        T tmp_factor = static_cast<T>(pr_trial_func(mf, cc));
        if (tmp_factor >= 2) {    // we found a good factor (maybe prime)
            // Try to factor the factor (it may or may not be prime)
            iter = factorize_pollard_rho(iter, tmp_factor, is_prime_pr,
                                   threshold_always_prime, static_cast<T>(c+1));
            // Next try to factor the quotient.
            // since 1 < tmp_factor < x, we know 1 < (x/tmp_factor) < x
            T quotient = static_cast<T>(x/tmp_factor);  
            iter = factorize_pollard_rho(iter, quotient, is_prime_pr,
                                   threshold_always_prime, static_cast<T>(c+1));
            return iter;
        }
        else {
            // tmp_factor < 2 indicates pr_trial_func failed to find a factor.
            // This is a low probability, but expected possibility.  We want
            // to retry the trial_func in the next loop iteration using an
            // incremented value of 'c', and keep doing so until we get a good
            // factor, or until 'c' gets to be as large as x (which is so
            // unlikely to happen that it's nearly impossible).
            // Meanwhile, we don't need to do anything here.
        }
    }
    // We went through every allowed value of c, but didn't find a factor.  This
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
    return factorize_wheel210(iter, x);
}

} // end namespace prf_detail



// the following macro is intended only for testing purposes
//#define HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH 1


// Dispatch function to the fastest factorize_pr template function instantiation
// available for x (parameterized on the fastest MontgomeryForm for x).
template <class PrimalityFunctor, class OutputIt, typename T>
OutputIt factorize_pollard_rho(OutputIt iter, T x,
                               const PrimalityFunctor& is_prime_pr,
                               T threshold_always_prime, T base_c)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    HPBC_PRECONDITION2(x % 2 == 1);  // x odd is required for montgomery form

    using S = sized_uint<HURCHALLA_TARGET_BIT_WIDTH>::type;

    // factor using a native integer type if possible
    if (ut_numeric_limits<T>::digits <= ut_numeric_limits<S>::digits ||
        x <= ut_numeric_limits<S>::max()) {
        using U = typename std::conditional<
                    ut_numeric_limits<T>::digits < ut_numeric_limits<S>::digits,
                    T, S>::type;
#if defined(HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH)
        using MF = MontgomeryStandardMathWrapper<U>;
        return prf_detail::factorize_pr<MF>(iter, x, is_prime_pr,
                                                threshold_always_prime, base_c);
#else
        U Udiv4 = static_cast<U>(static_cast<U>(1) <<
                                            (ut_numeric_limits<U>::digits - 2));
        if (x < Udiv4) {
            using MF = MontgomeryQuarter<U>;
            return prf_detail::factorize_pr<MF>(iter, x, is_prime_pr,
                                                threshold_always_prime, base_c);
        } else {
            using MF = MontgomeryForm<U>;
            return prf_detail::factorize_pr<MF>(iter, x, is_prime_pr,
                                                threshold_always_prime, base_c);
        }
#endif
    }
    // If we reach this point, the following clause will be true.  We explicitly
    // use an 'if' anyway, so that we can be sure the compiler will not generate
    // any code for it when T's bit width <= HURCHALLA_TARGET_BIT_WIDTH.
#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif
    if constexpr (ut_numeric_limits<T>::digits > HURCHALLA_TARGET_BIT_WIDTH) {
#if defined(HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH)
        using MF = MontgomeryStandardMathWrapper<T>;
        return prf_detail::factorize_pr<MF>(iter, x, is_prime_pr,
                                                threshold_always_prime, base_c);
#else
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 2));
        if (x < Rdiv4) {
            using MF = MontgomeryQuarter<T>;
            return prf_detail::factorize_pr<MF>(iter, x, is_prime_pr,
                                                threshold_always_prime, base_c);
        } else {
            using MF = MontgomeryForm<T>;
            return prf_detail::factorize_pr<MF>(iter, x, is_prime_pr,
                                                threshold_always_prime, base_c);
        }
#endif
    } else {
        unreachable();
        HPBC_ASSERT2(false);  // it ought to be impossible to reach this code.
        return iter;
    }
#if defined(_MSC_VER)
#  pragma warning(pop)
#endif
#if defined(__INTEL_COMPILER)
    // for some reason icc doesn't see that all paths return prior to this point
    // so we'll uselessly return iter to quiet icc's warning/error
    return iter;
#endif
}


}}  // end namespace

#endif
