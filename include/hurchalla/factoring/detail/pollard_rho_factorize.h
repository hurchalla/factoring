// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_FACTORIZE_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_FACTORIZE_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/detail/factorize_wheel210.h"
#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// defined lower in this file
template <template<class,int> class TTD, template<class> class Functor,
          class OutputIt, typename T>
OutputIt pollard_rho_factorize(OutputIt iter, T x, T threshold_always_prime = 0,
                               T base_c = 1, T* pIterations_performed = nullptr);


namespace prf_detail {

template <template<class,int> class TTD, template<class> class Functor,
          class MF, class OutputIt, typename T>
OutputIt pr_factorize(OutputIt iter, T x, T threshold_always_prime, T base_c,
                                                       T* pIterations_performed)
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

    if (pIterations_performed) *pIterations_performed = 0;

    MF mf(static_cast<S>(x));

    if (x < threshold_always_prime || is_prime_miller_rabin(mf)) {
        *iter++ = x;  // x is prime
        return iter;
    }
    Functor<MF> pr_trial_func;

    C cc = (base_c < x) ? mf.getCanonicalValue(
                     mf.convertIn(static_cast<S>(base_c))) : mf.getUnityValue();
    for (T c = base_c; c < x; ++c) {
        S iterations;
        T tmp_factor = static_cast<T>(pr_trial_func(mf, cc, &iterations));
        if (pIterations_performed)
            *pIterations_performed =
                            static_cast<T>(*pIterations_performed + iterations);

        if (tmp_factor >= 2) // we found a good factor (we don't know if prime)
        {
            // Try to factor the factor (it may or may not be prime)
            T sub_iterations;
            iter = pollard_rho_factorize<TTD, Functor>(iter, tmp_factor,
                  threshold_always_prime, static_cast<T>(c+1), &sub_iterations);
            if (pIterations_performed)
                *pIterations_performed = static_cast<T>(*pIterations_performed
                                                              + sub_iterations);

            // Next try to factor the quotient.
            // since 1 < tmp_factor < x, we know 1 < (x/tmp_factor) < x
            T quotient = static_cast<T>(x/tmp_factor);  
            iter = pollard_rho_factorize<TTD, Functor>(iter, quotient,
                  threshold_always_prime, static_cast<T>(c+1), &sub_iterations);
            if (pIterations_performed)
                *pIterations_performed = static_cast<T>(*pIterations_performed
                                                              + sub_iterations);
            return iter;
        }
        else
        {
            // tmp_factor < 2 indicates the trial_func failed to find a factor.
            // This is a low probability, but expected possibility.  We want
            // to retry the trial_func in the next loop iteration using an
            // incremented value of 'c', and keep doing so until we get a good
            // factor, or until 'c' gets to be as large as x (which is so
            // unlikely to happen that it's nearly impossible).
            // Meanwhile, we don't need to do anything here.
        }
        cc = mf.getCanonicalValue(mf.add(cc, mf.getUnityValue()));
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

    // set the number of iterations to max possible, since Pollard Rho failed to
    // find any factor at all.
    if (pIterations_performed)
        *pIterations_performed = ut_numeric_limits<T>::max();

    // Since we didn't find a factor, fall back to slow trial division
    T q;
    iter = factorize_wheel210<TTD>(iter, q, x);
    // factorize_wheel210 should always completely factor x (and set q = 1)
    HPBC_ASSERT2(q == 1);
    return iter;
}

} // end namespace prf_detail



// the following macro is intended just for testing purposes
//#define HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH 1



// Dispatch function to the fastest pr_factorize template function instantiation
// available for x (parameterized on the fastest MontgomeryForm for x).
// The template-template param TTD should be either TrialDivisionWarren or
// TrialDivisionMayer.
template <template<class,int> class TTD, template<class> class Functor,
          class OutputIt, typename T>
OutputIt pollard_rho_factorize(OutputIt iter, T x, T threshold_always_prime,
                                             T base_c, T* pIterations_performed)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
    HPBC_PRECONDITION2(x % 2 == 1);  // x odd is required for montgomery form

    using S = sized_uint<HURCHALLA_TARGET_BIT_WIDTH>::type;
    // factor using a native integer type whenever possible
    if (x <= ut_numeric_limits<S>::max()) {
#if defined(HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH)
        using MF = MontgomeryStandardMathWrapper<S>;
        return prf_detail::pr_factorize<TTD, Functor, MF>(iter, x,
                         threshold_always_prime, base_c, pIterations_performed);
#else
        S Sdiv2 = static_cast<S>(static_cast<S>(1) <<
                                              (HURCHALLA_TARGET_BIT_WIDTH - 1));
        S Sdiv4 = static_cast<S>(Sdiv2 / 2);
        if (x < Sdiv4) {
            using MF = MontgomeryQuarter<S>;
            return prf_detail::pr_factorize<TTD, Functor, MF>(iter, x,
                         threshold_always_prime, base_c, pIterations_performed);
        } else {
            using MF = MontgomeryFull<S>;
            return prf_detail::pr_factorize<TTD, Functor, MF>(iter, x,
                         threshold_always_prime, base_c, pIterations_performed);
        }
#endif
    }
    // If we reach this point, the following clause will be true.  We explicitly
    // use an 'if' anyway, so that we can be sure the compiler will not generate
    // any code for it when T digits <= HURCHALLA_TARGET_BIT_WIDTH.
    if (ut_numeric_limits<T>::digits > HURCHALLA_TARGET_BIT_WIDTH) {
#if defined(HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH)
        using MF = MontgomeryStandardMathWrapper<T>;
        return prf_detail::pr_factorize<TTD, Functor, MF>(iter, x,
                         threshold_always_prime, base_c, pIterations_performed);
#else
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                            (ut_numeric_limits<T>::digits - 1));
        T Rdiv4 = static_cast<T>(Rdiv2 / 2);
        if (x < Rdiv4) {
            using MF = MontgomeryQuarter<T>;
            return prf_detail::pr_factorize<TTD, Functor, MF>(iter, x,
                         threshold_always_prime, base_c, pIterations_performed);
        } else {
            using MF = MontgomeryFull<T>;
            return prf_detail::pr_factorize<TTD, Functor, MF>(iter, x,
                         threshold_always_prime, base_c, pIterations_performed);
        }
#endif
    } else {
        HPBC_ASSERT2(false);  // it ought to be impossible to reach this code.
        return iter;
    }
}


}}  // end namespace

#endif
