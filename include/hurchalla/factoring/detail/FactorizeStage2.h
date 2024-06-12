// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FACTORIZE_STAGE2_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_STAGE2_H_INCLUDED


#include "hurchalla/factoring/detail/microecm.h"
#include "hurchalla/factoring/detail/PollardRhoBrentSwitchingTrial.h"
#ifndef HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME
#  define HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME PollardRhoBrentSwitchingTrial
#else
#  include "hurchalla/factoring/detail/PollardRhoBrentTrial.h"
#  include "hurchalla/factoring/detail/PollardRhoBrentTrialParallel.h"
#  include "hurchalla/factoring/detail/experimental/PollardRhoTrial.h"
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
#include <cstdint>


namespace hurchalla { namespace detail {


#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#  pragma warning(disable : 4702)
#endif

// FactorizeStage2 is a functor that checks the size of x, and dispatches to its
// helper function (factorize2) using the fastest Montgomery form for the size
// of x.  The helper function first checks if x is prime, and if so just writes
// x to the output iterator.  If not, it calls a factor finding function (ECM or
// Pollard-Rho, depending on the size of x), and then recursively calls the
// dispatch function again for both the found_factor and x/factor, which ensures
// all factors get written to the output iterator.


template <int EcmMinBits, int MaxBitsX, class T>
class FactorizeStage2 {
    const T always_prime_limit;
    bool expect_arbitrary_size_factors;
    T base_c;               // for Pollard-Rho
    T expected_iterations;  // for Pollard-Rho
    std::uint64_t loc_lcg;  // for ECM

public:
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);

    FactorizeStage2(T alwaysprime_limit, bool expect_arbitrarysize_factors)
        : always_prime_limit(alwaysprime_limit),
          expect_arbitrary_size_factors(expect_arbitrarysize_factors),
          base_c(1), expected_iterations(0), loc_lcg(0)
    {
    }

    template <class OutputIt, class PrimalityFunctor>
    OutputIt operator()(OutputIt iter, const PrimalityFunctor& is_prime_func,
                        T x)
    {
        base_c = 1;
        expected_iterations = 0;
        return dispatch(iter, is_prime_func, x);
    }


private:
    template <class OutputIt, class PrimalityFunctor>
    OutputIt dispatch(OutputIt iter, const PrimalityFunctor& is_prime_func, T x)
    {
        HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
        HPBC_PRECONDITION2(x % 2 == 1); // x odd is required for montgomery form

        using S = sized_uint<HURCHALLA_TARGET_BIT_WIDTH>::type;

        // factor using a native integer type or smaller, if possible
        if (ut_numeric_limits<T>::digits <= ut_numeric_limits<S>::digits ||
                                       x <= ut_numeric_limits<S>::max()) {
            using U = typename std::conditional<
                  (ut_numeric_limits<T>::digits < ut_numeric_limits<S>::digits),
                  T, S>::type;
#if defined(HURCHALLA_FACTORIZE_NEVER_USE_MONTGOMERY_MATH)
            using MF = MontgomeryStandardMathWrapper<U>;
            return factorize2<MF>(iter, is_prime_func, x);
#else
            static_assert(ut_numeric_limits<U>::digits >= 2);
            constexpr U Udiv4 = static_cast<U>(static_cast<U>(1) <<
                                            (ut_numeric_limits<U>::digits - 2));
            constexpr U Udiv2 = static_cast<U>(static_cast<U>(1) <<
                                            (ut_numeric_limits<U>::digits - 1));
            if (x < Udiv4) {
                using MF = MontgomeryQuarter<U>;
                return factorize2<MF>(iter, is_prime_func, x);
            } else if (x < Udiv2) {
                using MF = MontgomeryHalf<U>;
                return factorize2<MF>(iter, is_prime_func, x);
            } else {
                constexpr bool limited_to_halfrange =
                                 (MaxBitsX == ut_numeric_limits<T>::digits - 1);
                // Let's not generate code under conditions where it is
                // impossible for the code to ever be executed.
                // limited_to_halfrange applies to T, so if T and U are the same
                // size and limited_to_halfrange is true, then we know that
                // x < Udiv2.  Thus under those conditions, we would never enter
                // this 'else' clause - we would have entered the clause above
                // instead.
                if constexpr (limited_to_halfrange &&
                 ut_numeric_limits<U>::digits == ut_numeric_limits<T>::digits) {
                    ::hurchalla::unreachable();
                    HPBC_ASSERT2(false); //it ought to be impossible to get here
                    return iter;
                } else {
                    using MF = MontgomeryForm<U>;
                    return factorize2<MF>(iter, is_prime_func, x);
                }
            }
#endif
        }
        // If we reach this point, the following clause will be true.  We
        // explicitly use an 'if' anyway, so that we can be sure the compiler
        // will not generate any code for it when T's bit width is
        // <= HURCHALLA_TARGET_BIT_WIDTH.
        if constexpr (ut_numeric_limits<T>::digits > HURCHALLA_TARGET_BIT_WIDTH) {
#if defined(HURCHALLA_FACTORIZE_NEVER_USE_MONTGOMERY_MATH)
            using MF = MontgomeryStandardMathWrapper<T>;
            return factorize2<MF>(iter, is_prime_func, x);
#else
            static_assert(ut_numeric_limits<T>::digits >= 2);
            T Rdiv4 = static_cast<T>(
                       static_cast<T>(1) << (ut_numeric_limits<T>::digits - 2));
            if (x < Rdiv4) {
                using MF = MontgomeryQuarter<T>;
                return factorize2<MF>(iter, is_prime_func, x);
            } else {
                using MF = MontgomeryForm<T>;
                return factorize2<MF>(iter, is_prime_func, x);
            }
#endif
        } else {
            ::hurchalla::unreachable();
            HPBC_ASSERT2(false);  // it ought to be impossible to reach this
            return iter;
        }
#if defined(__INTEL_COMPILER)
        // for some reason icc doesn't see that all paths return prior to this
        // point, so we'll uselessly return iter to quiet icc's warning/error.
        return iter;
#endif
    }



    template <class MF, class OutputIt, class PrimalityFunctor>
    OutputIt factorize2(OutputIt iter,
                        const PrimalityFunctor& is_prime_func, T x)
    {
        using U = typename MF::IntegerType;
        using P = typename safely_promote_unsigned<T>::type;
        using C = typename MF::CanonicalValue;
        static_assert(ut_numeric_limits<U>::is_integer);
        static_assert(!ut_numeric_limits<U>::is_signed);
        HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations
        HPBC_PRECONDITION2(x % 2 == 1); // odd x is required for montgomery form
        HPBC_PRECONDITION2(x <= ut_numeric_limits<U>::max());

        if (x < always_prime_limit) {
            *iter++ = x;  // x is prime
            return iter;
        }
        MF mf(static_cast<U>(x));

        if (is_prime_func(mf)) {
            *iter++ = x;  // x is prime
            return iter;
        }

        // Try to use microECM if we can
        static_assert(0 <= EcmMinBits);
        if constexpr (EcmMinBits < ut_numeric_limits<T>::digits) {
            int ecm_crossover_bits = (expect_arbitrary_size_factors)
                                     ? EcmMinBits + 6 : EcmMinBits;
            if (ecm_crossover_bits < ut_numeric_limits<T>::digits) {
                T ecm_threshold = static_cast<T>(1) << ecm_crossover_bits;
                if (x >= ecm_threshold) {   // use ecm only if x is 'big'
                    T tmp_factor = micro_ecm::get_ecm_factor(mf,
                                        expect_arbitrary_size_factors, loc_lcg);
                    if (tmp_factor >= 2) {   // we found a factor
                        // Try to factor the factor (fyi, it is usually prime)
                        iter = dispatch(iter, is_prime_func, tmp_factor);
                        // Next try to factor the quotient.
                        // since 1 < tmp_factor < x,
                        // we know 1 < (x/tmp_factor) < x
                        T quotient = static_cast<T>(x/tmp_factor);
                        iter = dispatch(iter, is_prime_func, quotient);
                        return iter;
                    }
                }
            }
        }
        // If we don't use ECM (or if it found no factors), use Pollard-Rho.

        HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME<MF> pr_trial_func;

        // Pollard-Rho factoring:
        // we don't want to use a sequence  x[i+1] = x[i]*x[i] + c  where c == 0
        // or c == -2.  See JM Pollard "A Monte Carlo method for factorization".
        // The following only guarantees we avoid those sequences on the first
        // iteration of the loop, but unless x is very small, it's extremely
        // unlikely base_c would ever grow large enough during loop iterations
        // to reach x-2.  If we do end up with one of those sequences, it's
        // fairly harmless - they just have low success rate, making them
        // inefficient since we may need another loop iteration after it.
        if (static_cast<P>(base_c) >= static_cast<P>(x)-2 || base_c == 0)
            base_c = 1;
        C unity = mf.getUnityValue();

        C cc = mf.getCanonicalValue(mf.convertIn(static_cast<U>(base_c)));
        for (T i = 0; i < x; ++i) {
            HPBC_ASSERT2(expected_iterations <= ut_numeric_limits<U>::max());
            U tmp_expected = static_cast<U>(expected_iterations);
            T tmp_factor = static_cast<T>(pr_trial_func(mf, tmp_expected, cc));
            HPBC_ASSERT2(tmp_expected <= ut_numeric_limits<T>::max());
            expected_iterations = static_cast<T>(tmp_expected);
            if (tmp_factor >= 2) {    // we found a factor.
                // Base_c could overflow, but that's okay.  We'd prefer for
                // efficiency that it didn't, but any T value would be valid.
                base_c = static_cast<T>(base_c + static_cast<P>(i) + 1);
                // Try to factor the factor (fyi, it is usually prime)
                iter = dispatch(iter, is_prime_func, tmp_factor);
                // Next try to factor the quotient.
                // since 1 < tmp_factor < x, we know 1 < (x/tmp_factor) < x
                T quotient = static_cast<T>(x/tmp_factor);
                iter = dispatch(iter, is_prime_func, quotient);
                return iter;
            }
            else {
                // tmp_factor < 2 indicates pr_trial_func failed to find a
                // factor.  This is a low probability, but expected possibility.
                // We want to retry the trial_func in the next loop iteration
                // using an incremented value of 'i', and keep doing so until we
                // get a good factor, or until 'i' gets to be as large as x
                // (which is so unlikely to happen that it's nearly impossible).
                // Meanwhile, we don't need to do anything here.
            }
            cc = mf.add(cc, unity);
        }
        // We went through every allowed value of i, but didn't find a factor.
        // This is so unlikely that we could assert it never happens, since it's
        // more likely that a coding error would result in reaching this point
        // than that we could ever legitimately reach this point.
        // If the assert below gets triggered, we might investigate whether it
        // was legit that the algorithm should fail for the current value of x,
        // or whether a coding error was responsible, such as trying to use
        // pollard rho trials on a prime number.
        HPBC_ASSERT2(false);

        // Since we didn't find a factor, fall back to slow trial division.
        // factorize_wheel210 in principle will always completely factor x,
        // though for values of x > (1<<50), it may effectively be too slow to
        // be usable.
        return factorize_wheel210::call(iter, x);
    }
};

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif


}}  // end namespace

#endif
