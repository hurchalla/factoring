// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_GET_SINGLE_FACTOR_H_INCLUDED
#define HURCHALLA_FACTORING_GET_SINGLE_FACTOR_H_INCLUDED


#include "hurchalla/factoring/is_prime.h"
#include "hurchalla/factoring/detail/experimental/microecm_cpp.h"
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
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {

namespace detail {

// For internal use, do not call directly.
struct PollardRhoGetFactor {
    template <class MF>
    typename MF::IntegerType operator()(const MF& mf) const
    {
        using T = typename MF::IntegerType;
        using C = typename MF::CanonicalValue;
        T x = mf.getModulus();
        HPBC_PRECONDITION(!is_prime(x));
        HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME<MF> pr_trial;
        T expected_iterations = 0;
        C unity = mf.getUnityValue();
        C cc = unity;
        for (T i = 0; i < x; ++i) {
            T result = static_cast<T>(pr_trial(mf, expected_iterations, cc));
            if (result >= 2)
                return result;
            cc = mf.add(cc, unity);
        }
        return 0;
    }
};
// For internal use, do not call directly.
struct EcmGetFactor {
    template <class MF>
    typename MF::IntegerType operator()(const MF& mf) const
    {
        HPBC_PRECONDITION(!is_prime(mf.getModulus()));
        uint64_t loc_lcg = 0;
        return micro_ecm::get_ecm_factor(mf, loc_lcg);
    }
};
// For internal use, do not call directly.
template <typename T, class TrialFunctor>
T internal_get_single_factor(T x, const TrialFunctor& trial_functor)
{
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    if (x % 2 == 0)
        return 2;
    // Montgomery arithmetic needs an odd x, which we have.

    constexpr T Udiv4 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 2);
    constexpr T Udiv2 = static_cast<T>(1) << (ut_numeric_limits<T>::digits - 1);
    if (x < Udiv4) {
        using MF = MontgomeryQuarter<T>;
        return trial_functor(MF(x));
    }
    if constexpr (ut_numeric_limits<T>::digits <= 64) {
        if (x < Udiv2) {
            using MF = MontgomeryHalf<T>;
            return trial_functor(MF(x));
        }
    }
    using MF = MontgomeryForm<T>;
    return trial_functor(MF(x));
}

} // end hurchalla::detail namespace


// ------- API --------

// Returns a single factor of x, using the ECM algorithm.
// Has the precondition that x must be composite.
// Factoring is by ECM only; there is no trial division stage, which makes it
// relatively slow at finding very small factors.
template <typename T>
T get_single_factor_ecm(T x)
{
    static_assert(ut_numeric_limits<T>::digits <= 128);  //factor up to 128 bits
    HPBC_PRECONDITION(!is_prime(x));
    return detail::internal_get_single_factor(x, detail::EcmGetFactor());
}

// Returns a single factor of x, using the Pollard-Rho Brent algorithm.
// Has the precondition that x must be composite.
// Factoring is by Pollard-Rho only; there is no trial division stage, which
// makes it relatively slow at finding very small factors.
template <typename T>
T get_single_factor_pollard_rho(T x)
{
    static_assert(ut_numeric_limits<T>::digits <= 128);  //factor up to 128 bits
    HPBC_PRECONDITION(!is_prime(x));
    return detail::internal_get_single_factor(x, detail::PollardRhoGetFactor());
}


} // end hurchalla namespace

#endif
