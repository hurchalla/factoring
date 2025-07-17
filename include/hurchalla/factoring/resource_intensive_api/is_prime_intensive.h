// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FUNCTION_IS_PRIME_INTENSIVE_H_INCLUDED
#define HURCHALLA_FACTORING_FUNCTION_IS_PRIME_INTENSIVE_H_INCLUDED


#if !defined(HURCHALLA_FACTORING_DISALLOW_INLINE_ASM) && \
        !defined(HURCHALLA_ALLOW_INLINE_ASM_ALL)
#  define HURCHALLA_ALLOW_INLINE_ASM_ALL
#endif

#include "hurchalla/factoring/detail/impl_is_prime_intensive.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// is_prime_intensive(T x) returns true if x is prime, and otherwise returns
// false.
//
// T can be any integral type <= 128 bits.
//
// You can optionally specifiy a value for TRIAL_DIVISION_SIZE, which is the
// number of primes that this function will use for trial division primality
// testing, before it begins using more complex algorithms.  Generally speaking,
// it is a good idea to not specify any value at all (just use the default).
// However, if you know beforehand that you will be testing numbers that are
// very likely to be prime, you may be able to improve this function's
// performance by choosing a small value for TRIAL_DIVISION_SIZE.  Or, if you
// know beforehand that you will be testing numbers that are very likely to be
// composite, you may be able to improve performance by choosing a large value
// for TRIAL_DIVISION_SIZE.
template <typename T, unsigned int TRIAL_DIVISION_SIZE =
                 detail::impl_is_prime_intensive::defaultTrialDivisionSize<T>()>
bool is_prime_intensive(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::digits <= 128, "");
    HPBC_PRECONDITION2(x >= 0);

    return detail::impl_is_prime_intensive::call<TRIAL_DIVISION_SIZE>(x);
}


}  // end namespace

#endif
