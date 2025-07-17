// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IS_PRIME_ULTIMATE_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_ULTIMATE_H_INCLUDED


#if !defined(HURCHALLA_FACTORING_DISALLOW_INLINE_ASM) && \
        !defined(HURCHALLA_ALLOW_INLINE_ASM_ALL)
#  define HURCHALLA_ALLOW_INLINE_ASM_ALL
#endif

#include "hurchalla/factoring/resource_intensive_api/IsPrimeTable.h"
#include "hurchalla/factoring/detail/impl_is_prime_intensive.h"
#include "hurchalla/factoring/resource_intensive_api/is_prime_intensive.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// is_prime_ultimate(T x) returns true if x is prime, and otherwise returns
// false.
//
// T can be any integral type <= 128 bits.
//
// This function is intended for use when you plan to repeatedly and intensively
// perform primality testing.  It requires you to pass a reference to an
// IsPrimeTable<uint32_t> object, which may take a few seconds for you to
// construct prior to your calls to this function.  Excluding the time it takes
// for you to construct the IsPrimeTable<uint32_t> object, this function
// provides the fastest primality testing available in this library.  It uses a
// mix of IsPrimeTable (the fastest test of values that can fit in 32 bits) and
// is_prime_intensive (the fastest test for values above 32 bits).  Note that
// IsPrimeTable uses 256MB of memory, and so this function indirectly uses about
// 256MB - thus it has a large impact on memory (and CPU cache) use.  In
// contrast, the normal API function is_prime() uses an amount of memory (less
// than 1 KB) that is negligible for most purposes.
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
bool is_prime_ultimate(T x, const IsPrimeTable<uint32_t>& table)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::digits <= 128, "");
    HPBC_PRECONDITION2(x >= 0);

    using U = typename extensible_make_unsigned<T>::type;
    if (static_cast<U>(x) <= ut_numeric_limits<uint32_t>::max())
        return table(static_cast<uint32_t>(x));
    else
        return is_prime_intensive<T, TRIAL_DIVISION_SIZE>(x);
}


}  // end namespace

#endif
