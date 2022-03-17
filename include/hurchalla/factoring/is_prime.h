// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/impl_is_prime.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// A general purpose function for testing primality of an integer x.
// (As an alternative, if you are doing intensive repeated primality testing,
// you might prefer /resource_intensive_api/IsPrimeIntensive.h)

// T can be any unsigned integral type <= 128 bits.
//
// Returns true if x is prime, and false if is x is not prime.
template <typename T>
bool is_prime(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits <= 128, "");

    return ::hurchalla::detail::impl_is_prime::call(x);
}


}  // end namespace

#endif
