// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/impl_is_prime.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


template <typename T>
bool is_prime(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    return detail::impl_is_prime(x);
}


}  // end namespace

#endif
