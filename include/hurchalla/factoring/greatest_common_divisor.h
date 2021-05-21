// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_GREATEST_COMMON_DIVISOR_H_INCLUDED
#define HURCHALLA_FACTORING_GREATEST_COMMON_DIVISOR_H_INCLUDED


#include "hurchalla/factoring/detail/impl_greatest_common_divisor.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// For info on GCD see https://en.wikipedia.org/wiki/Greatest_common_divisor
template <typename T>
T greatest_common_divisor(T a, T b)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(a > 0 || b > 0);

    T gcd = hurchalla::detail::impl_greatest_common_divisor(a, b);

    HPBC_POSTCONDITION2(gcd > 0);
    return gcd;
}


}  // end namespace

#endif
