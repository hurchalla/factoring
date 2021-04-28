// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_GREATEST_COMMON_DIVISOR_H_INCLUDED
#define HURCHALLA_FACTORING_GREATEST_COMMON_DIVISOR_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


// For details, see https://en.wikipedia.org/wiki/Greatest_common_divisor
//
// Note that it is most efficient to provide inputs a<=b (otherwise the first
// loop of the algorithm effectively performs swap(a,b) via a % operation).
//
template <typename T>
T greatest_common_divisor(T a, T b)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(a > 0 || b > 0);

    while (a != 0) {
        T tmp = a;
        a = static_cast<T>(b % a);
        b = tmp;
    }

    HPBC_POSTCONDITION2(b > 0);
    return b;
}


}  // end namespace

#endif
