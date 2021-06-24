// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IS_PRIME_BRUTEFORCE_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_BRUTEFORCE_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla {


template <typename T>
constexpr bool is_prime_bruteforce(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    if (x < 2)
        return false;
    if (x%2 == 0)
        return (x==2);

    T sqrtT = static_cast<T>(1) << (ut_numeric_limits<T>::digits / 2);

    // skip even factors- we already checked them
    for (T f = 3; (f < sqrtT) && (f*f <= x); f = static_cast<T>(f + 2)) {
        if (x%f == 0) {
            HPBC_CONSTEXPR_ASSERT(f < x);
            return false;
        }
    }
    return true;
}


}  // end namespace

#endif
