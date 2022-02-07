// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FACTORIZE_BRUTEFORCE_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_BRUTEFORCE_H_INCLUDED


#include "integer_sqrt.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <vector>
#include <limits>

namespace hurchalla {


// Postconditions:
// Returns a vector that contains all the prime factors of x.  If x is prime
//   then the vector will contain only the single element x.  The vector will
//   never be empty.
//
// factorize_bruteforce() guarantees it will trial all potential prime factors
// for x.  Note that a large input value for x can pose an intractable problem
// for this algorithm, effectively causing the function to never complete.


template <typename T>
std::vector<T> factorize_bruteforce(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations

    std::vector<T> factors;
    T q = x;
    HPBC_ASSERT2(q > 1);

    while (q % 2 == 0) {
        factors.push_back(2);
        q = static_cast<T>(q >> 1);
        if (q == 1)
            return factors;
    }
    HPBC_ASSERT2(q > 1);

    T s = integer_sqrt(q);  // no factor can be > sqrt(q)
    // skip all even trial factors- we've already tested for divisibility by 2
    for (T f = 3; f <= s; f = static_cast<T>(f + 2)) {
        HPBC_ASSERT2(q > 1);
        while (q % f == 0) {
            factors.push_back(f);
            q = static_cast<T>(q/f);
            if (q == 1)
                return factors;
            s = integer_sqrt(q);
        }
    }
    HPBC_ASSERT2(q > 1);
    // since no factors exist < sqrt(q), q must be prime
    factors.push_back(q);
    return factors;
}


}  // end namespace

#endif
