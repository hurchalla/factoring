// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_IS_PRIME_BRUTEFORCE_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_BRUTEFORCE_H_INCLUDED


#include "integer_sqrt.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace factoring {


template <typename T>
bool is_prime_bruteforce(T x)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!std::numeric_limits<T>::is_signed, "");

    if (x < 2)
        return false;
    if (x%2 == 0)
        return (x==2);

    // normally I'd use s = integer_sqrt(x), and test f <= s in the loop, but
    // gcc complains about a missed loop optimization doing it that way because
    // it can't prove overflow of f is impossible (even though it indeed would
    // be impossible).  The following is equivalent and it silences gcc.
    T s = static_cast<T>(integer_sqrt(x) + 1);

    // skip even factors- we already checked them
    for (T f = 3; f < s; f = static_cast<T>(f + 2)) {
        if (x%f == 0) {
            HPBC_ASSERT2(f < x);
            return false;
        }
    }
    return true;
}


}}  // end namespace

#endif
