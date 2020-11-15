
#ifndef HURCHALLA_FACTORING_GREATEST_COMMON_DIVISOR_H_INCLUDED
#define HURCHALLA_FACTORING_GREATEST_COMMON_DIVISOR_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace factoring {


// For details, see https://en.wikipedia.org/wiki/Greatest_common_divisor
//
// Note that it is most efficient to provide inputs a<=b (otherwise the first
// loop of the algorithm effectively performs swap(a,b) via a % operation).
//
template <typename T>
T greatest_common_divisor(T a, T b)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(a > 0 || b > 0);

    while (a != 0) {
        T tmp = a;
        a = static_cast<T>(b % a);
        b = tmp;
    }

    HPBC_POSTCONDITION2(b > 0);
    return b;
}


}}  // end namespace

#endif
