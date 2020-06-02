
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
    HPBC_PRECONDITION2(x >= 0);

    if (x < 2)
        return false;
    if (x%2 == 0)
        return (x==2);

    T s = integer_sqrt(x);
    for (T f=3; f<=s; f+=2) {  // skip even factors- we already checked them
        if (x%f == 0) {
            HPBC_ASSERT2(f < x);
            return false;
        }
    }
    return true;
}


}}  // end namespace

#endif
