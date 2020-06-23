// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/impl_is_prime.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace factoring {


template <typename T>
bool is_prime(T x)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION2(x >= 0);

    return impl_is_prime(x);
}


}}  // end namespace

#endif
