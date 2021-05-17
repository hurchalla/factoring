// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_GCD_FUNCTOR_H_INCLUDED
#define HURCHALLA_FACTORING_GCD_FUNCTOR_H_INCLUDED


#include "hurchalla/factoring/greatest_common_divisor.h"

namespace hurchalla { namespace detail {


template <typename T>
struct GcdFunctor {
    HURCHALLA_FORCE_INLINE HURCHALLA_FLATTEN T operator()(T a, T b)
    {
        return greatest_common_divisor(a, b);
    }
};


}}  // end namespace

#endif
