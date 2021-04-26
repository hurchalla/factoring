// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_INTEGER_SQRT_H_INCLUDED
#define HURCHALLA_FACTORING_INTEGER_SQRT_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>
#include <limits>
#include <cstdint>

namespace hurchalla {


namespace detail_integer_sqrt {
// this stops the compiler from complaining about "shift count >= width of type"
// in cases where the shift wouldn't have actually been executed
template <int shift, typename T>
HURCHALLA_FORCE_INLINE typename std::enable_if<!(shift <
                                         ut_numeric_limits<T>::digits), T>::type
safe_right_shift(T)
{
    return 0;
}
template <int shift, typename T>
HURCHALLA_FORCE_INLINE  typename std::enable_if<(shift <
                                         ut_numeric_limits<T>::digits), T>::type
safe_right_shift(T x1)
{
    return static_cast<T>(x1 >> shift);
}
}


// The integer sqrt algorithm here is based on Newton's method.  For a detailed
// description and proof of correctness of the algorithm, see "Hacker's Delight"
// 2nd edition, chapter 11 section 1, by Henry S. Warren.  This function is
// adapted from Warren's Figure 11-1 function isqrt(unsigned int).
template <typename T>
T integer_sqrt(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::digits <= 128, "");
    HPBC_PRECONDITION2(x >= 0);

    if (x <= 1)
        return x;

    using namespace detail_integer_sqrt;
    // set the first guess g0 equal to the least power of 2 that is >= sqrt(x)
    int s = 1;
    T x1 = static_cast<T>(x - 1);
    if (x1 > UINT64_C(18446744073709551615))
                                   { s += 32; x1 = safe_right_shift<64>(x1); }
    // note x1 now fits in a uint64_t, and after next line it fits uint32_t, etc
    if (x1 > UINT32_C(4294967295)) { s += 16; x1 = safe_right_shift<32>(x1); }
    if (x1 > UINT16_C(65535))      { s +=  8; x1 = safe_right_shift<16>(x1); }
    if (x1 > 255)                  { s +=  4; x1 = safe_right_shift<8>(x1);  }
    if (x1 > 15)                   { s +=  2; x1 = safe_right_shift<4>(x1);  }
    if (x1 > 3)                    { s +=  1; }
    T g0 = static_cast<T>(static_cast<T>(1) << s);

    T g1 = static_cast<T>((g0 + (x >> s)) >> 1);  // same as g1 = (g0 + x/g0)/2
    while (g1 < g0) {            // Loop until approximations no longer decrease
        g0 = g1;
        g1 = static_cast<T>((g0 + (x/g0)) >> 1);
    }
    return g0;
}


}  // end namespace

#endif
