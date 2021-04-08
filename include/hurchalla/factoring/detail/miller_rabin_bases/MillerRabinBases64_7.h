// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES64_7_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES64_7_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// You can use these 7 bases to determine the primality of any 64 bit unsigned
// int number, via 7 base miller-rabin primality testing.  The 7 bases are
// constants (no hash tables are used to set them).
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<64, 7, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
    // This function returns 7 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of any 64
    // bit unsigned integer.
    static std::array<std::uint32_t, 7> get(std::uint64_t)
    {
        // These 7 bases were discovered by Jim Sinclair, according to
        // http://miller-rabin.appspot.com
        // I verified they are correct.  See README.TXT.
        const std::array<std::uint32_t, 7> bases = { UINT32_C(2), UINT32_C(325),
                              UINT32_C(9375), UINT32_C(28178), UINT32_C(450775),
                              UINT32_C(9780504), UINT32_C(1795265022) };
        return bases;
    }
};


}}  // end namespace

#endif
