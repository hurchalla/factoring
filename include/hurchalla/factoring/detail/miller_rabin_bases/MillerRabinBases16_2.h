// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES16_2_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES16_2_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// You can use these 2 bases to determine the primality of any 16 bit unsigned
// int number, via 2 base miller-rabin primality testing.  The 2 bases are
// constants (no hash tables are used to set them).
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<16, 2, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
    // This function returns 2 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of any 16
    // bit unsigned integer.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint8_t, 2> get(std::uint16_t)
    {
        // I generated/verified these bases.  See README.TXT
        const std::array<std::uint8_t, 2> bases = { 2, 3 };
        return bases;
    }
};


}}  // end namespace

#endif
