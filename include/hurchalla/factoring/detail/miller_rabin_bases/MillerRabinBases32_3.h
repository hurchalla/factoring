// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES32_3_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES32_3_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// You can use these 3 bases to determine the primality of any 32 bit unsigned
// int number, via 3 base miller-rabin primality testing.  The 3 bases are
// constants (no hash tables are used to set them).
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<32, 3, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
    // This function returns 3 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of any 32
    // bit unsigned integer.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint8_t, 3> get(std::uint32_t)
    {
        // These 3 bases were discovered by Gerhard Jaeschke, according to
        // http://miller-rabin.appspot.com
        // I verified they are correct.  See README.TXT.
        const std::array<std::uint8_t, 3> bases = { 2, 7, 61 };
        return bases;
    }
};


}}  // end namespace

#endif
