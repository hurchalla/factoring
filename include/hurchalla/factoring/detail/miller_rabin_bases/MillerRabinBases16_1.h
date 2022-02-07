// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES16_1_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES16_1_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 8 byte hash table lets you determine the primality of any 16 bit
// unsigned int number, via miller-rabin primality testing using 1 single base.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<16, 1, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
    // 'num' is the unsigned 16 bit number being tested for primality.
    // This function returns a single base that can be used by miller-rabin
    // testing to correctly determine (non-probabilistically) the primality of
    // num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 1> get(std::uint16_t num)
    {
        // I generated/verified the hash table for the base.  See README.TXT
        std::array<std::uint16_t, 1> bases;
        std::uint16_t hash_bucket = static_cast<std::uint16_t>((num >> 1) & 3);
        const std::array<std::uint16_t, 4> table = { 2779, 2, 19203, 3027 };
        bases[0] = table[hash_bucket];
        return bases;
    }
};


}}  // end namespace

#endif
