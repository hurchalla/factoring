// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES32_2_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES32_2_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 16 byte hash table lets you determine the primality of any 32 bit
// unsigned int number, via miller-rabin primality testing using 2 bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<32, 2, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 32 bit number being tested for primality.
    // This function returns 2 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 2> get(std::uint32_t num)
    {
        // I generated/verified the hash table and bases.  See README.TXT
        std::array<std::uint16_t, 2> bases;
        bases[0] = 15;
        const std::array<std::uint16_t, 8> table =
                         { 21232, 56020, 5049, 7595, 13468, 54516, 5681, 4473 };
        std::uint8_t hash_bucket = static_cast<std::uint8_t>((num >> 10) & 7u);
        bases[1] = table[hash_bucket];
        return bases;
    }
};


}}  // end namespace

#endif
