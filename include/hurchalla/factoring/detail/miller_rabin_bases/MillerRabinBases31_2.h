// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES31_2_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES31_2_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 12 byte hash table lets you determine the primality of any unsigned int
// number less than (1<<31), via miller-rabin primality testing using 2 bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<31, 2, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 32 bit number being tested for primality.
    // As a precondition, 'num' must be less than (1<<31).
    // This function returns 2 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 2> get(std::uint32_t num)
    {
        HPBC_PRECONDITION2(num < (static_cast<std::uint32_t>(1) << 31));
        // I generated/verified the hash table and bases.  See README.TXT
        std::array<std::uint16_t, 2> bases;
        bases[0] = 41334;
        const std::array<std::uint16_t, 6> table =
                         { 554, 29078, 61981, 25681, 44173, 28415 };
        std::uint32_t hash_bucket = ((num & 4095) * 3) >> 11;
        bases[1] = table[hash_bucket];
        return bases;
    }
};


}}  // end namespace

#endif
