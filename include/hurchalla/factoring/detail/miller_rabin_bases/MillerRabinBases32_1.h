// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES32_1_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES32_1_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 512 byte hash table lets you determine the primality of any 32 bit
// unsigned int number, via miller-rabin primality testing using 1 single base.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<32, 1, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 32 bit number being tested for primality.
    // This function returns a single base that can be used by miller-rabin
    // testing to correctly determine (non-probabilistically) the primality of
    // num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 1> get(std::uint32_t num)
    {
        // I generated/verified the hash table for the base.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 1> bases;
        uint32_t hash_bucket = static_cast<uint32_t>(num ^ (num << 22)) >> 24;
        bases[0] = table[hash_bucket];
        return bases;
    }
private:
    static constexpr std::size_t SIZE = 256;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE] = {
#else
    static constexpr std::uint16_t table[] = {
#endif
        2678,
        158,
        42718,
        2707,
        1460,
        1882,
        22456,
        9875,
        4067,
        731,
        25935,
        7163,
        27002,
        6472,
        2787,
        6026,
        3263,
        10124,
        12084,
        7424,
        52928,
        15830,
        40810,
        2332,
        2968,
        1528,
        7941,
        4004,
        18135,
        30084,
        10413,
        1273,
        491,
        10196,
        18536,
        16743,
        5143,
        1478,
        28,
        18695,
        5074,
        13106,
        4310,
        1526,
        1976,
        16494,
        7084,
        16142,
        1078,
        26784,
        12909,
        8697,
        2244,
        244,
        19851,
        7116,
        23346,
        11009,
        40992,
        5338,
        325,
        27860,
        17208,
        28700,
        11620,
        19706,
        5819,
        630,
        10693,
        1548,
        33843,
        594,
        739,
        8715,
        435,
        7942,
        1150,
        3828,
        2320,
        337,
        567,
        4458,
        1194,
        3232,
        510,
        34937,
        430,
        8328,
        11504,
        1856,
        55944,
        6884,
        9088,
        9402,
        3540,
        830,
        240,
        41356,
        32698,
        11850,
        8410,
        11500,
        1826,
        4357,
        13672,
        1469,
        2102,
        21769,
        4797,
        71,
        1698,
        750,
        56262,
        2998,
        2191,
        10498,
        6681,
        1899,
        9211,
        1272,
        9506,
        16800,
        688,
        1130,
        29096,
        10251,
        4798,
        354,
        24800,
        1956,
        21844,
        3988,
        7492,
        42422,
        11830,
        26248,
        35170,
        1268,
        814,
        25230,
        5631,
        6,
        2560,
        852,
        2623,
        6053,
        7306,
        1975,
        5772,
        6828,
        916,
        7922,
        407,
        48408,
        7610,
        154,
        9217,
        396,
        3553,
        1259,
        19605,
        844,
        9907,
        15598,
        6590,
        2559,
        17430,
        5509,
        4009,
        14312,
        26544,
        27778,
        9689,
        6312,
        12135,
        11642,
        3186,
        7695,
        5125,
        3706,
        3538,
        1818,
        7232,
        26620,
        3676,
        2208,
        8120,
        9979,
        5478,
        7652,
        23028,
        26,
        15868,
        982,
        53285,
        7644,
        26922,
        3592,
        17938,
        9914,
        11793,
        1178,
        34751,
        2128,
        555,
        597,
        4414,
        23707,
        7364,
        2752,
        1834,
        818,
        7930,
        6505,
        5447,
        3588,
        4896,
        3196,
        477,
        8616,
        420,
        2723,
        6682,
        2435,
        22074,
        5500,
        50595,
        3094,
        25979,
        1624,
        3236,
        6179,
        6815,
        340,
        1834,
        705,
        4888,
        1262,
        10380,
        1746,
        23688,
        9763,
        9122,
        4231,
        2634,
        4028,
        1392,
        24759,
        8253,
        20811,
        1365,
        1975,
        610,
        6187,
        224,
        2280
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
};
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<32, 1, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<32, 1, DUMMY>::table[];


}}  // end namespace

#endif
