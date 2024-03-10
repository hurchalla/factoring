// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES31_1_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES31_1_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 384 byte hash table lets you determine the primality of any unsigned int
// number less than (1<<31), via miller-rabin primality testing using 1 single
// base.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<31, 1, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 32 bit number being tested for primality.
    // As a precondition, 'num' must be less than (1<<31).
    // This function returns a single base that can be used by miller-rabin
    // testing to correctly determine (non-probabilistically) the primality of
    // num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 1> get(std::uint32_t num)
    {
        HPBC_PRECONDITION2(num < (static_cast<std::uint32_t>(1) << 31));
        // I generated/verified the hash table for the base.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 1> bases;
        // 46073 is significant only in that it happened to be a multiplier that
        // let me successfully generate a working table.  It took many tries
        // with different multipliers to get one that worked for all 192 bins.
        std::uint16_t hashednum = static_cast<std::uint16_t>(num *
                                             static_cast<std::uint16_t>(46073));
        uint32_t hash_bucket = (static_cast<uint32_t>(hashednum) * 3u) >> 10;
        bases[0] = table[hash_bucket];
        return bases;
    }
private:
    static constexpr std::size_t SIZE = 192;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE] = {
#else
    static constexpr std::uint16_t table[] = {
#endif
        11235,
        3035,
        78,
        2594,
        21908,
        37794,
        16254,
        59146,
        383,
        8636,
        7590,
        3465,
        4959,
        12434,
        11549,
        7288,
        21493,
        11667,
        1092,
        21379,
        6790,
        13793,
        4218,
        33338,
        271,
        1666,
        15106,
        3554,
        5297,
        1314,
        467,
        744,
        1975,
        1508,
        33656,
        53666,
        18056,
        24917,
        45789,
        15974,
        5392,
        22184,
        1791,
        594,
        14438,
        4613,
        2598,
        1066,
        3230,
        5750,
        10504,
        9808,
        6787,
        16750,
        2093,
        56788,
        850,
        4119,
        10178,
        6705,
        38500,
        1048,
        11685,
        1692,
        714,
        353,
        6394,
        5816,
        12254,
        344,
        8084,
        26350,
        4884,
        4807,
        2226,
        2293,
        4221,
        9206,
        13726,
        583,
        51131,
        8986,
        3196,
        60698,
        5379,
        41888,
        2817,
        9917,
        4732,
        13866,
        1942,
        5250,
        1071,
        812,
        695,
        23474,
        12719,
        8167,
        1055,
        7227,
        4104,
        4020,
        6669,
        15587,
        18752,
        5296,
        12820,
        2650,
        6527,
        15132,
        1311,
        3950,
        1519,
        24442,
        5520,
        7592,
        9050,
        936,
        43870,
        33322,
        12248,
        13064,
        2427,
        8044,
        16874,
        30882,
        9920,
        20160,
        63239,
        2947,
        8087,
        4584,
        2960,
        23501,
        8600,
        2226,
        16320,
        1164,
        21097,
        19448,
        60950,
        526,
        322,
        5496,
        16484,
        21303,
        18590,
        1060,
        6309,
        14199,
        2061,
        18075,
        444,
        4734,
        8897,
        4347,
        42504,
        5073,
        434,
        43085,
        7574,
        1044,
        23166,
        12006,
        15015,
        5670,
        6798,
        21796,
        4530,
        1108,
        934,
        1738,
        55236,
        8518,
        122,
        45214,
        4100,
        2866,
        15924,
        5379,
        7634,
        7930,
        13262,
        6174,
        383,
        2668,
        10790,
        5402,
        11360,
        22854,
        4130,
        10875
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
};
// This section is only needed prior to C++17, and can cause deprecation
// warnings if enabled after C++17
#if __cplusplus < 201703L
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<31, 1, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<31, 1, DUMMY>::table[];
#endif


}}  // end namespace

#endif
