// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES30_1_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES30_1_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 320 byte hash table lets you determine the primality of any unsigned int
// number less than (1<<30), via miller-rabin primality testing using 1 single
// base.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<30, 1, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 32 bit number being tested for primality.
    // As a precondition, 'num' must be less than (1<<30).
    // This function returns a single base that can be used by miller-rabin
    // testing to correctly determine (non-probabilistically) the primality of
    // num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 1> get(std::uint32_t num)
    {
        HPBC_PRECONDITION2(num < (static_cast<std::uint32_t>(1) << 30));
        // I generated/verified the hash table for the base.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 1> bases;
        uint32_t mask = (static_cast<uint32_t>(1) << 14) - 1u;
        uint32_t hash_bucket = (((9u*num) & mask) * 5u) >> 9;
        bases[0] = table[hash_bucket];
        return bases;
    }
private:
    static constexpr std::size_t SIZE = 160;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE] = {
#else
    static constexpr std::uint16_t table[] = {
#endif
        3982,
        18006,
        21494,
        2005,
        1838,
        5822,
        1251,
        13402,
        7804,
        930,
        224,
        3949,
        41695,
        552,
        312,
        2734,
        212,
        15580,
        1920,
        454,
        312,
        1169,
        2760,
        3182,
        4672,
        1488,
        9144,
        11158,
        11020,
        6880,
        8516,
        1128,
        6118,
        706,
        3040,
        35,
        14331,
        1144,
        278,
        2180,
        20618,
        5212,
        4326,
        7432,
        206,
        8784,
        38591,
        9342,
        2073,
        5372,
        1841,
        1244,
        5552,
        5684,
        5739,
        5795,
        1563,
        399,
        2622,
        5282,
        3501,
        24732,
        7735,
        3814,
        515,
        115,
        2402,
        1072,
        6295,
        11424,
        11524,
        18347,
        2240,
        20722,
        2087,
        2857,
        1651,
        4775,
        7453,
        733,
        1503,
        7791,
        171,
        264,
        377,
        4043,
        1526,
        3787,
        370,
        7280,
        518,
        2062,
        7523,
        7310,
        571,
        4076,
        7847,
        1664,
        2779,
        407,
        7358,
        7714,
        1120,
        852,
        5079,
        982,
        18669,
        3344,
        9934,
        2688,
        16602,
        1352,
        19686,
        2283,
        1085,
        21702,
        2198,
        268,
        394,
        17184,
        13168,
        21192,
        4702,
        1392,
        2067,
        3897,
        8190,
        311,
        341,
        2532,
        5093,
        3972,
        184,
        411,
        4655,
        3398,
        45648,
        1869,
        2418,
        2392,
        1518,
        42028,
        9582,
        4606,
        2135,
        3981,
        15390,
        9040,
        2574,
        2426,
        2122,
        10222,
        3183,
        1338,
        27412,
        562,
        3170,
        8256,
        5887,
        20538
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
};
// This section is only needed prior to C++17, and can cause deprecation
// warnings if enabled after C++17
#if __cplusplus < 201703L
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<30, 1, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<30, 1, DUMMY>::table[];
#endif


}}  // end namespace

#endif
