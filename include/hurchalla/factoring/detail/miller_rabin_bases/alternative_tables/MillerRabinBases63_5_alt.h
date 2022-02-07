// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES63_5_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES63_5_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This alternative version uses slightly smaller tables in total than the
// normal version, but its hash function is more complicated, and its get()
// function must read two cache lines (one for each of the two tables) instead
// of the single cache line access needed by the normal version.

// These 16+208 byte hash tables let you determine the primality of any unsigned
// int number less than (1<<63), via miller-rabin primality testing using 5
// bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<63, 5, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 64 bit number being tested for primality.
    // As a precondition, 'num' must be less than (1<<63).
    // This function returns 5 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 5> get(std::uint64_t num)
    {
        // I generated/verified these hash tables and bases.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 5> bases;
        bases[0] = 2;
        bases[1] = 15;
        bases[2] = 15925;
        constexpr std::uint16_t table1[] = { 9767, 8690, 54574, 824,
                                             37691, 2714, 16397, 35794 };
        uint32_t mask = (static_cast<uint32_t>(1) << 28) - 1u;
        uint32_t hash_bucket2 = ((static_cast<uint32_t>(num) & mask) * 13) >>25;
        std::uint32_t hash_bucket1 = hash_bucket2 & 7;
        bases[3] = table1[hash_bucket1];
        bases[4] = table2[hash_bucket2];
        return bases;
    }
private:
    static constexpr std::size_t TABLE2_SIZE = 104;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table2[TABLE2_SIZE] = {
#else
    static constexpr std::uint16_t table2[] = {
#endif
        36879,
        22602,
        47414,
        6861,
        65417,
        53109,
        26902,
        8030,
        45297,
        18462,
        22587,
        38269,
        15854,
        13417,
        21101,
        12738,
        45835,
        12702,
        9773,
        45329,
        7063,
        20290,
        38251,
        30353,
        38050,
        19858,
        10021,
        32498,
        8961,
        24611,
        12609,
        24770,
        61302,
        14102,
        6194,
        35463,
        14653,
        20187,
        12247,
        1577,
        33438,
        31718,
        23233,
        369,
        15933,
        53449,
        34397,
        2845,
        14110,
        14037,
        39625,
        21229,
        51459,
        40813,
        18749,
        6753,
        654,
        6611,
        58113,
        27846,
        4353,
        23325,
        10835,
        62511,
        513,
        6006,
        46685,
        14945,
        58706,
        6743,
        10197,
        13330,
        20097,
        24485,
        61803,
        15734,
        349,
        6182,
        49279,
        52251,
        40439,
        38869,
        39387,
        7094,
        34389,
        56793,
        54234,
        15338,
        22790,
        8779,
        25251,
        14694,
        40839,
        18953,
        33183,
        16179,
        20729,
        5458,
        50121,
        54775,
        8943,
        26086,
        15194,
        10518
    };
    static_assert(sizeof(table2)/sizeof(table2[0]) == TABLE2_SIZE, "");
};
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<63, 5, DUMMY>::TABLE2_SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<63, 5, DUMMY>::table2[];


}}  // end namespace

#endif
