// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES62_5_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES62_5_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 192 byte hash table lets you determine the primality of any unsigned int
// number less than (1<<62), via miller-rabin primality testing using 5 bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<62, 5, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 64 bit number being tested for primality.
    // As a precondition, 'num' must be less than (1<<62).
    // This function returns 5 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 5> get(std::uint64_t num)
    {
        HPBC_PRECONDITION2(num < (static_cast<std::uint64_t>(1) << 62));
        // I generated/verified the hash table and bases.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 5> bases;
        bases[0] = 2;
        bases[1] = 15;
        bases[2] = 15925;
        uint32_t mask = (static_cast<uint32_t>(1) << 10) - 1u;
        uint32_t hash_bucket = ((static_cast<uint32_t>(num) & mask) * 3u) >> 6;
        bases[3] = table[hash_bucket][0];
        bases[4] = table[hash_bucket][1];
        return bases;
    }
private:
    static constexpr std::size_t SIZE = 48;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE][2] = {
#else
    static constexpr std::uint16_t table[][2] = {
#endif
        {  7400, 12248 },
        { 18142, 46722 },
        {  7105,  6558 },
        {  8517, 30461 },
        { 40909, 36007 },
        { 53428, 63769 },
        { 45870, 49355 },
        { 48490, 28834 },
        { 23125, 26491 },
        { 14913, 39742 },
        { 54899, 27103 },
        { 11956, 26386 },
        { 25327, 38708 },
        { 14437, 57042 },
        {  8829, 26770 },
        {  9692, 31894 },
        { 53428, 24244 },
        {   783,  9426 },
        {   111, 14519 },
        { 40916, 38390 },
        { 14376, 23291 },
        {  2402, 65414 },
        { 21756, 51357 },
        { 21312, 41342 },
        { 50207, 17271 },
        { 33901, 60307 },
        { 27750, 60415 },
        { 14703, 35020 },
        { 44506, 24395 },
        {  1305, 60869 },
        { 40695, 55758 },
        {  4313, 55119 },
        { 54900, 35982 },
        { 17393, 38355 },
        { 19252,  9599 },
        { 23441, 17641 },
        { 33958, 53674 },
        { 45983, 14910 },
        { 16470,  9666 },
        { 29970, 30487 },
        { 48979, 41983 },
        {   549, 60222 },
        { 47530,  5934 },
        { 47593, 39383 },
        { 62197, 24797 },
        { 46481, 49763 },
        { 20940, 59710 },
        { 18759, 58580 }
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
    static_assert(sizeof(table[0])/sizeof(table[0][0]) == 2u, "");
    static_assert(sizeof(table)/sizeof(table[0][0]) == SIZE*2, "");
};
// This section is only needed prior to C++17, and can cause deprecation
// warnings if enabled after C++17
#if __cplusplus < 201703L
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<62, 5, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<62, 5, DUMMY>::table[][2];
#endif


}}  // end namespace

#endif
