// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES44_3_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES44_3_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 64 byte hash table lets you determine the primality of any unsigned
// integer less than (1<<44), via miller-rabin primality testing using 3 bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<44, 3, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 64 bit integer being tested for primality.
    // As a precondition, 'num' must be less than (1<<44).
    // This function returns 3 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 3> get(std::uint64_t num)
    {
        HPBC_PRECONDITION2(num < (static_cast<std::uint64_t>(1) << 44));
        // I generated/verified the hash table and bases.  See README.TXT
        using std::uint16_t;
        std::array<std::uint16_t, 3> bases;
        bases[0] = 2;
        uint16_t hash_bucket = static_cast<uint16_t>(
                                              static_cast<uint16_t>(num) >> 12);
        bases[1] = table[hash_bucket][0];
        bases[2] = table[hash_bucket][1];
        return bases;

//shift == 12
//uint64_t mask = (static_cast<uint64_t>(1) << 16) - 1u;
//uint32_t hash_bucket2 = static_cast<uint32_t>((psp & mask) >> shift);
    }
private:
    static constexpr std::size_t SIZE = 16;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE][2] = {
#else
    static constexpr std::uint16_t table[][2] = {
#endif
        { 16689, 28124 },
        {  2885, 25943 },
        { 32559, 42064 },
        { 65533, 34531 },
        { 52215,  5179 },
        {    60, 46618 },
        { 33135, 19528 },
        { 42237, 57254 },
        {   936, 25451 },
        {   540, 51929 },
        { 38381, 45834 },
        {  4290, 38983 },
        {  4499, 26821 },
        { 62529, 40449 },
        { 12933, 45642 },
        { 40101, 61627 }
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
    static_assert(sizeof(table[0])/sizeof(table[0][0]) == 2u, "");
    static_assert(sizeof(table)/sizeof(table[0][0]) == SIZE*2, "");
};
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<44, 3, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<44, 3, DUMMY>::table[][2];


}}  // end namespace

#endif
