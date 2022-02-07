// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

// purposely placed outside the include guard
#ifdef HURCHALLA_CHOOSE_ALTERNATE_MILLER_RABIN_BASES63_5
#include "hurchalla/factoring/detail/miller_rabin_bases/alternative_tables/MillerRabinBases63_5_alt.h"
#endif
#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES63_5_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES63_5_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 256 byte hash table lets you determine the primality of any unsigned
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
        HPBC_PRECONDITION2(num < (static_cast<std::uint64_t>(1) << 63));
        // I generated/verified the hash table and bases.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 5> bases;
        bases[0] = 2;
        bases[1] = 15;
        bases[2] = 15925;
        uint32_t hash_bucket = static_cast<uint32_t>(num) >> 26;
        bases[3] = table[hash_bucket][0];
        bases[4] = table[hash_bucket][1];
        return bases;
    }
private:
    static constexpr std::size_t SIZE = 64;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE][2] = {
#else
    static constexpr std::uint16_t table[][2] = {
#endif
        {  8954, 64619 },
        { 13431, 23949 },
        { 38035, 53037 },
        { 15400, 29509 },
        {   783, 46411 },
        { 26714, 49375 },
        { 61831, 58599 },
        { 59284, 42994 },
        { 39960, 43365 },
        { 56133, 34301 },
        { 61831, 56043 },
        { 23414, 54204 },
        { 13636, 46459 },
        { 60034, 60075 },
        { 64260, 49575 },
        {  7105, 51693 },
        { 17655, 51362 },
        { 19105, 29981 },
        { 41688, 62154 },
        {  8833, 19721 },
        { 60173, 50469 },
        { 19526, 31593 },
        { 22337,  9596 },
        { 34498, 64439 },
        { 10722, 22905 },
        { 18717, 39939 },
        { 63945, 64334 },
        { 46582, 31239 },
        { 45325, 53241 },
        { 50786,  7131 },
        { 43628, 54314 },
        {  1155, 29495 },
        { 53367, 19929 },
        { 34420, 34207 },
        { 35845, 45484 },
        {  4239, 19491 },
        { 37590, 30139 },
        { 21312, 36231 },
        { 52732, 12266 },
        { 21433, 39019 },
        { 18779, 58255 },
        { 49410, 38033 },
        { 10693, 53154 },
        { 30011, 29948 },
        { 34608, 33656 },
        { 27925, 61446 },
        { 11224, 43074 },
        { 35073, 37870 },
        { 59200, 25927 },
        {  1744, 57208 },
        { 59940, 52786 },
        { 63936, 52475 },
        { 39146, 40245 },
        { 21024, 54244 },
        { 53428, 18989 },
        { 25951, 20419 },
        { 31117, 29434 },
        { 16398,  5726 },
        { 55635, 53716 },
        { 17818, 34641 },
        { 62197,  2798 },
        { 20206, 43233 },
        {  9760, 30961 },
        { 32269, 30035 }
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
    static_assert(sizeof(table[0])/sizeof(table[0][0]) == 2u, "");
    static_assert(sizeof(table)/sizeof(table[0][0]) == SIZE*2, "");
};
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<63, 5, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<63, 5, DUMMY>::table[][2];


}}  // end namespace

#endif
