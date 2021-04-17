// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES64_5_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES64_5_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 320 byte hash table lets you determine the primality of any 64 bit
// unsigned int number, via miller-rabin primality testing using 5 bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<64, 5, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 64 bit number being tested for primality.
    // This function returns 5 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 5> get(std::uint64_t num)
    {
        // I generated/verified the hash table and bases.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 5> bases;
        bases[0] = 2;
        bases[1] = 15;
        bases[2] = 52;
        uint32_t mask = (static_cast<uint32_t>(1) << 19) - 1u;
        uint32_t hash_bucket = ((static_cast<uint32_t>(num) & mask) * 5u) >> 15;
        bases[3] = table[hash_bucket][0];
        bases[4] = table[hash_bucket][1];
        return bases;
    }
private:
    static constexpr std::size_t SIZE = 80;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE][2] = {
#else
    static constexpr std::uint16_t table[][2] = {
#endif
        {  5477,  8321 },
        { 14011, 36773 },
        {  1313, 55018 },
        {  1765, 28577 },
        { 15481, 31351 },
        {  7205, 31499 },
        {  2909, 47379 },
        {  1219, 40643 },
        {   747, 51098 },
        {  3611, 42294 },
        {   481,  2157 },
        { 13918, 19993 },
        {   929,  9723 },
        {   641, 46917 },
        {  7885, 22273 },
        { 10311, 35865 },
        {  3318, 57993 },
        { 10862, 11902 },
        { 12749, 48323 },
        {  5433, 39434 },
        {   979, 21139 },
        {  1114, 22514 },
        {   279, 61962 },
        {  9957, 60319 },
        {   571, 27607 },
        {  1205,  9145 },
        {  4450, 57955 },
        {  1870, 22181 },
        {   978, 28059 },
        {  2474, 18093 },
        {   209, 47862 },
        {   710, 53001 },
        {    46,  1522 },
        {  1869, 44347 },
        {   199, 36241 },
        {  2201, 62787 },
        {   873, 43378 },
        {   818, 18998 },
        {  8782, 61778 },
        {  3457, 42963 },
        {  3947, 37219 },
        {  3671, 41738 },
        {  8553, 61969 },
        {  2550, 25363 },
        {  2339, 19691 },
        {  6207,  9719 },
        {   749, 32245 },
        {   233, 45209 },
        {  6779, 29361 },
        {  3611, 28679 },
        {  2230, 44379 },
        { 11983, 38026 },
        {   973, 47706 },
        {  3495, 27177 },
        {  8846, 41123 },
        { 17326, 19729 },
        {  4297, 30885 },
        {   761, 59549 },
        {  9099, 55766 },
        { 16734, 55099 },
        {  3778, 49627 },
        {  1658, 41254 },
        {   109, 19265 },
        {  2015, 17813 },
        {    62, 44441 },
        {  7411, 62994 },
        {   262, 56374 },
        { 15394, 45407 },
        {  1114, 40359 },
        {   505, 33886 },
        {  2073, 59671 },
        {  3818, 32269 },
        {  8517, 20673 },
        {   569, 32277 },
        {  2622, 18819 },
        {   835, 64263 },
        {  1671, 39043 },
        {   119, 23294 },
        {  1046, 63910 },
        {  1643, 56118 }
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
    static_assert(sizeof(table[0])/sizeof(table[0][0]) == 2u, "");
    static_assert(sizeof(table)/sizeof(table[0][0]) == SIZE*2, "");
};
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<64, 5, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<64, 5, DUMMY>::table[][2];


}}  // end namespace

#endif
