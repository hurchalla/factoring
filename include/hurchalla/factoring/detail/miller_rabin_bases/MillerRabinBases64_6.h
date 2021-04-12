// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES64_6_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES64_6_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 40 byte hash table lets you determine the primality of any 64 bit
// unsigned int number, via miller-rabin primality testing using 6 bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<64, 6, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 64 bit number being tested for primality.
    // This function returns 6 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 6> get(std::uint64_t num)
    {
        // I generated/verified these hash tables and bases.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 6> bases;
        bases[0] = 2;
        bases[1] = 15;
        bases[2] = 925;
        bases[3] = 28717;
        // Either version below should work fine.  I verified both are correct.
        // Typically compilers seem to produce slightly smaller machine code
        // with the first version.
#if 1
        constexpr std::uint16_t table1[] = { 29385, 23886, 3029, 31566 };
        constexpr std::uint16_t table2[] = {
                        16307, 15537, 17271, 56697, 46229, 16499, 2098, 62203,
                        38699, 27257, 60153, 44757, 29306, 3257, 42962, 24003 };
        uint32_t hash_bucket1 = static_cast<uint32_t>(num) >> 30;
        uint32_t hash_bucket2 = static_cast<uint32_t>(num) >> 28;
        bases[4] = table1[hash_bucket1];
        bases[5] = table2[hash_bucket2];
#else
        bases[4] = 52226;
        constexpr std::uint16_t table[] = { 58031, 20842, 58103, 42934, 52270,
                                            57725,  5873, 60766, 64507, 47847,
                                            38101, 44554, 32042,  2603, 29374,
                                            22234, 19658,  5829, 37303, 41046 };
        static_assert(sizeof(table)/sizeof(table[0]) == 20, "");
        uint32_t mask = (static_cast<uint32_t>(1) << 22) - 1u;
        uint32_t hash_bucket = ((static_cast<uint32_t>(num) & mask) * 5) >> 20;
        bases[5] = table[hash_bucket];
#endif
        return bases;
    }
};


}}  // end namespace

#endif
