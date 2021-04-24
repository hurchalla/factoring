// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES62_6_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES62_6_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 20 byte hash table lets you determine the primality of any unsigned int
// number less than (1<<62), via miller-rabin primality testing using 6 bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<62, 6, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 64 bit number being tested for primality.
    // As a precondition, 'num' must be less than (1<<62).
    // This function returns 6 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 6> get(std::uint64_t num)
    {
        HPBC_PRECONDITION2(num < (static_cast<std::uint64_t>(1) << 62));
        // I generated/verified these hash tables and bases.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 6> bases;
        bases[0] = 2;
        bases[1] = 15;
        bases[2] = 925;
        bases[3] = 28717;
        bases[4] = 3727;
        constexpr std::uint16_t table[] = { 65186, 1983, 2557, 49382, 19999,
                                            6218, 51695, 6637, 43774, 14137 };
        static_assert(sizeof(table)/sizeof(table[0]) == 10, "");
        uint32_t mask = (static_cast<uint32_t>(1) << 17) - 1u;
        uint32_t hash_bucket = ((static_cast<uint32_t>(num) & mask) * 5) >> 16;
        bases[5] = table[hash_bucket];
        return bases;
    }
};


}}  // end namespace

#endif
