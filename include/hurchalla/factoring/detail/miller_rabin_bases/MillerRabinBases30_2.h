// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES30_2_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES30_2_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 8 byte hash table lets you determine the primality of any unsigned int
// number less than (1<<30), via miller-rabin primality testing using 2 bases.
// This function can of course test an odd number for primality; it can also
// test an even number if the macro HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS is
// defined.  Note however that montgomery arithmetic (which is one way to
// implement the miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<30, 2, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 32 bit number being tested for primality.
    // As a precondition, 'num' must be less than (1<<30).
    // This function returns 2 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 2> get(std::uint32_t num)
    {
        HPBC_PRECONDITION2(num < (static_cast<std::uint32_t>(1) << 30));
        // I generated/verified the hash table and bases.  See README.TXT
        std::array<std::uint16_t, 2> bases;
#ifdef HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS
        bases[0] = 30;
        const std::array<std::uint16_t, 4> table =
                         { 4418, 54365, 18, 52797 };
        std::uint32_t hash_bucket = (num ^ (num >> 1)) & 3;
#else
        // This section uses a simpler hash function, and is correct for all odd
        // numbers.  However it fails for the number 4 (it's correct for all
        // other even numbers).  Since HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS
        // is not defined, this isn't a problem - we only allow odd inputs.
        HPBC_PRECONDITION2(num % 2 == 1);
        bases[0] = 42685;
        const std::array<std::uint16_t, 4> table =
                         { 38165, 50768, 59722, 23646 };
        std::uint32_t hash_bucket = (num >> 7) & 3;
#endif
        bases[1] = table[hash_bucket];
        return bases;
    }
};


}}  // end namespace

#endif
