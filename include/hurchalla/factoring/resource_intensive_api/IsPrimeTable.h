// Copyright (c) 2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IS_PRIME_TABLE_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_TABLE_H_INCLUDED


#if !defined(HURCHALLA_FACTORING_DISALLOW_INLINE_ASM) && \
        !defined(HURCHALLA_ALLOW_INLINE_ASM_ALL)
#  define HURCHALLA_ALLOW_INLINE_ASM_ALL
#endif


#include "hurchalla/factoring/detail/SieveOfEratosthenes.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla {


// operator()(x) returns true if x is prime, and otherwise returns false.
//
// T can be any integral type <= 32 bits.
//
// This functor is intended for use when you plan to repeatedly and intensively
// perform primality testing.  The constructor creates a table (a Sieve of
// Eratosthenes) in memory, and depending on the size of this table it may take
// a few seconds to construct.  Excluding the time it takes to construct, this
// functor provides what is most likely the fastest primality testing possible
// for any type T <= 32 bits, at the cost of using far more memory (and CPU
// cache) than the normal API function is_prime() would have used.
//
// The amount of memory used by the table is as follows for these types T
// uint8_t: 16 Bytes
// uint16_t: 4 Kilobytes
// uint32_t: 256 Megabytes

template <typename T>
struct IsPrimeTable {
private:
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::digits <= 32, "");

    const detail::SieveOfEratosthenes sieve;

public:
    IsPrimeTable() :
              sieve(static_cast<uint64_t>(1) << ut_numeric_limits<T>::digits) {}

    bool operator()(T x) const
    {
        HPBC_PRECONDITION2(x >= 0);
        return sieve.isPrime(x);
    }
};


}  // end namespace

#endif
