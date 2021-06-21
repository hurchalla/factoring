// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IS_PRIME_INTENSIVE_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_INTENSIVE_H_INCLUDED


#include "hurchalla/factoring/detail/ImplIsPrimeIntensive.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"

namespace hurchalla {


// T can be any unsigned integral type <= 128 bits.
//
// This functor is intended for use when you plan to repeatedly and intensively
// perform primality testing.  For some values of T (only uint32_t at the
// moment) the constructor will likely take a couple of seconds to complete.
// Once the object has been constructed, the functor will provide the fastest
// possible primality testing, at the cost of often using far more memory (and
// CPU cache) than the alternative function is_prime() would have used.

// At present, the memory usage of this functor is as follows for different
// template types T -
// uint8_t: 16 Bytes
// uint16_t: 4 Kilobytes
// uint32_t: 256 Megabytes
// uint64_t: ~1.3 Kilobytes (with OPTIMIZE_PRIMES == false)
// uint64_t: 448 Kilobytes (with OPTIMIZE_PRIMES == true)
// __uint128_t: ~2 Kilobytes

// The template parameter OPTIMIZE_PRIMES controls whether the functor will be
// optimized for the situation where the numbers that you test will most often
// be prime.  If this is not the case for you, use the default false.
template <typename T, bool OPTIMIZE_PRIMES = false>
class IsPrimeIntensive {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits <= 128, "");
    const detail::ImplIsPrimeIntensive<T, OPTIMIZE_PRIMES> isprime_impl;
public:
    IsPrimeIntensive() : isprime_impl() {}

    // returns true if x is prime, otherwise returns false
    bool operator()(T x) const
    {
        return isprime_impl(x);
    }
};


}  // end namespace

#endif
