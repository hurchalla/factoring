// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_SIEVE_OF_ERATOSTHENES_H_INCLUDED
#define HURCHALLA_FACTORING_SIEVE_OF_ERATOSTHENES_H_INCLUDED


#include "hurchalla/util/programming_by_contract.h"
#include <vector>
#include <cstdint>
#include <limits>

namespace hurchalla { namespace detail {


// For details on the sieve of Eratosthenes, see
// https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes

// To help make the later code easier to understand, I have left earlier/simple
// versions here.  These disabled functions provide a progression from the
// simplest and easiest to understand version up to the harder to understand and
// ordinarily fastest version (which is init_sieve_odd_primes).
#if 0
// Initialize a bit vector that marks all prime indices as true up to size,
// and marks all non-prime indices as false.
inline std::vector<bool> init_primes(std::uint64_t size)
{
    using std::uint64_t;
    std::vector<bool> primes(size, true);
    primes[0] = false;
    primes[1] = false;
    for (uint64_t j=2*2; j<size; j+=2)
        primes[j] = false;
    for (uint64_t i=3; i*i<size; i+=2) {
        if (primes[i]) {
            for (uint64_t j=i*i; j<size; j+=(2*i))
                primes[j] = false;
        }
    }
    return primes;
}
#endif
#if 0
// Similar to the above, except the bit vector indices represent odd numbers.
// Thus the resulting bit vector is half the size of the above.
inline std::vector<bool> init_odd_primes_simple(std::uint64_t size)
{
    using std::uint64_t;
    HPBC_PRECONDITION(2 <= size);
    uint64_t size_odds = size/2;
    std::vector<bool> primes(size_odds, true);
    primes[1/2] = false;

    for (uint64_t i=3; i*i<size; i+=2) {
        if (primes[i/2]) {
            for (uint64_t j=i*i; j<size; j+=(2*i)) {
                HPBC_ASSERT2(j/2 < size_odds);
                primes[j/2] = false;
            }
        }
    }
    return primes;
}
#endif

// The version below avoids the main performance problem of the above versions,
// which is that they skip through memory and rarely reuse CPU cache lines prior
// to cache evictions.  To avoid this problem, this version works on a range of
// memory smaller than the CPU cache, and writes all values that will ever be
// needed for that block, and then moves on to the next block of memory and
// writes all values needed for that block, etc.
// I measured ~4x performance improvement on Intel Haswell CPU using this
// version, compared to the above versions.


// Returns a bit vector with indices that represent odd numbers.  Every true
// entry in the bit vector means that the odd number represented by the index is
// prime, and every false entry means that the odd number represented by the
// index is not prime.

// *Note that the bit vector this function creates and returns will be (size/16)
// bytes large.  E.g. for a size of 1<<32, the vector will take up 256 MB.
// If you go too much above 1<<32, the bit vector is going to start taking an
// insane amount of memory (or run out of memory), and calculating the sieve
// here will take a very long time (e.g. for size 1<< 36 you might expect over a
// minute).
// Also note that on some systems, the vector's size_type might be 32 bit, which
// would mean by one of function's preconditions below you'd be limited to a
// size < (1<<33).  In general I'd recommend staying below that limit for
// portability and to keep memory usage/computation time reasonable.

// we use DUMMY because we want the function to be inline to avoid ODR errors,
// but gcc warns (-Winline) that it can't inline the function if we explicitly
// make it inline.  Using a template implicitly makes it inline, without any
// warning.
template <typename DUMMY=void>
std::vector<bool> init_sieve_odd_primes(std::uint64_t size,
                                        std::uint64_t cache_blocking_size)
{
    using vec_size_type = std::vector<bool>::size_type;
    HPBC_PRECONDITION(2 <= size);
    HPBC_PRECONDITION(size/2 <= std::numeric_limits<vec_size_type>::max());
    HPBC_PRECONDITION(cache_blocking_size > 0);

    using std::uint64_t;
    using std::uint32_t;
    vec_size_type size_odds = size/2;
    std::vector<bool> primes_bitvec(size_odds, true);
    primes_bitvec[1/2] = false;

    std::vector<uint64_t> prime_multiple_vec;
    std::vector<uint32_t> prime_doubled_vec;

    // the primes < sqrt(size) are special, in that they are all that we need to
    // filter out all composites >= sqrt(size).  We store info related to them
    // in prime_multiple_vec and prime_doubled_vec, and use that info later to
    // mark the composites in the bit vector.
    // Later, for each block of memory we will process, prime_multiple_vec[i]
    // will tell us the value of the first odd multiple of its associated prime
    // (the associated prime == prime_doubled_vec[i]/2).  prime_multiple_vec[i]
    // gets updated every block.  prime_doubled_vec[i] remains unchanged.
    uint64_t i=3;
    // find all primes and composites < sqrt(size)
    for (; i*i<size; i+=2) {
        if (primes_bitvec[i/2]) {
            uint64_t j=i*i;
            for (; j*j<size; j+=(2*i)) {
                HPBC_ASSERT2(j/2 < size_odds);
                primes_bitvec[j/2] = false;
            }
            HPBC_ASSERT2(2*i < (static_cast<uint64_t>(1) << 32));
            prime_doubled_vec.push_back(static_cast<uint32_t>(2*i));
            prime_multiple_vec.push_back(j);
        }
    }
    HPBC_ASSERT2(prime_doubled_vec.size() == prime_multiple_vec.size());

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunsafe-loop-optimizations"
#endif
    // mark all the composites that are >= sqrt(size).  After loop exit, any
    // entries in primes_bitvec that are still unmarked (i.e. that are still
    // left as true) represent prime numbers.
    for (; i<size; i+=cache_blocking_size) {
        uint64_t next = i+cache_blocking_size;
        if (next > size)
            next = size;
        for (uint64_t j=0; j<prime_doubled_vec.size(); ++j) {
            uint64_t prime_multiple = prime_multiple_vec[j];
            HPBC_ASSERT2(prime_multiple >= i);
            uint64_t prime_doubled = prime_doubled_vec[j];
            for (; prime_multiple<next; prime_multiple+=prime_doubled) {
                HPBC_ASSERT2(prime_multiple/2 < size_odds);
                primes_bitvec[prime_multiple/2] = false;
            }
            prime_multiple_vec[j] = prime_multiple;
        }
    }
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic pop
#endif
    return primes_bitvec;
}


class SieveOfEratosthenes {
    const std::vector<bool> oddprimes;
    const std::uint64_t length;

public:
    SieveOfEratosthenes(std::uint64_t size,
                        std::uint64_t cache_blocking_size = 262144)
            : oddprimes(init_sieve_odd_primes(size, cache_blocking_size)),
              length(size)
    {
        HPBC_ASSERT2(length/2 == oddprimes.size());
    }

    std::uint64_t size() const { return length; }

    bool isPrime(std::uint64_t value) const
    {
        HPBC_PRECONDITION2(value < length);
        if (value % 2 == 0)
            return (value == 2);
        else
            return oddprimes[value/2];
    }

    bool operator[](std::uint64_t index) const { return isPrime(index); }
};


}}

#endif
