// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_SIEVE_OF_ERATOSTHENES_H_INCLUDED
#define HURCHALLA_FACTORING_SIEVE_OF_ERATOSTHENES_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <vector>
#include <cstdint>
#include <memory>

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

// The version further below avoids the main performance problem of the above
// versions, which is that they skip through memory and rarely reuse CPU cache
// lines prior to cache evictions.  To avoid this problem, the next version
// works on a range of memory smaller than the CPU cache, and writes all values
// that will ever be needed for that block, and then moves on to the next block
// of memory and writes all values needed for that block, etc.
// I measured ~4x performance improvement on Intel Haswell CPU using the next
// version, compared to the above versions.



class SieveBitVector {
// In most respects, having a std::vector<bool> member do this work is
// preferable to having custom bit vector logic, but unfortunately on 32 bit
// architectures std::vector<bool> doesn't always permit a large enough vector
// to be created for our uses (std::vector<bool>::max_size() might be 2^31 - 1,
// which isn't large enough.  Note that 2^31 - 1 is a much smaller limit to the
// vector<bool> length than a 32bit architecture would actually require, given
// that vector<bool> uses approximately vectorlength/8 bytes of memory; for
// example, a length of 2^31 would need only ~256MB).
#if HURCHALLA_TARGET_BIT_WIDTH > 32
    std::vector<bool> vb;
public:
    SieveBitVector(std::uint32_t count, bool value) : vb(count, value) {}
    bool get(std::uint32_t index) const { return vb[index]; }
    void clear(std::uint32_t index) { vb[index] = false; }
#else
    std::uint32_t size;
    std::unique_ptr<unsigned char[]> membytes;
public:
    SieveBitVector(std::uint32_t count, bool value) :
                               size((count%8 == 0) ? count/8 : (count/8)+1),
                               membytes(std::make_unique<unsigned char[]>(size))
    {
        unsigned char val = static_cast<unsigned char>((value) ? 255 : 0);
        std::fill(membytes.get(), membytes.get() + size, val);
    }
    bool get(std::uint32_t index) const
    {
        HPBC_PRECONDITION(index < static_cast<std::uint64_t>(size)*8);
        std::uint32_t bytenum = index/8;
        std::uint8_t offset = static_cast<std::uint8_t>(index % 8);
        return (membytes[bytenum] >> offset) & 1;
    }
    void clear(std::uint32_t index)
    {
        using U = std::uint8_t;
        HPBC_PRECONDITION(index < static_cast<std::uint64_t>(size)*8);
        std::uint32_t bytenum = index/8;
        U offset = static_cast<U>(index % 8);
        U on_mask = static_cast<U>(static_cast<U>(1) << offset);
        U off_mask = static_cast<U>(~on_mask);
        membytes[bytenum] = static_cast<U>(membytes[bytenum] & off_mask);
    }
#endif
};


// Returns a bit vector with indices that represent odd numbers.  Every true
// entry in the bit vector means that the odd number represented by the index is
// prime, and every false entry means that the odd number represented by the
// index is not prime.
//
// *Note that the bit vector this function creates and returns will use
// size_odds/8 bytes of memory.  E.g. if size_odds == 1<<31, the vector will
// take up 256 MB.
//
// we use DUMMY because we want the function to be inline to avoid ODR errors,
// but gcc warns (-Winline) that it can't inline the function if we explicitly
// make it inline.  Using a template implicitly makes it inline, without any
// warning.
template <typename DUMMY=void>
SieveBitVector init_sieve_odd_primes(std::uint32_t size_odds,
                                       std::uint64_t cache_blocking_size)
{
    using std::uint64_t;
    using std::uint32_t;
    uint64_t size = static_cast<uint64_t>(size_odds)*2;
    HPBC_PRECONDITION(cache_blocking_size > 0);
    HPBC_PRECONDITION(2 <= size);

    SieveBitVector primes_bitvec(size_odds, true);
    primes_bitvec.clear(1/2);  // the value 1 is not a prime

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
        HPBC_ASSERT2(i/2 < size_odds);
        if (primes_bitvec.get(static_cast<uint32_t>(i/2))) {
            uint64_t j=i*i;
            for (; j*j<size; j+=(2*i)) {
                HPBC_ASSERT2(j/2 < size_odds);
                primes_bitvec.clear(static_cast<uint32_t>(j/2));
            }
            HPBC_ASSERT2(2*i <= ut_numeric_limits<uint32_t>::max());
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
        using pd_size_type = decltype(prime_doubled_vec)::size_type;
        for (pd_size_type j=0; j<prime_doubled_vec.size(); ++j) {
            uint64_t multiple = prime_multiple_vec[j];
            HPBC_ASSERT2(multiple >= i);
            uint64_t prime_doubled = prime_doubled_vec[j];
            for (; multiple<next; multiple+=prime_doubled) {
                HPBC_ASSERT2(multiple/2 < size_odds);
                primes_bitvec.clear(static_cast<uint32_t>(multiple/2));
            }
            prime_multiple_vec[j] = multiple;
        }
    }
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic pop
#endif
    return primes_bitvec;
}


class SieveOfEratosthenes {
    const SieveBitVector oddprimes;
    const std::uint64_t length;

public:
    SieveOfEratosthenes(std::uint64_t size,
                        std::uint64_t cache_blocking_size = 262144)
          : oddprimes(init_sieve_odd_primes(static_cast<std::uint32_t>(size/2),
                                            cache_blocking_size)),
            length(size)
    {
        HPBC_PRECONDITION(size % 2 == 0);
        HPBC_PRECONDITION(2 <= size);
        HPBC_PRECONDITION(size/2 <= ut_numeric_limits<uint32_t>::max());
    }

    std::uint64_t size() const { return length; }

    bool isPrime(std::uint64_t value) const
    {
        HPBC_PRECONDITION2(value < length);
        if (value % 2 == 0)
            return (value == 2);
        else {
            HPBC_ASSERT2(value/2 <= ut_numeric_limits<std::uint32_t>::max());
            return oddprimes.get(static_cast<std::uint32_t>(value/2));
        }
    }

    bool operator[](std::uint64_t index) const { return isPrime(index); }
};


}}

#endif
