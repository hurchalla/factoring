
#ifndef HURCHALLA_FACTORING_IS_PRIME_MONTY_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_MONTY_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#include <type_traits>
#include <limits>

namespace hurchalla { namespace factoring {


// Generally for performance, you should test for small prime factors before
// calling this function.  For example:
//    if (num%2 == 0)
//        return (num == 2);
//    if (num%3 == 0)
//        return (num == 3);
//    if (num%5 == 0)
//        return (num == 5);
//    if (num%7 == 0)
//        return (num == 7);


template <typename T, typename M>
bool is_prime_monty_trials(const T* bases, size_t total_bases, const M& mont)
{
    static_assert(std::is_same<T, typename M::T_type>::value, "");
    using V = typename M::V;
    T num = mont.getModulus();

    HPBC_PRECONDITION2(num % 2 == 1);

    V unity = mont.getUnityValue();
    HPBC_ASSERT2(unity == mont.getCanonicalForm(unity));
    V negativeOne = mont.getNegativeOneValue();
    HPBC_ASSERT2(negativeOne == mont.getCanonicalForm(negativeOne));

    // write num−1 as 2^s * d by factoring powers of 2 from num−1
    T d = static_cast<T>(num - 1);
    int s = 0;
    while (d % 2 == 0)
    {
        ++s;
        d /= 2;
    }

    for (size_t i=0; i<total_bases; ++i)
    {
        if (bases[i] == 0)
            continue;
        V result = mont.pow(mont.convertIn(bases[i]), d);
        V canonicalResult = mont.getCanonicalForm(result);
        if (canonicalResult == unity || canonicalResult == negativeOne)
            continue;
        for (int r=1; r<s; r++)
        {
            result = mont.square(result);
            canonicalResult = mont.getCanonicalForm(result);
            if (canonicalResult == negativeOne)
                break;
        }
        if (canonicalResult != negativeOne)
            return false;
    }
    return true;
}



// 34359738368 == 2^35; So far I verified the correctness of isPrime64() for all
// input values from 0 to 34359738368, looking for any mismatch compared to
// isPrime_bruteforce_optimized()
static constexpr uint64_t MAX_VERIFIED_NUMBER_MONTGOMERY_ISPRIME =
                                                          UINT64_C(34359738368);

template <typename M>
bool is_prime_monty(const M& mont)
{
    using T = typename M::T_type;
    T num = mont.getModulus();
    HPBC_PRECONDITION2(num % 2 == 1);
    HPBC_PRECONDITION2(num >= 3);
    HPBC_PRECONDITION2(std::numeric_limits<T>::digits <= 64 ||
                         num <= (T)(std::numeric_limits<uint64_t>::max()));

    uint64_t num64 = static_cast<uint64_t>(num); // the line above makes this ok

    if (num64 > std::numeric_limits<uint32_t>::max()) {
        // All bases here are smaller than num (we know num is > UINT32_MAX), so
        // there's no need to reduce any bases mod num, prior to trials.
        if (num64 <= MAX_VERIFIED_NUMBER_MONTGOMERY_ISPRIME) {
            static_assert(MAX_VERIFIED_NUMBER_MONTGOMERY_ISPRIME <
                                                    UINT64_C(154639673381), "");
            // According to http://miller-rabin.appspot.com, the 3 bases below
            // are sufficient for num up to 154639673381.  I'll use this only up
            // to the max I've verified so far, which is given by
            // MAX_VERIFIED_NUMBER_MONTGOMERY_ISPRIME.
            const T bases[] = { 15, (T)176006322, (T)4221622697 };
            size_t total_bases = sizeof(bases)/sizeof(bases[0]);
            return is_prime_monty_trials<T, M>(bases, total_bases, mont);
        }
        else {
            // According to http://miller-rabin.appspot.com, the 7 bases below
            // are sufficient for num up to 2^64.  Practically speaking, I
            // can't verify this.  Note that an alternative (which I also can't
            // verify) is to use the first 12 primes as bases for num < 2^64;
            // for details see http://oeis.org/A014233
            //
            // **TODO**
            // I've combined these two sets of bases to make it extra likely
            // that they will be sufficient for any num < 2^64.  I'll use these
            // bases for any number range I haven't brute force verified.
            // Combining the bases is almost certainly overkill for safety.
            const T bases[] = { 2, (T)325, (T)9375, (T)28178, (T)450775,
                                (T)9780504, (T)1795265022, 3, 5, 7, 11, 13, 17,
                                19, 23, 27, 31, 37 };
            size_t total_bases = sizeof(bases)/sizeof(bases[0]);
            return is_prime_monty_trials<T, M>(bases, total_bases, mont);
        }
    }
    else {   // num < 2^32
        uint32_t num32 = static_cast<uint32_t>(num);
        T bases[3];
        size_t total_bases;

#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  error HURCHALLA_TARGET_BIT_WIDTH must be defined
#endif
#if HURCHALLA_TARGET_BIT_WIDTH >= 64
        // all initialization values for bases and ranges of validity obtained
        // from http://miller-rabin.appspot.com/
        // In order to isolate the 64 bit math I've moved the (now 64 bit) mods
        // of the bases to the initialization here.  Since we know num <= 2^32,
        // any mod result is < 2^32 and thus fits in uint32_t.
        if (num32 >= UINT32_C(1050535501)) {
            total_bases = 3;
            // mod isn't needed here, since num is greater than all these bases.
            bases[0] = 2;
            bases[1] = 7;
            bases[2] = 61;    
        }
        else if (num32 >= UINT32_C(19471033)) {
            total_bases = 2;
            bases[0] = (T)(UINT64_C(336781006125) % num32);
            bases[1] = (T)(UINT64_C(9639812373923155) % num32);
        }
        else if (num32 >= UINT32_C(341531)) {
            total_bases = 2;
            // mod isn't needed here, since num is greater than all these bases.
            bases[0] = 2;
            bases[1] = (T)299417;
        }
        else {   // from the precondition we know num is still >= 3
            total_bases = 1;
            bases[0] = (T)(UINT64_C(9345883071009581737) % num32);
        }
#else
        // all initialization values for bases and ranges of validity obtained
        // from http://miller-rabin.appspot.com/ 
        if (num32 >= UINT32_C(360018361) {
            total_bases = 3;
            // these bases should cover up to 4,759,123,141 (all uint32_t).
            // mod isn't needed here, since num is greater than all these bases.
            bases[0] = 2;
            bases[1] = 7;
            bases[2] = 61;
        }
        else if (num32 >= UINT32_C(49141) {
            total_bases = 2;
            uint32_t tmp = UINT32_C(1143370);
            if (tmp >= num32)
                bases[0] = (T)(tmp % num32);
            else
                bases[0] = (T)tmp;
            bases[1] = (T)(UINT32_C(2350307676) % num32);
        }
        else {  // from the precondition we know num is still >= 3
            total_bases = 1;
            bases[0] = (T)(UINT32_C(921211727) % num32);
        }
#endif
        return is_prime_monty_trials<T, M>(bases, total_bases, mont);
    }
}


}}  // end namespace

#endif
