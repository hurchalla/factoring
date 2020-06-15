
#ifndef HURCHALLA_FACTORING_IS_PRIME_MILLER_RABIN_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_MILLER_RABIN_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#include <type_traits>
#include <limits>
#include <memory>

namespace hurchalla { namespace factoring {


// This algorithm is the deterministic variant of the Miller-Rabin primailty
// test.  See:
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants
// subheading "Testing against small sets of bases".
// See also:
// https://miller-rabin.appspot.com/
template <typename T, typename MontType>
bool is_prime_mr_trials(const T* bases,size_t total_bases, const MontType& mont)
{
    static_assert(std::is_same<T, typename MontType::T_type>::value, "");
    T num = mont.getModulus();

    HPBC_ASSERT2(num % 2 == 1);

    auto unity = mont.getUnityValue();
    HPBC_ASSERT2(unity == mont.getCanonicalForm(unity));
    auto negativeOne = mont.getNegativeOneValue();
    HPBC_ASSERT2(negativeOne == mont.getCanonicalForm(negativeOne));

    // write num−1 as 2^r * d by factoring powers of 2 from num−1
    T d = static_cast<T>(num - 1);
    int r = 0;
    while (d % 2 == 0) {
        ++r;
        d = static_cast<T>(d/2);
    }
    HPBC_ASSERT2(r > 0);

    for (size_t i=0; i<total_bases; ++i) {
        if (bases[i] == 0)
            continue;
        auto result = mont.pow(mont.convertIn(bases[i]), d);
        auto canonicalResult = mont.getCanonicalForm(result);
        if (canonicalResult == unity || canonicalResult == negativeOne)
            continue;
        for (int j=0; j<r-1; ++j) {
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



// So far I've verified the correctness of is_prime() for all odd input values
// from 0 to 34359738368 (== 2^35), looking for any mismatch compared to
// is_prime_wheel210() in  test/is_prime_wheel210.h
static constexpr uint64_t MAX_VERIFIED_NUMBER_MONTGOMERY_ISPRIME =
                                                          UINT64_C(34359738368);

template <typename MontType>
bool is_prime_miller_rabin(const MontType& mont)
{
    using T = typename MontType::T_type;
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(ma::ma_numeric_limits<T>::digits <= 128, "");

    T num = mont.getModulus();

    HPBC_PRECONDITION2(num >= 3);
    HPBC_ASSERT2(num % 2 == 1);

    if (num > std::numeric_limits<uint64_t>::max()) {
        // Other algorithms for primality testing are likely to be far more
        // suitable for large 128 bit numbers than Miller-Rabin, but we can
        // still use Miller-Rabin.
        //
        // Using Miller-Rabin for large 128 bit numbers, we don't have a known
        // deterministic test.  Instead we'll depend upon setting up an
        // incredibly low probability that the trials will ever declare that a
        // composite number is prime.  We can do this by using a very large
        // number of bases.  According to
        // https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Miller%E2%80%93Rabin_test
        // if k is the number of bases, then the probability is 4^(-k) that
        // the trials will declare that a particular composite number is prime
        // (it will never declare a prime is composite).  So with 64 bases, if
        // we tested every number from 1 to 2^128, we might expect very roughly
        // about 1 false declaration of prime.  With an additional 64 bases (for
        // a total of 128), the probability that this single expected false
        // declaration would still remain would be very roughly 1/2^128, or 1
        // chance in 340282366920938463463374607431768211456.  See the Wikipedia
        // link's discussion of the conditional probability and bayes theorem
        // for some info on why I wouldn't be surprised if the true probability
        // is 100x more likely.  Yet even assuming it's 100000 times more
        // likely, it would have odds of 1 in a decillion (roughly one over the
        // number of grains of sand on earth, squared) that any number exists
        // between 1 and 1^128 for which it would deliver a wrong answer.
        // We'll use all primes for the bases, simply following the example of
        // http://oeis.org/A014233 (I don't know why it uses primes).
        static const uint16_t bases[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
            31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
            103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
            173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
            241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
            317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397,
            401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
            479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569,
            571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643,
            647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719 };
        constexpr size_t total_bases = 128;
        static_assert(sizeof(bases)/sizeof(bases[0]) == total_bases, "");
        auto basesT = std::unique_ptr<T[]>(new T[total_bases]);
        for (size_t i=0; i<total_bases; ++i)
            basesT[i] = static_cast<T>(bases[i]);
        return is_prime_mr_trials<T, MontType>(basesT.get(), total_bases, mont);

    } else if (num > std::numeric_limits<uint32_t>::max()) {
        uint64_t num64 = static_cast<uint64_t>(num);

        // All bases here are smaller than num (we know num is > UINT32_MAX), so
        // we don't need to reduce any bases mod num, prior to trials.
        if (num64 <= MAX_VERIFIED_NUMBER_MONTGOMERY_ISPRIME) {
            static_assert(MAX_VERIFIED_NUMBER_MONTGOMERY_ISPRIME <
                                                    UINT64_C(154639673381), "");
            // According to http://miller-rabin.appspot.com, the 3 bases below
            // are sufficient for num up to 154639673381.  We'll use this only
            // up to the max I've verified so far, which is given by
            // MAX_VERIFIED_NUMBER_MONTGOMERY_ISPRIME.
            const T bases[] = { 15, static_cast<T>(176006322),
                                static_cast<T>(4221622697) };
            size_t total_bases = sizeof(bases)/sizeof(bases[0]);
            return is_prime_mr_trials<T, MontType>(bases, total_bases, mont);
        } else {
            // According to http://miller-rabin.appspot.com, the first 7 bases
            // below are sufficient for num up to 2^64.  Practically speaking, I
            // can't verify this.  Note that an alternative (which I also can't
            // verify) is to use the smallest 12 primes (2,3,5,7,11, etc) as
            // bases for num < 2^64; for details see http://oeis.org/A014233
            //
            // For now I've combined these two sets of bases to make it extra
            // likely that they'll be sufficient for any num < 2^64.  We'll use
            // these bases for any number range I haven't brute force verified.
            // Combining the bases is almost certainly overkill.
            //
            // **TODO - use a less paranoid set of bases here**
            //
            const T bases[] = { 2, static_cast<T>(325), static_cast<T>(9375),
                                static_cast<T>(28178), static_cast<T>(450775),
                                static_cast<T>(9780504),
                                static_cast<T>(1795265022), 3, 5, 7, 11, 13, 17,
                                19, 23, 27, 31, 37 };
            size_t total_bases = sizeof(bases)/sizeof(bases[0]);
            return is_prime_mr_trials<T, MontType>(bases, total_bases, mont);
        }
    } else {   // num < 2^32
        // All initialization values for bases, and ranges of validity come from
        // http://miller-rabin.appspot.com/

        uint32_t num32 = static_cast<uint32_t>(num);
        T bases[3];
        size_t total_bases;

#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  error HURCHALLA_TARGET_BIT_WIDTH must be defined
#endif
#if HURCHALLA_TARGET_BIT_WIDTH >= 64
        if (num32 >= UINT32_C(1050535501)) {
            total_bases = 3;
            // mod isn't needed here, since num is greater than all these bases.
            bases[0] = 2;
            bases[1] = 7;
            bases[2] = 61;    
        } else if (num32 >= UINT32_C(19471033)) {
            total_bases = 2;
            bases[0] = static_cast<T>(UINT64_C(336781006125) % num32);
            bases[1] = static_cast<T>(UINT64_C(9639812373923155) % num32);
        } else if (num32 >= UINT32_C(341531)) {
            total_bases = 2;
            // mod isn't needed here, since num is greater than all these bases.
            bases[0] = 2;
            bases[1] = static_cast<T>(299417);
        } else {   // From the precondition we know num is still >= 3.
            total_bases = 1;
            bases[0] = static_cast<T>(UINT64_C(9345883071009581737) % num32);
        }
#else
        if (num32 >= UINT32_C(360018361) {
            total_bases = 3;
            // These bases should cover up to 4,759,123,141 (all uint32_t).
            // mod isn't needed here, since num is greater than all these bases.
            bases[0] = 2;
            bases[1] = 7;
            bases[2] = 61;
        } else if (num32 >= UINT32_C(49141) {
            total_bases = 2;
            uint32_t tmp = UINT32_C(1143370);
            if (tmp >= num32)
                bases[0] = static_cast<T>(tmp % num32);
            else
                bases[0] = static_cast<T>(tmp);
            bases[1] = static_cast<T>(UINT32_C(2350307676) % num32);
        } else {  // From the precondition we know num is still >= 3.
            total_bases = 1;
            bases[0] = static_cast<T>(UINT32_C(921211727) % num32);
        }
#endif
        return is_prime_mr_trials<T, MontType>(bases, total_bases, mont);
    }
}


}}  // end namespace

#endif
