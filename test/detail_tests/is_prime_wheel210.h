// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IS_PRIME_WHEEL210_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_WHEEL210_H_INCLUDED


#include "hurchalla/factoring/detail/trial_divide_mayer.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>

namespace hurchalla { namespace detail {


// This algorithm is an adaptation of Wheel factorization, where we return false
// upon the first factor found, or we return true if we determine that the input
// is prime.  This function will always be able to successfully determine
// primality.
// For more info see:
// https://en.wikipedia.org/wiki/Wheel_factorization
// Wheel factorization is a slight optimization (~2x) over is_prime_bruteforce()
// but it's only a constant time improvement and remains a brute force approach.

template <typename T>
bool is_prime_wheel210(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    using std::size_t;
    using std::uint8_t;
    // avoid integral promotion hassles
    using P = typename safely_promote_unsigned<T>::type;
    static constexpr int bitsT = ut_numeric_limits<T>::digits;
    static_assert(bitsT % 2 == 0, "");
    static constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);

    P q = static_cast<P>(x);
    if (q < 2) return false;
    // We test divisors to 13, to cover all possible factors for type uint8_t,
    // meaning that uint8_t never needs to use the wheel (it might result in
    // overflow if it did).  Note this is a bit overkill for larger types (which
    // with a modified wheel below would only need to test up to 7 here), but it
    // makes essentially no difference on performance.
    if (q%2 == 0) return (q==2);
    if (q%3 == 0) return (q==3);
    if (q%5 == 0) return (q==5);
    if (q%7 == 0) return (q==7);
    if (q%11 == 0) return (q==11);
    if (q%13 == 0) return (q==13);
    // We know by the fact that we got to this point that x has no
    // factors <= 13.  Thus if x is type uint8_t, it has no prime factors that
    // are < sqrtR, and therefore x must be prime.
    if constexpr (std::is_same<T, uint8_t>::value)
        return true;

    // The wheel spans the 210 number range [17, 227), skipping all multiples
    // of 2,3,5,7 within that range
    static constexpr uint8_t wheel[] = { 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
        59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 121, 127,
        131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187,
        191, 193, 197, 199, 209, 211, 221, 223 };
    static constexpr size_t wheel_len = sizeof(wheel)/sizeof(wheel[0]);
    static constexpr uint8_t cycle_len = 210;

    // Use the wheel to skip division by any multiple of 2,3,5,7 in the loop
    // below - we already tested 2,3,5,7 and found they were not factors.
    // Note the wheel pattern cycles every 2*3*5*7 == 210 numbers.
    P maybe_factor;
    for (P start=0; ; start=static_cast<P>(start + cycle_len)) {
        maybe_factor = start + wheel[0];
        if (maybe_factor >= sqrtR || maybe_factor*maybe_factor > q) {
            // note: since (maybe_factor*maybe_factor) is evaluated only when
            // maybe_factor < sqrtR, it will never overflow.
            // Since R > q, we know sqrtR > sqrt(q).  Therefore
            // maybe_factor >= sqrtR  implies  maybe_factor > sqrt(q).
            // And  maybe_factor * maybe_factor > q  obviously implies
            // the same.  So inside this clause,  maybe_factor > sqrt(q).
            // This means we have tried all potential prime factors
            // less than or equal to sqrt(q), so q must be prime.
            return true;
        }
        HPBC_ASSERT2(maybe_factor < sqrtR);
        // This inner loop will usually trial a few maybe_factor(s) that are
        // greater than sqrtR or sqrt(q), which is unneeded but harmless.
        for (size_t i=0; i < wheel_len; ++i) {
            // Assert that  start + wheel[i]  will never overflow.  Let
            // S = ut_numeric_limits<P>::max().  Overflow would mean that
            // start + wheel[i] > S, or equivalently  S - wheel[i] < start.  We
            // asserted above that  start + wheel[0] == maybe_factor < sqrtR.
            // So overflow woule mean  S - wheel[i] < start < sqrtR - wheel[0],
            // which would mean  S - sqrtR < wheel[i] - wheel[0] < cycle_len ==
            // 210.  But S - sqrtR < 210 is impossible because P is always at
            // at least as large as uint16_t (due to promotion rules it's based
            // upon), and thus S >= 65535, and sqrtR is always <= sqrt(S+1).
            static_assert(cycle_len == 210, ""); //support the preceding comment
            HPBC_ASSERT2(start <= ut_numeric_limits<P>::max() - wheel[i]);
            maybe_factor = start + wheel[i];
            P div_result;
            // test whether  maybe_factor divides q  without any remainder.
            if (trial_divide_mayer(div_result, q, maybe_factor))
                return false;
        }
    }
}


}}  // end namespace

#endif
