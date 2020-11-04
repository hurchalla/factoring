// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_WHEEL_FACTORIZATION210_H_INCLUDED
#define HURCHALLA_FACTORING_WHEEL_FACTORIZATION210_H_INCLUDED


#include "hurchalla/factoring/detail/small_trial_division256.h"
#include "hurchalla/factoring/detail/trial_divide.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <numeric>
#include <functional>

namespace hurchalla { namespace factoring {


// This algorithm is Wheel factorization
// https://en.wikipedia.org/wiki/Wheel_factorization
// but with testing for potential factors only up to (and including) max_factor.

// Postconditions:
// (note these are the same as for small_trial_division256())
// 1) The return value is an output iterator to the position one past the last
//   factor that the function wrote to the destination range (iterated by the
//   function's parameter 'iter').  The destination range consists of all the
//   prime factors that the function was able to find for x.  This range will
//   be empty if the function could not find any factors for x and could not
//   determine whether x was prime.  The range will consist of the single
//   element x if the function determined that x was prime.
// 2) q will be set to the quotient of x divided by all the elements written to
//   the destination range (iterated by the function's parameter 'iter').  If
//   nothing was written to the range (if it's empty), then q will be set to x.
//   There are specific details that naturally follow from these facts and from
//   Postcondition 1: if q gets set to 1, then this indicates the function was
//   able to completely factor x and the destination range consists of all the
//   factors.  If q > 1, then this indicates the function was not able to
//   completely factorize x, and q represents the value still remaining to be
//   factored.  q will never be set to zero (or a value < 0).
//
// wheel_factorization210 guarantees it will trial all potential prime factors
// less than or equal to max_factor.


template <class OutputIt, typename T>
OutputIt wheel_factorization210(OutputIt iter, T& q, T x,
      T max_factor = hurchalla::modular_arithmetic::ma_numeric_limits<T>::max())
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations

    // the maximum possible number of factors for a type T variable occurs when
    // all factors are 2.  Thus bitsT is >= to the number of factors of x.
    static constexpr int bitsT = ma::ma_numeric_limits<T>::digits;
    static_assert(bitsT % 2 == 0, "");
    static constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);
    if (max_factor >= sqrtR)
        max_factor = sqrtR - 1;

    // factor out all primes < 256
    iter = small_trial_division256(iter, q, x);
    HPBC_ASSERT2(q >= 1);  // small_trial_division256 guarantees this
    if (q == 1)   // if small_trial_division256 completely factored x
        return iter;

    using std::size_t;
    using std::uint8_t;
    // avoid integral promotion hassles
    using P = typename montgomery_arithmetic::safely_promote_unsigned<T>::type;
    P q2 = q;
    HPBC_ASSERT2(q2 > 1);

    // The wheel spans the 210 number range [47, 257), skipping all multiples
    // of 2,3,5,7 within that range
    static constexpr uint8_t wheel[] = { 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151,
        157, 163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209, 211,
        221, 223, 227, 229, 233, 239, 241, 247, 251, 253 };
    static constexpr size_t wheel_len = sizeof(wheel)/sizeof(wheel[0]);
    static const uint8_t cycle_len = 210;

    P cycle_start = 210;
    // 257 is the first prime > 256, so we need to resume trial division there
    HPBC_ASSERT2(cycle_start + wheel[0] == 257);

    // Use the wheel to skip division by any multiple of 2,3,5,7 in the loop
    // below - we already tested 2,3,5,7 via small_trial_division256().
    P maybe_factor;
    for (P start=cycle_start; ; start=static_cast<P>(start + cycle_len)) {
        for (size_t i=0; i < wheel_len; ++i) {
            HPBC_ASSERT2(q2 > 1);
            maybe_factor = start + wheel[i];
            if (maybe_factor > max_factor)
                goto endloop;
            // Since the above clause leaves the loop, we know at this point
            // that  maybe_factor<=max_factor.  And since  max_factor<sqrtR,  we
            // know maybe_factor<sqrtR, and thus  maybe_factor*maybe_factor<R,
            // which means (maybe_factor*maybe_factor) doesn't overflow.
            HPBC_ASSERT2(maybe_factor < sqrtR);
            if (maybe_factor * maybe_factor > q2)
                goto endloop;   // no factors ever exist > sqrt(q)
            P div_result;
            // test whether  maybe_factor divides q2  without any remainder.
            while (trial_divide(div_result, q2, maybe_factor)) {
                *iter++ = static_cast<T>(maybe_factor);
                q2 = div_result;
                if (q2 == 1) {  // we completely factored x
                    q = static_cast<T>(q2);
                    return iter;
                }
            }
        }
    }
endloop:
    if (maybe_factor >= sqrtR || maybe_factor * maybe_factor > q2) {
        // Since R > q2, we know sqrtR > sqrt(q2).  Therefore
        // maybe_factor >= sqrtR  implies  maybe_factor > sqrt(q2).
        // And  maybe_factor * maybe_factor > q2  obviously implies
        // the same.  So inside this clause,  maybe_factor > sqrt(q2).
        // This means we have tried all potential prime factors
        // less than or equal to sqrt(q2), so q2 must be prime.
        *iter++ = static_cast<T>(q2);
        q2 = 1;
    }
    q = static_cast<T>(q2);
    HPBC_POSTCONDITION2(q > 0);
    return iter;
}


}}  // end namespace

#endif
