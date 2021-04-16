// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_FACTORIZE_WHEEL210_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_WHEEL210_H_INCLUDED


#include "hurchalla/factoring/detail/factorize_trialdivision.h"
#include "hurchalla/factoring/detail/trial_divide_mayer.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <numeric>
#include <functional>

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunsafe-loop-optimizations"
#endif

namespace hurchalla { namespace detail {


// This algorithm is Wheel factorization
// https://en.wikipedia.org/wiki/Wheel_factorization
// but with testing for potential factors only up to (and including) max_factor.

// Postconditions:
// (note these are the same as for factorize_trialdivision())
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
//   factored.  q will never be set to zero (or a value < 0).  If max_factor
//   is defaulted or set to a value >= sqrt(x), then the function will always
//   completely factor x and set q = 1.
//
// factorize_wheel210 guarantees it will trial all potential prime factors
// less than or equal to max_factor.  It will usually also trial a few factors
// that are larger than max_factor.

// for discussion purposes inside the function, let the theoretical constant
// R == 1 << ut_numeric_limits<T>::digits.  For example, for a type T that is
// uint16_t, R would equal 65536 and sqrtR would equal 256.

// the template-template param TTD should be either TrialDivisionWarren or
// TrialDivisionMayer.
template <template<class,int> class TTD, class OutputIt, typename T>
OutputIt factorize_wheel210(OutputIt iter, T& q, T x,
                                     T max_factor = ut_numeric_limits<T>::max())
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations

    // the maximum possible number of factors for a type T variable occurs when
    // all factors are 2.  Thus bitsT is >= to the number of factors of x.
    static constexpr int bitsT = ut_numeric_limits<T>::digits;
    static_assert(bitsT % 2 == 0, "");
    static constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);
    if (max_factor >= sqrtR)
        max_factor = sqrtR - 1;

    // factor out all primes < 256
    constexpr int SIZE=54;   // there are 54 primes < 256
    T next_prime;  // ignored for now
    iter = factorize_trialdivision<TTD, SIZE>(iter, q, next_prime, x, SIZE);
    HPBC_ASSERT2(q >= 1);  // factorize_trialdivision guarantees this
    if (q == 1)   // if factorize_trialdivision completely factored x
        return iter;

    using std::size_t;
    using std::uint8_t;
    // avoid integral promotion hassles
    using P = typename safely_promote_unsigned<T>::type;
    P q2 = q;
    HPBC_ASSERT2(q2 > 1);

    // The wheel spans the 210 number range [47, 257), skipping all multiples
    // of 2,3,5,7 within that range
    static constexpr uint8_t wheel[] = { 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97, 101, 103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151,
        157, 163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209, 211,
        221, 223, 227, 229, 233, 239, 241, 247, 251, 253 };
    static constexpr size_t wheel_len = sizeof(wheel)/sizeof(wheel[0]);
    static constexpr uint8_t cycle_len = 210;

    static constexpr P cycle_start = 210;
    // 257 is the first prime > 256, so we need to resume trial division there
    static_assert(cycle_start + wheel[0] == 257, "");

    // Use the wheel to skip division by any multiple of 2,3,5,7 in the loop
    // below - we already tested 2,3,5,7 via factorize_trialdivision().
    P maybe_factor;
    for (P start=cycle_start; ; start=static_cast<P>(start + cycle_len)) {
        maybe_factor = start + wheel[0];
        if (maybe_factor > max_factor)
            break;
        // Since the above clause leaves the loop, we know at this point
        // that  maybe_factor<=max_factor.  And since  max_factor<sqrtR,  we
        // know maybe_factor<sqrtR, and thus  maybe_factor*maybe_factor<R,
        // which means (maybe_factor*maybe_factor) doesn't overflow.
        HPBC_ASSERT2(maybe_factor < sqrtR);
        if (maybe_factor * maybe_factor > q2)
            break;   // no factors ever exist > sqrt(q)
        // This inner loop will usually trial a few maybe_factor(s) that are
        // greater than max_factor or sqrt(q2), which is unneeded but harmless.
        for (size_t i=0; i < wheel_len; ++i) {
            // Assert that  start + wheel[i]  will never overflow.  Let
            // S = ma_numeric_limits<P>::max().  Overflow would mean that
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
            HPBC_ASSERT2(q2 > 1);
            // test whether  maybe_factor divides q2  without any remainder.
            while (trial_divide_mayer(div_result, q2, maybe_factor)) {
                *iter++ = static_cast<T>(maybe_factor);
                q2 = div_result;
                if (q2 == 1) {  // we completely factored x
                    q = static_cast<T>(q2);
                    return iter;
                }
            }
        }
    }
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

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic pop
#endif

#endif
