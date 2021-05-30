// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_FACTORIZE_WHEEL210_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_WHEEL210_H_INCLUDED


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
#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4702)
#endif

namespace hurchalla { namespace detail {


// This algorithm is Wheel Factorization
// https://en.wikipedia.org/wiki/Wheel_factorization

// Preconditions for factorize_wheel210():
//   Requires x >= 2.
// Postconditions:
//   The return value is an output iterator to the position one past the last
//   factor that the function wrote to the destination range (iterated by the
//   function's parameter 'iter').  The destination range will consist of all
//   the prime factors of x; the function will always completely factor x,
//   though it may take an extremely long time if x is large.  The destination
//   range will consist of the single element x if the function determined that
//   x was prime.
//
// factorize_wheel210 will trial all potential prime factors less than or equal
// to sqrt(x).  It may also trial a few factors that are larger than sqrt(x).

// For discussion purposes inside the function, let the theoretical constant
// R == 1 << ut_numeric_limits<T>::digits.  For example, for a type T that is
// uint16_t, R would equal 65536 and sqrtR would equal 256.

template <class OutputIt, typename T>
OutputIt factorize_wheel210(OutputIt iter, T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations

    using std::size_t;
    using std::uint8_t;
    // avoid integral promotion hassles
    using P = typename safely_promote_unsigned<T>::type;
    // the maximum possible number of factors for a type T variable occurs when
    // all factors are 2.  Thus bitsT is >= to the number of factors of x.
    static constexpr int bitsT = ut_numeric_limits<T>::digits;
    static_assert(bitsT % 2 == 0);
    static constexpr T sqrtR = static_cast<T>(1) << (bitsT / 2);

    P q = static_cast<P>(x);
    HPBC_ASSERT2(q > 1);
    // We test divisors to 13, to cover all possible factors for type uint8_t,
    // meaning that uint8_t never needs to use the wheel (it might result in
    // overflow if it did).  Note this is a bit overkill for larger types (which
    // with a modified wheel below would only need to test up to 7 here), but it
    // makes essentially no difference on performance.
    while (q%2 == 0) {
        q = q/2;
        *iter++ = static_cast<T>(2);
        if (q == 1)  // we completely factored x
            return iter;
    }
    while (q%3 == 0) {
        q = q/3;
        *iter++ = static_cast<T>(3);
        if (q == 1)  // we completely factored x
            return iter;
    }
    while (q%5 == 0) {
        q = q/5;
        *iter++ = static_cast<T>(5);
        if (q == 1)  // we completely factored x
            return iter;
    }
    while (q%7 == 0) {
        q = q/7;
        *iter++ = static_cast<T>(7);
        if (q == 1)  // we completely factored x
            return iter;
    }
    while (q%11 == 0) {
        q = q/11;
        *iter++ = static_cast<T>(11);
        if (q == 1)  // we completely factored x
            return iter;
    }
    while (q%13 == 0) {
        q = q/13;
        *iter++ = static_cast<T>(13);
        if (q == 1)  // we completely factored x
            return iter;
    }
    HPBC_ASSERT2(q > 1);
    if constexpr (std::is_same<T, uint8_t>::value) {
        // if x is type uint8_t, we would have tested all potential factors
        // less than sqrtR above, and therefore q must be prime.
        *iter++ = static_cast<T>(q);
        return iter;
    }

    // The wheel spans the 210 number range [17, 227), skipping all multiples
    // of 2,3,5,7 within that range
    static constexpr uint8_t wheel[] = { 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
        59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 121, 127,
        131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187,
        191, 193, 197, 199, 209, 211, 221, 223 };
    static constexpr size_t wheel_len = sizeof(wheel)/sizeof(wheel[0]);
    static constexpr uint8_t cycle_len = 210;

    HPBC_ASSERT2(q > 1);
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
            *iter++ = static_cast<T>(q);
            return iter;
        }
        HPBC_ASSERT2(maybe_factor < sqrtR);
        // This inner loop will usually trial a few maybe_factor(s) that are
        // greater than sqrtR or sqrt(q), which is unneeded but harmless.
        for (size_t i=0; i < wheel_len; ++i) {
            // Assert that  start + wheel[i]  will never overflow.  Let
            // S = ma_numeric_limits<P>::max().  Overflow would mean that
            // start + wheel[i] > S, or equivalently  S - wheel[i] < start.  We
            // asserted above that  start + wheel[0] == maybe_factor < sqrtR.
            // So overflow would mean  S - wheel[i] < start < sqrtR - wheel[0],
            // which would mean  S - sqrtR < wheel[i] - wheel[0] < cycle_len ==
            // 210.  But S - sqrtR < 210 is impossible because P is always at
            // at least as large as uint16_t (due to promotion rules it's based
            // upon), and thus S >= 65535, and sqrtR is always <= sqrt(S+1).
            static_assert(cycle_len == 210); //support the preceding comment
            HPBC_ASSERT2(start <= ut_numeric_limits<P>::max() - wheel[i]);
            maybe_factor = start + wheel[i];
            P div_result;
            HPBC_ASSERT2(q > 1);
            // test whether  maybe_factor divides q  without any remainder.
            while (trial_divide_mayer(div_result, q, maybe_factor)) {
                *iter++ = static_cast<T>(maybe_factor);
                q = div_result;
                if (q == 1)  // we completely factored x
                    return iter;
            }
        }
    }
}


}}  // end namespace

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic pop
#endif

#endif
