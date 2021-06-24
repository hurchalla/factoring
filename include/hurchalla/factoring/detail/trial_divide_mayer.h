// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IS_TRIAL_DIVIDE_MAYER_H_INCLUDED
#define HURCHALLA_FACTORING_IS_TRIAL_DIVIDE_MAYER_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/inverse_mod_R.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstdint>
#include <type_traits>
#include <limits>

namespace hurchalla { namespace detail {


#if defined(HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE)
#  error "HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE must not be predefined"
#endif
#if defined(HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE) || \
    !defined(HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE)
#  define HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE true
#else
#  define HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE false
#endif

// This algorithm is "ALGORITHM A: IS_DIV_A" from section 3. A Fast Multiword
// Divisibility Test, in the paper "Efficient long division via Montgomery
// multiply" by Ernst W. Mayer (https://arxiv.org/abs/1303.0328).  The form of
// the algorithm we use here is the special case where the "unsigned b-bit
// integer array x[n]" (see the paper) has exactly one element (i.e. n == 1).
// Extra note on the algorithm: the special case we use is quite closely related
// to the very differently derived algorithm "Test for Zero Remainder after
// Division by a Constant" of Section 10-17 of Hacker's Delight 2nd edition by
// Henry Warren.
//
// Since we use only a special case of Mayer's algorithm, we can present an easy
// explanation/proof for it:
// Let T be an unsigned integer type, and let the value R = 2^(bit_width_of_T),
// where '^' is shorthand for exponentiation (not xor).
// For example, if T is uint64_t then R = 2^64.
// Given a type T value x and a type T odd value n, we know n is coprime to R
// because R is a power of 2 and n is odd.  Therefore the inverse of n (mod R)
// exists, which we will name inv_n.  Since x is type T, we know x < R, and
// since n is type T and odd, n > 0.  Thus x < n*R.  This fulfills all of the
// requirements to use the alternate REDC algorithm, which is described at
// https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/README_REDC.md
//
// Using the linked REDC algorithm, REDC(x, R, n, inv_n) outputs a value
// t ≡ x * R^(-1)  (mod n).  Multiplying both sides by R,  t*R ≡ x (mod n).
// Therefore if t ≡ 0 (mod n), then x ≡ 0 (mod n).
// If instead we assume x ≡ 0 (mod n), then  t ≡ 0 * R^(-1) ≡ 0  (mod n);  by
// the converse, if t !≡ 0 (mod n), then x !≡ 0 (mod n).
// Therefore, x ≡ 0 (mod n) iff t ≡ 0 (mod n).  Since REDC gives us  0 <= t < n,
// x ≡ 0 (mod n) iff t == 0.  This finishes the proof.  The remainder of the
// explanation fills in the details of using the REDC algorithm for our needs:
//
// Directly implementing the alternate REDC algorithm, we get
//    T x_lo = x;
//    T x_hi = 0;
//    T m = x_lo * inv_n;
//    T mn_hi = multiply_to_hiword(m, n);  // ignore the low word of the product
//    T t = x_hi - mn_hi;
//    if (x_hi < mn_hi)
//        t += n;
// Since we know x ≡ 0 (mod n) iff t == 0, we return "x is divisible by n" via
//    return (t == 0);
//
// Since we know x_hi == 0, we can simplify the implementation to
//    T m = x * inv_n;
//    T mn_hi = multiply_to_hiword(m, n);
//    T t;
//    if (mn_hi != 0)   // simplified using the fact that x_hi == 0
//        t = n - mn_hi;
//    else
//        t = -mn_hi;   // to get here, mn_hi must == 0; thus this sets t = 0.
//    return (t == 0);
//
// We know mn_hi < n, due to the proof of Assertion #1 in REDC_non_finalized()
// https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/low_level_api/detail/platform_specific/impl_REDC.h
// so we therefore know  n - mn_hi != 0.  So we can further simplify:
//    T m = x * inv_n;
//    T mn_hi = multiply_to_hiword(m, n);
//    if (mn_hi != 0)
//        return false;
//    else
//        return true;


// trial_divide_mayer():  returns true if n divides x, otherwise returns false.
// If n divides x, then the quotient is placed in div_result.  If n
// does not divide x, then the value of div_result is unspecified.
//
// Precondition: n must be odd.

template <typename T>
HURCHALLA_FLATTEN
typename std::enable_if<(HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE ||
                     ut_numeric_limits<T>::digits > HURCHALLA_TARGET_BIT_WIDTH),
                 bool>::type
// Use the special algorithm described above.
// Note 'x' is the dividend and 'n' is the divisor.
trial_divide_mayer(T& div_result, T x, T n)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(n%2 == 1);  // required for calling inverse_mod_R

    using P = typename safely_promote_unsigned<T>::type;

    T inv_n = inverse_mod_R(n);
    T m = static_cast<T>(static_cast<P>(x) * inv_n);
    T mn_lo;
    T mn_hi = unsigned_multiply_to_hilo_product(mn_lo, m, n);
    div_result = m;
    return (mn_hi == 0);
}

template <typename T>
typename std::enable_if<!(HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE ||
                     ut_numeric_limits<T>::digits > HURCHALLA_TARGET_BIT_WIDTH),
                 bool>::type
// Don't use the algorithm described above.  Instead, just use plain standard
// integer division.
// Note 'x' is the dividend and 'n' is the divisor.
trial_divide_mayer(T& div_result, T x, T n)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(n > 0);  // disallow division by 0
    HPBC_PRECONDITION2(n%2 == 1); // for consistency with trial_divide_mayer
                                 // above, though this version doesn't need odds
    using P = typename safely_promote_unsigned<T>::type;
    div_result = static_cast<T>(static_cast<P>(x) / n);
    // test whether the remainder (x - n*div_result) equals 0
    return (x == n*div_result);
}

#undef HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE


}}  // end namespace

#endif
