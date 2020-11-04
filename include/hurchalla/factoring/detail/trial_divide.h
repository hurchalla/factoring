// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_IS_TRIAL_DIVIDE_H_INCLUDED
#define HURCHALLA_FACTORING_IS_TRIAL_DIVIDE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/montgomery_arithmetic/detail/inverse_mod_r.h"
#include "hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h"
#include <cstdint>
#include <type_traits>
#include <limits>

namespace hurchalla { namespace factoring {


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
// Let T be an unsigned integer type, and let the value R = 2^(bit_width_of_T).
// For example, if T is uint64_t then R = 2^64.
// Given a type T value x and a type T odd value n, we know n is coprime to R
// because R is a power of 2 and n is odd.  Therefore the inverse of n (mod R)
// exists, which we will name inv_n.  Since x is type T, we know x < R, and
// since n is type T and odd, n > 0.  Thus x < n*R.  This fulfills all of the
// requirements to use the alternate REDC algorithm, which is described at
// https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/platform_specific/README_REDC.md
//
// Using the linked REDC algorithm, REDC(x, R, n, inv_n) outputs a value
// t ≡ x * R^(-1)  (mod n).  Multiplying both sides by R,  t*R ≡ x (mod n).
// Therefore if t ≡ 0 (mod n), then x ≡ 0 (mod n).
// And if we assume x ≡ 0 (mod n), then  t ≡ 0 * R^(-1) ≡ 0  (mod n);  by the
// converse, if t !≡ 0 (mod n), then x !≡ 0 (mod n).
// Therefore, x ≡ 0 (mod n) iff t ≡ 0 (mod n).  Since REDC gives us  0 <= t < n,
// x ≡ 0 (mod n) iff t == 0.  This finishes the proof.  The remainder of the
// explanation fills in the details of using the REDC algorithm for our needs:
//
// Directly implementing the alternate REDC algorithm, we get
//    T x_lo = x;   T x_hi = 0;
//    T m = x_lo*inv_n;
//    T mn_hi = multiply_to_hiword(m, n);  // ignore the low word of the product
//    T t = x_hi - mn_hi;
//    if (x_hi < mn_hi)
//        t += n;
// Since we know x ≡ 0 (mod n) iff t == 0, we return (x is divisible by n) via
//    return (t == 0);
//
// Since we know x_hi == 0, we can simplify the implementation to
//    T m = x*inv_n;
//    T mn_hi = multiply_to_hiword(m, n);
//    T t;
//    if (mn_hi != 0)
//        t = n - mn_hi;
//    else
//        t = -mn_hi;  // to get here, mn_hi must == 0; thus this sets t = 0.
//    return (t == 0);
//
// We know mn_hi < n, due to the proof of Assertion #1 in REDC_non_finalized()
// https://github.com/hurchalla/modular_arithmetic/blob/master/montgomery_arithmetic/include/hurchalla/montgomery_arithmetic/detail/platform_specific/Redc.h
// so we therefore know  n - mn_hi != 0.  So we can further simplify:
//    T m = x*inv_n;
//    T mn_hi = multiply_to_hiword(m, n);
//    if (mn_hi != 0)
//        return false;
//    else
//        return true;

// trial_divide:  returns true if n divides x, otherwise returns false.
// If n divides x, then the quotient is placed in div_result.  If n
// does not divide x, then the value of div_result is undefined.

template <typename T>
typename std::enable_if<(HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE ||
                   hurchalla::modular_arithmetic::ma_numeric_limits<T>::digits >
                     HURCHALLA_TARGET_BIT_WIDTH), bool>::type
// Use our special algorithm described above.
// Note 'x' is the numerator and 'n' is the denominator.
trial_divide(T& div_result, T x, T n)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(n%2 == 1);  // required for calling inverse_mod_r

    namespace mont = hurchalla::montgomery_arithmetic;
    using P = typename mont::safely_promote_unsigned<T>::type;

    T inv_n = mont::inverse_mod_r(n);
    T m = static_cast<T>(static_cast<P>(x) * inv_n);
    T mn_lo;
    T mn_hi = mont::unsigned_multiply_to_hilo_product(&mn_lo, m, n);
    div_result = m;
    return (mn_hi == 0);
}

template <typename T>
typename std::enable_if<!(HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE ||
                   hurchalla::modular_arithmetic::ma_numeric_limits<T>::digits >
                     HURCHALLA_TARGET_BIT_WIDTH), bool>::type
// Don't use the algorithm described above.  Instead, just use plain standard
// integer division.
// Note 'x' is the numerator and 'n' is the denominator.
trial_divide(T& div_result, T x, T n)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(n > 0);  // disallow division by 0

    namespace mont = hurchalla::montgomery_arithmetic;
    using P = typename mont::safely_promote_unsigned<T>::type;
    div_result = static_cast<T>(static_cast<P>(x) / n);
    // test whether the remainder (x - n*div_result) equals 0
    return (x == n*div_result);
}

#undef HURCHALLA_USE_TRIAL_DIVIDE_VIA_INVERSE




#ifdef HURCHALLA_ENABLE_EXPERIMENTAL_TRIAL_DIVIDE_UINT32_VIA_DOUBLE
// This version of the function is unproven to be correct (or incorrect).  It
// divides two uint32_t variables by casting to double and performing the
// division with doubles.  In theory the compiler should use SIMD instructions
// instructions for it (such as SSE) if available, which for x86_64 would
// provide a divide with latency somewhere around 15-22 cycles (depending on
// architecture) and throughput somewhere around 4-10 cycles per instruction
// (depending on architecture).  According to uops.info this could provide
// similar latency to integer division for uint32_t, but either similar or up to
// 1.5x better throughput.  I expect it to provide similar or slightly better
// performance than the generic version of trial_divide above that uses the REDC
// technique.  Note that if this function is proven correct, it opens up the
// ability for a SIMD version of this function to be created, and used by
// clients, which could offer 2-4x the throughput of this version (for 128 bit
// to 256 bit vectors).  Some quesions about performance that I don't have
// answers to on this function's performance: depending on compiler flags could
// there be slowdown due to the compiler handling the possibility of divide by
// zero, or handling possibility of other floating point exceptions, or NaNs or
// infinities?  In theory, the compiler should at least see NaNs/inifinites are
// impossible in this function.
// This function is experimental because I have not proven that calculating
// the quotient of two uint32_t variables via a double (floating point) division
// produces the same result as a straightforward divide of the two uint32_t
// variables.  I suspect that it is always a correct operation, but I have
// neither made a proof nor found example input values that disprove it.
// A further complication is that the compiler (if fast-math optimization flag
// is enabled) might in some circumstances first calculate the double value
// for the reciprocal of n, and then multiply x by that reciprocal to get the
// quotient - this introduces additional rounding error, but I suspect that
// getting the results via the reciprocal will always produce a correct answer
// too, though it becomes an even more uncertain guess.  There are three ways I
// can think to prevent the compiler from using the reciprocal, but none of them
// are good ideas.  First way would be to use inline asm, but for x86 I would
// want to use SSE asm instructions, and if the calling code happens to use AVX
// there is the potential for the upper bits to be dirty which means a
// performance penalty:
// https://stackoverflow.com/questions/41303780/why-is-this-sse-code-6-times-slower-without-vzeroupper-on-skylake
// and I'd prefer not to use AVX instructions in case the processor doesn't
// support AVX.  Regardless, even with a perfect solution to that, inline asm
// would require different code for ARM and x86 and whatever other ISA, and
// carries with it all the downsides of inline asm:
// https://gcc.gnu.org/wiki/DontUseInlineAsm
// A second way would be for this function to explicitly disable
// -freciprocal-math (the flag associated with fast-math that causes the
// problem), but there appears to be no production worthy way to do this.  Gcc
// allows us to do
//#pragma GCC push_options
//#pragma GCC optimize ("-fno-reciprocal-math")
// ... insert function code here ...
//#pragma GCC pop_options
// but the gcc documentation states this is intended for debugging only.
// Alternatively gcc lets us use a function attribute:
// bool trial_divide(args) __attribute__ ((optimize(no-reciprocal-math)))
// but again the docs state it is intended for debugging only.  Additionally
// clang, icc, and msvc may all have different ways to disable the reciprocal.
// And a possible reason why this is intended as debug only, is that link-time
// optimization perhaps might reintroduce the optimization despite it being
// locally disabled for the function.
// A third way would work, which is to require clients of this library to
// compile their code without using -freciprocal-math (hence no fast-math which
// enables it), or explicitly disable it via -fno-reciprocal-math (and require
// similar measures for MSVC, etc).  Regardless of compiler, placing this
// requirement onto clients would be unacceptable because they might not know
// it is a requirement and might therefore compile this code in a broken form.
//
// In the end, what I would need is proof of correctness for this function
// both using standard double precision fp division, and with using the
// alternative of double precision fp reciprocal followed by double precision
// multiply (which might be generated due to -freciprocal-math).
// Until I have created a proof, or found a proof elsewhere, this function will
// remain experimental.
#include <immintrin.h>
bool trial_divide(std::uint32_t& div_result, std::uint32_t x, std::uint32_t n)
{
    // make sure type double is IEEE double precision fp (with 52 bit mantissa)
    static_assert(std::numeric_limits<double>::is_iec559, "");
    HPBC_PRECONDITION2(n > 0);   // disallow division by 0
    using std::uint32_t;

    double numerator = static_cast<double>(x);
    double denom = static_cast<double>(n);
#if 0
    // SSE intrinsics: very likely the compiler will utilize SSE (and any other
    // SIMD instruction set) automatically if it is available for the target and
    // if it can both meet the specifications for type double and improve
    // performance over a non-SIMD floating point unit.  Therefore I believe
    // there is unlikely to be any benefit in creating a dependence on a
    // particular SIMD ISA, and so I've disabled this SSE-specific section.
    double quotient;
    __m128d nv = _mm_set_sd(numerator);
    __m128d dv = _mm_set_sd(denom);
    __m128d qv = _mm_div_sd(nv, dv);
    _mm_store_sd(&quotient, qv);
#else
    double quotient = numerator / denom;
#endif
    div_result = static_cast<uint32_t>(quotient);

    // test whether the remainder (x - n*div_result) equals 0
    return (x == n*div_result);
}

#endif


}}  // end namespace

#endif
