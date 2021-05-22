// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_BRENT_TRIAL_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_BRENT_TRIAL_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// See PollardRhoTrial.h for descriptions of the input parameters, template
// parameters, and the meaning of the return value.
// This functor's input and output is exactly the same.  Only the performance
// characteristics should be different from PollardRhoTrial.


#ifndef HURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD
#  define HURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD 400
#endif
#ifndef HURCHALLA_POLLARD_RHO_BRENT_PRE_CYCLE_SIZE
#  define HURCHALLA_POLLARD_RHO_BRENT_PRE_CYCLE_SIZE 40
#endif
#ifndef HURCHALLA_POLLARD_RHO_BRENT_INITIAL_CYCLE_SIZE
#  define HURCHALLA_POLLARD_RHO_BRENT_INITIAL_CYCLE_SIZE 20
#endif


/*
// The following version of the Pollard Rho Brent function is a reference
// implementation using standard integers (not montgomery form) for type T.
// The Montgomery form version of this function that comes after was based upon
// this function.  Although this version of the function should work fine, it
// is intended as reference documentation, and the Montgomery Form version
// should be more efficient for nearly any CPU (though if you need to be
// absolutely certain you need to benchmark).
// Hence it currently exists only within the comments here.

// You'll need to move these #includes to the top of this file, and uncomment
// them.  They won't work if placed here.
//
//#include "hurchalla/factoring/greatest_common_divisor.h"
//#include "hurchalla/modular_arithmetic/modular_multiplication.h"
//#include "hurchalla/modular_arithmetic/modular_addition.h"

template <class T>
T pollard_rho_brent_trial(T num, T c)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    HPBC_PRECONDITION2(num > 2);
    HPBC_PRECONDITION2(c < num);
    // Require num odd, strictly for compatibility with the montgomery version
    // of this function below.  Regardless, we should be extremely surprised if
    // num is even; at the least it would signify an efficiency error, if not a
    // coding error - it would mean easy factors of 2 had not been factored out
    // prior to calling this more intensive factoring function.
    HPBC_PRECONDITION2(num % 2 == 1);

    constexpr T INITIAL_CYCLE_SIZE = 8;  // potentially any power of 2
    // If we used INITIAL_VALUE>2, we'd need to mod it by num (or we'd need to
    // require num > INITIAL_VALUE)
    constexpr T INITIAL_VALUE = 2;
    HPBC_ASSERT2(INITIAL_VALUE < num);

    T b = INITIAL_VALUE;
    T cycle_size = INITIAL_CYCLE_SIZE;
    T gcd_threshold = HURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD;

    T product = 1;
    HPBC_ASSERT2(product < num);

    while (true) {
        T a_fixed = b;
        for (T i = 0; i < cycle_size; ++i) {
            // set b := (b*b + c) % num, while ensuring overflow doesn't happen.
            b = modular_multiplication_prereduced_inputs(b, b, num);
            b = modular_addition_prereduced_inputs(b, c, num);
        }
        for (T i = 0; i < cycle_size; i += gcd_threshold) {
            T gcd_loop_size = (gcd_threshold < cycle_size - i) ? gcd_threshold
                                                               : cycle_size - i;
            T j, absValDiff;
            for (j = 0; j < gcd_loop_size; ++j) {
                HPBC_INVARIANT2(1 <= product && product < num);
                // set b := (b*b + c) % num, while ensuring no overflow
                b = modular_multiplication_prereduced_inputs(b, b, num);
                b = modular_addition_prereduced_inputs(b, c, num);
                absValDiff = (a_fixed > b) ? a_fixed - b : b - a_fixed;
                T result = modular_multiplication_prereduced_inputs(product,
                                                               absValDiff, num);
                if (result == 0) {
                    // Since result == 0, we know that absValDiff == 0 -or-
                    // product and absValDiff together had all the factors of
                    // num (and likely more than just that), though neither
                    // could have contained all factors since they're both
                    // always mod num.  It's fairly obvious product could have a
                    // factor when absValDiff == 0 (they're uncorrelated) and
                    // the above shows product must have a factor when
                    // absValDiff != 0.  So we need to test product for a factor
                    // before checking for absValDiff == 0.
                    HPBC_INVARIANT2(1 <= product && product < num);
                    break;
                }
                product = result;
            }
            T p = greatest_common_divisor(product, num);
            // Since product is in the range [1,num) and num is required to
            // be > 1, GCD will never return num or 0.  So we know the GCD will
            // be in the range [1, num).
            HPBC_ASSERT2(1<=p && p<num);
            if (p > 1)
                return p;
            if (absValDiff == 0)
                return 0; // the sequence cycled before we could find a factor
        }
        cycle_size *= 2;
    }
}
*/


// The following functor is an adaptation of the function above, using a
// Montgomery Form type M, and Montgomery domain values and arithmetic.
// For type M, ordinarily you'll use a template class instantiation of
// hurchalla/montgomery_arithmetic/MontgomeryForm.h
template <class M>
struct PollardRhoBrentTrial {
typename M::T_type operator()(const M& mf, typename M::CanonicalValue c) const
{
    using T = typename M::T_type;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    T num = mf.getModulus();
    HPBC_PRECONDITION2(num > 2);
    HPBC_PRECONDITION2(!is_prime_miller_rabin_integral(num));

    constexpr T gcd_threshold = HURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD;
    constexpr T pre_cycle_size = HURCHALLA_POLLARD_RHO_BRENT_PRE_CYCLE_SIZE;

    T cycle_size = HURCHALLA_POLLARD_RHO_BRENT_INITIAL_CYCLE_SIZE;

    V b = mf.getUnityValue();
    b = mf.add(b, b);   // sets b = mf.convertIn(2)
    // negate c so that we can use fmsub inside the loop instead of fmadd (fmsub
    // is very slightly more efficient).
    C negative_c = mf.negate(c);

    for (T i = 0; i < pre_cycle_size; ++i)
        b = mf.fmsub(b, b, negative_c);

    V product = mf.getUnityValue();
    while (true) {
        V a_fixed = b;
        for (T i = 0; i < cycle_size; ++i) {
            b = mf.fmsub(b, b, negative_c);
        }

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunsafe-loop-optimizations"
#endif
        // Because num is not prime, there are at least two factors, each of
        // which is <= sqrt(ut_numeric_limits<T>::max()).  Since the worst case
        // length for the smallest hidden cycle that we will eventually uncover
        // is a length no larger than the smallest factor, the variable
        // cycle_size will never grow larger than double that factor.  So since
        // cycle_size <= 2*smallestfactor <= 2*sqrt(ut_numeric_limits<T>::max())
        // we know cycle_size <= 2*sqrt(ut_numeric_limits<T>::max())
        //
        // Our only problem would be if gcd_threshold is (significantly) larger
        // than 1 << (ut_numeric_limits<T>::digits - 1), which perhaps could
        // happen if T is some tiny integral type.  In such a scenario, the
        // addition of (i + gcd_threshold) in the loop below maybe could
        // overflow.  Let's use static_assert to cause a compile error in that
        // situation, so we have no worry about overflow:
        static_assert(gcd_threshold <
                      (static_cast<T>(1) << (ut_numeric_limits<T>::digits -1)));

        for (T i = 0; i < cycle_size; i = static_cast<T>(i + gcd_threshold)) {

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic pop
#endif
            T gcd_loop_size = (gcd_threshold < static_cast<T>(cycle_size - i))
                               ? gcd_threshold : static_cast<T>(cycle_size - i);
            V absValDiff;
            T j;
            for (j = 0; j < gcd_loop_size; ++j) {
                HPBC_INVARIANT2(mf.convertOut(product) > 0);
                // modular unordered subtract isn't the same as absolute value
                // of a subtraction, but it works equally well for pollard rho.
                absValDiff = mf.unorderedSubtract(a_fixed, b);
                bool isZero;
                // do montgomery multiplication of product*absValDiff, and set
                // isZero to (mf.getCanonicalValue(result) == mf.getZeroValue())
                V result = mf.multiplyIsZero(product, absValDiff, isZero);
                if (isZero) {
                    // Since result == 0, we know that absValDiff == 0 -or-
                    // product and absValDiff together had all the factors of
                    // num (and likely more than just that), though neither
                    // could have contained all factors since they're both
                    // always mod num.  It's fairly obvious product could have a
                    // factor when absValDiff == 0 (they're uncorrelated) and
                    // the above shows product must have a factor when
                    // absValDiff != 0.  So we need to test product for a factor
                    // before checking for absValDiff == 0.
                    break;
                }
                product = result;
                b = mf.fmsub(b, b, negative_c);
            }
            // The following is a more efficient way to compute
            // p = greatest_common_divisor(mf.convertOut(product), num)
            T p = mf.gcd_with_modulus(product, [](auto x, auto y)
                                    { return greatest_common_divisor(x, y); } );
            // Since product is in the range [1,num) and num is required to
            // be > 1, GCD will never return num or 0.  So we know the GCD will
            // be in the range [1, num).
            HPBC_ASSERT2(1<=p && p<num);
            if (p > 1) 
                return p;
            if (mf.getCanonicalValue(absValDiff) == mf.getZeroValue())
                return 0; // the sequence cycled before we could find a factor
        }
        cycle_size = static_cast<T>(2 * cycle_size);
    }
}
};


}}  // end namespace

#endif
