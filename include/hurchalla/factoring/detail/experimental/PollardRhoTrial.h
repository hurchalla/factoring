// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_TRIAL_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_TRIAL_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


#ifndef HURCHALLA_POLLARD_RHO_GCD_THRESHOLD
#  define HURCHALLA_POLLARD_RHO_GCD_THRESHOLD 72
#endif


// Note: We generally expect Pollard Rho Brent factorization will be faster than
// Pollard Rho (this function).  However, this function has some advantages in
// instruction level parallelism that might make it faster than Pollard Rho
// Brent.  As usual, you'll need to benchmark to know.
// In any event, this function serves as a reference implementation of Pollard
// Rho factorization.

// The input parameter 'num' is the number we are factoring.  When using the
// Montgomery Form version of this function, 'num' is an implicit parameter that
// is set to the modulus of the montgomery form parameter 'mf'.  Num must be odd
// and must be > 2.  Num must not be prime (and thus the modulus of 'mf' when
// using Montgomery Form must not be prime); otherwise this function will always
// fail to find a factor (returning 0) after a potentially extremely long
// computation.
// Since this function can also fail for composite values of num, this function
// can not be used to determine primality.
//
// The input parameter 'c' is the constant addend in the generated sequence:
// a[n+1] = a[n]*a[n] + c
// typically a caller would first call this function with c == 1 and re-call the
// function with an incremented value of c any time the function failed to find
// a factor.
//
// A return of 0 indicates the function failed to find a factor.  This is not an
// error, but is instead a low probability, expected possibility inherent to the
// algorithm when factoring a composite number.  Whenever the function returns
// 0, we typically expect the caller will keep retrying this function using
// incremented values of 'c', until the function returns a non-zero value.

// If the function returns a non-zero value, then that return value is a factor
// of num.  Note that the return value is not guaranteed to be prime.

// This function was originally created based on the description at
// http://www.cs.colorado.edu/~srirams/classes/doku.php/pollard_rho_tutorial
// I added the trick of multiplying (modular mult) each iteration's absolute
// difference result with a loop carried product, as described in Variants at
// https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
// This allows us to calculate the gcd of num and product after some set number
// of loop iterations, rather than calling gcd once every iteration.
//
// basic algorithm
// ---------------
//  a = 2;
//  b = f(a);
//  while (a != b) {
//      p = GCD( | b - a | , num);
//      if ( p > 1)
//          return "Found factor: p";
//      a = f(a); // a runs once
//      b = f(f(b)); // b runs twice as fast.
//  }
//  return "Failed. :-("



/*
// The following version of the Pollard Rho function is a reference
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
T pollard_rho_trial(T num, T c)
{
    static_assert(ut_numeric_limits<U>::is_integer, "");
    static_assert(!(ut_numeric_limits<U>::is_signed), "");
    HPBC_PRECONDITION2(num > 2);
    HPBC_PRECONDITION2(c < num);

    // Require num odd, strictly for compatibility with the montgomery version
    // of this function below.  Regardless, we should be extremely surprised if
    // num is even; at the least it would signify an efficiency error, if not a
    // coding error - it would mean easy factors of 2 had not been factored out
    // prior to calling this more intensive factoring function.
    HPBC_PRECONDITION2(num % 2 == 1);

    T a = 2;
    T b = a;
    HPBC_ASSERT2(a < num);
    HPBC_ASSERT2(b < num);

    int gcd_threshold = HURCHALLA_POLLARD_RHO_GCD_THRESHOLD;

    T product = 1;
    HPBC_ASSERT2(product < num);

    namespace hc = ::hurchalla;
    while (true) {
        T absValDiff;
        for (int i = 0; i < gcd_threshold; ++i) {
            HPBC_INVARIANT2(1 <= product && product < num);

            // set a := (a*a + c) % num, while ensuring overflow doesn't happen.
            a = hc::modular_multiplication_prereduced_inputs(a, a, num);
            a = hc::modular_addition_prereduced_inputs(a, c, num);
            // set b := (b*b + c) % num, while ensuring overflow doesn't happen.
            b = hc::modular_multiplication_prereduced_inputs(b, b, num);
            b = hc::modular_addition_prereduced_inputs(b, c, num);
            // set b := (b*b + c) % num, while ensuring overflow doesn't happen.
            b = hc::modular_multiplication_prereduced_inputs(b, b, num);
            b = hc::modular_addition_prereduced_inputs(b, c, num);

            absValDiff = (a>b) ? a-b : b-a;
            T result = hc::modular_multiplication_prereduced_inputs(product,
                                                               absValDiff, num);
            if (result == 0) {
                // Since result == 0, we know that absValDiff == 0 -or-
                // product and absValDiff together had all the factors of
                // num (and likely more than just that), though neither could
                // have contained all factors since they're both always mod num.
                // It's fairly obvious product could have a factor when
                // absValDiff == 0 (they're uncorrelated) and the above shows
                // product might (in fact does given result == 0) have a
                // factor when absValDiff != 0.  So we need to test product
                // for a factor before checking for absValDiff == 0.
                HPBC_INVARIANT2(1 <= product && product < num);
                break;
            }
            product = result;
        }
        T p = hc::greatest_common_divisor(product, num);
        // Since product is in the range [1,num) and num is required to be > 1,
        // GCD will never return num or 0.  So we know the GCD will be in the
        // range [1, num).
        HPBC_ASSERT2(1<=p && p<num);
        if (p > 1)
            return p;
        if (absValDiff == 0)
            break;
    }
    return 0;  // the sequence cycled before we could find a factor
}
*/



// The following functor is an adaptation of the function above, using a
// Montgomery Form type M, and Montgomery domain values and arithmetic.
// For type M, ordinarily you'll use a template class instantiation of
// hurchalla/montgomery_arithmetic/MontgomeryForm.h

template <class M>
struct PollardRhoTrial {
  using T = typename M::IntegerType;
  using V = typename M::MontgomeryValue;
  using C = typename M::CanonicalValue;
  T operator()(const M& mf, T& expected_iterations, C c) const
  {
    T num = mf.getModulus();
    HPBC_PRECONDITION2(num > 2);
    HPBC_PRECONDITION2(!is_prime_miller_rabin::call(num));

    constexpr int gcd_threshold = HURCHALLA_POLLARD_RHO_GCD_THRESHOLD;

    // negate c so that we can use fusedSquareSub inside the loop instead of
    // fusedSquareAdd (fusedSquareSub may be slightly more efficient).
    C negative_c = mf.negate(c);

    V a = mf.getUnityValue();
    a = mf.add(a, a);   // sets a = mf.convertIn(2).

#if 0
// Using a pre-cycle doesn't seem to improve runtimes for this plain Pollard-Rho
// Trial (although it does help Pollard-Rho Brent Trials).
    constexpr int PRE_CYCLE_SIZE = 48;
    for (int i = 0; i < PRE_CYCLE_SIZE; ++i)
        a = mf.fusedSquareSub(a, negative_c);
#endif
    (void)expected_iterations;   // silence unused variable warning

    V b = a;
    V product = mf.getUnityValue();
    while (true) {
        V absValDiff;
        for (int i = 0; i < gcd_threshold; ++i) {
            HPBC_INVARIANT2(mf.convertOut(product) > 0);
            b = mf.fusedSquareSub(b, negative_c);
            b = mf.fusedSquareSub(b, negative_c);
            a = mf.fusedSquareSub(a, negative_c);

            // modular unordered subtract isn't the same as absolute value of a
            // subtraction, but it works equally well for pollard rho.
            absValDiff = mf.unorderedSubtract(a, b);
            bool isZero;
            // perform montgomery multiplication of product*absValDiff, and set
            // isZero to (mf.getCanonicalValue(result) == mf.getZeroValue()).
            V result = mf.multiply(product, absValDiff, isZero);
            if (isZero) {
                // Since result == 0, we know that absValDiff == 0 -or-
                // product and absValDiff together had all the factors of
                // num (and likely more than just that), though neither could
                // have contained all factors since they're both always mod num.
                // It's fairly obvious product could have a factor when
                // absValDiff == 0 (they're uncorrelated) and the above shows
                // that product must have a factor when absValDiff != 0.  So we
                // need to test product for a factor before checking whether
                // absValDiff == 0.
                break;
            }
            product = result;
        }

        // The following is a more efficient way to compute
        // p = greatest_common_divisor(mf.convertOut(product), num)
        T p = mf.gcd_with_modulus(product, [](auto x, auto y)
                       { return ::hurchalla::greatest_common_divisor(x, y); } );
        // Since product is in the range [1,num) and num is required to be >1,
        // GCD will never return num or 0.  So we know the GCD will be in the
        // range [1, num).
        HPBC_ASSERT2(1<=p && p<num);
        if (p > 1)
            return p;
        if (mf.getCanonicalValue(absValDiff) == mf.getZeroValue())
            break;
    }
    return 0;  // the sequence cycled before we could find a factor
  }
};


}}  // end namespace

#endif
