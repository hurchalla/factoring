// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_TRIAL_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_TRIAL_H_INCLUDED


#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace factoring {


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
//  while (a != b)
//  {
//      p = GCD( | b - a | , num);
//      if ( p > 1)
//          return "Found factor: p";
//      a = f(a); // a runs once
//      b = f(f(b)); // b runs twice as fast.
//  }
//  return "Failed. :-("



// TODO determine a good number for HURCHALLA_POLLARD_RHO_TRIAL_GCD_THRESHOLD
#ifndef HURCHALLA_POLLARD_RHO_TRIAL_GCD_THRESHOLD
#  define HURCHALLA_POLLARD_RHO_TRIAL_GCD_THRESHOLD 70
#endif


/*
// The following version of the Pollard Rho function is a reference
// implementation using standard integers (not montgomery form) for type T.
// The Montgomery form version of this function that comes after was based upon
// this function.  Although this version of the function should work fine, it
// is intended as reference documentation, and the Montgomery Form version
// should be more efficient for nearly any CPU (though if you need to be
// absolutely certain you need to benchmark).
// Hence it currently exists only within the comments here.

#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"

template <class T>
T pollard_rho_trial(T num, T c, T* pIterationsPerformed = nullptr)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<U>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<U>::is_signed), "");
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

    int gcd_threshold = HURCHALLA_POLLARD_RHO_TRIAL_GCD_THRESHOLD;

    T product = 1;
    HPBC_ASSERT2(product < num);

    T current_iteration = 0;
    while (true)
    {
        T absValDiff;
        for (int i = 0; i < gcd_threshold; ++i)
        {
            HPBC_INVARIANT2(1 <= product && product < num);
            ++current_iteration;

            // set a := (a*a + c) % num, while ensuring overflow doesn't happen.
            a = ma::modular_multiplication_prereduced_inputs(a, a, num);
            a = ma::modular_addition_prereduced_inputs(a, c, num);
            // set b := (b*b + c) % num, while ensuring overflow doesn't happen.
            b = ma::modular_multiplication_prereduced_inputs(b, b, num);
            b = ma::modular_addition_prereduced_inputs(b, c, num);
            // set b := (b*b + c) % num, while ensuring overflow doesn't happen.
            b = ma::modular_multiplication_prereduced_inputs(b, b, num);
            b = ma::modular_addition_prereduced_inputs(b, c, num);

            absValDiff = (a>b) ? a-b : b-a;
            T result = ma::modular_multiplication_prereduced_inputs(product,
                                                               absValDiff, num);
            if (result == 0)
            {
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
        T p = greatest_common_divisor(product, num);
        // Since product is in the range [1,num) and num is required to be > 1,
        // GCD will never return num or 0.  So we know the GCD will be in the
        // range [1, num).
        HPBC_ASSERT2(1<=p && p<num);

        if (pIterationsPerformed != nullptr)
            *pIterationsPerformed = current_iteration;
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
// For type M, ordinarily you'll use a template class instantiation from
// hurchalla/montgomery_arithmetic/MontgomeryForm.h

template <class M>
struct PollardRhoTrial {
typename M::T_type operator()(const M& mf, typename M::CanonicalValue c,
                             typename M::T_type* pIterationsPerformed = nullptr)
{
    using T = typename M::T_type;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;

    T num = mf.getModulus();
    HPBC_ASSERT2(num > 2);
    HPBC_ASSERT2(num % 2 == 1);  // montgomery form guarantees odd modulus.

    V a = mf.getUnityValue();
    a = mf.add(a, a);   // sets a = mf.convertIn(2).
    V b = a;
    // negate c so that we can use fmsub inside the loop instead of fmadd (fmsub
    // is very slightly more efficient).
    C negative_c = mf.negate(c);

    int gcd_threshold = HURCHALLA_POLLARD_RHO_TRIAL_GCD_THRESHOLD;

    V product = mf.getUnityValue();
    T current_iteration = 0;
    while (true)
    {
        V absValDiff;
        for (int i = 0; i < gcd_threshold; ++i)
        {
            ++current_iteration;
            a = mf.fmsub(a, a, negative_c);
            b = mf.fmsub(b, b, negative_c);
            b = mf.fmsub(b, b, negative_c);

            // modular unordered subtract isn't the same as absolute value of a
            // subtraction, but it works equally well for pollard rho.
            absValDiff = mf.unorderedSubtract(a, b);
            bool isZero;
            // perform montgomery multiplication of product*absValDiff, and set
            // isZero to (mf.getCanonicalValue(result) == mf.getZeroValue()).
            V result = mf.multiplyIsZero(product, absValDiff, isZero);
            if (isZero)
            {
                // Since result == 0, we know that absValDiff == 0 -or-
                // product and absValDiff together had all the factors of
                // num (and likely more than just that), though neither could
                // have contained all factors since they're both always mod num.
                // It's fairly obvious product could have a factor when
                // absValDiff == 0 (they're uncorrelated) and the above shows
                // product might (in fact does, given result == 0) have a
                // factor when absValDiff != 0.  So we need to test product
                // for a factor before checking for absValDiff == 0.
                break;
            }
            product = result;
        }

        // For discussion purposes, let R = 2^(ma_numeric_limits<T>::digits).
        // for example if T is uint64_t, then R = 2^64.
        // Let the integer b = mf.convertOut(product).  We can note that b is
        // the standard integer domain representation of product, or
        // equivalently that product is the montgomery domain representation of
        // b.
        // Since product is a MontgomeryValue, if we use mathematical integers
        // (infinite precision and no overflow), we know product satisfies:
        // product ≡ b*R (mod num);  therefore there exists some integer k such
        // that  product == b*R + k*num.  Since num is the montgomery form's
        // modulus, num is odd, and therefore num and R are coprime.

        // Proof that gcd(product, num) == gcd(b, num)
        // -------------------------------------------
        // If some integer d is a divisor of num, then obviously d divides
        // k*num.  Since num and R are coprime, d can not be a divisor of R
        // (unless d==1);  therefore d divides b*R if and only if d divides b.
        // Thus if d divides b, d divides b*R, and since d always divides k*num,
        // d divides b*R + k*num.  Thus d divides b implies d divides product.
        // Likewise if d divides product, then d divides b*R.  And since num and
        // R are coprime and d is a divisor of num, we know d does not divide R
        // (unless d == 1), and thus d must be a divisor of b.
        // Therefore d divides product if and only if d divides b.
        // Let p be the greatest common divisor of product and num.  Since p
        // divides product, p must divide b.  Let  q = gcd(b, num).  Since q
        // divides b, q must divide product.  Since q also divides num, q is a
        // common divisor of product and num, and thus q <= p.  Since p divides
        // both b and num, p is a common divisor of b and num.  Since q is the
        // gcd(b,num) and q<=p,  we therefore know q == p.  Thus:
        // gcd(product, num) == gcd(b, num).

        // We want to set  p = gcd(mf.convertOut(product), num).  By the proof
        // above, we can use the more efficient  p = gcd(product, num).
        T p = greatest_common_divisor(product.get(), num);
        // Since product is in the range [1,num) and num is required to be >1,
        // GCD will never return num or 0.  So we know the GCD will be in the
        // range [1, num).
        HPBC_ASSERT2(1<=p && p<num);

        if (pIterationsPerformed != nullptr)
            *pIterationsPerformed = current_iteration;
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
