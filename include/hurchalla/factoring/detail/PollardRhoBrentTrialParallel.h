// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_BRENT_TRIAL_PARALLEL2_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_BRENT_TRIAL_PARALLEL2_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// See PollardRhoTrial.h for descriptions of the input parameters, template
// parameters, and the meaning of the return value.
// This functor's input and output is exactly the same.  Only the performance
// characteristics should be different from PollardRhoTrial.


#ifndef HURCHALLA_PRB_PARALELL2_GCD_THRESHOLD
#  define HURCHALLA_PRB_PARALELL2_GCD_THRESHOLD 608
#endif
#ifndef HURCHALLA_PRB_PARALLEL2_STARTING_LENGTH
#  define HURCHALLA_PRB_PARALLEL2_STARTING_LENGTH 19
#endif


#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4701)
#endif

// The following functor is an adaptation of the function above, using a
// Montgomery Form type M, and Montgomery domain values and arithmetic.
// For type M, ordinarily you'll use a template class instantiation of
// hurchalla/montgomery_arithmetic/MontgomeryForm.h
template <class M>
struct PollardRhoBrentTrialParallel {
  using T = typename M::IntegerType;
  using V = typename M::MontgomeryValue;
  using C = typename M::CanonicalValue;
  T operator()(const M& mf, T& expected_iterations, C c) const
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    T num = mf.getModulus();
    HPBC_PRECONDITION2(num > 2);
    HPBC_PRECONDITION2(!is_prime_miller_rabin::call(num));

    constexpr T gcd_threshold = HURCHALLA_PRB_PARALELL2_GCD_THRESHOLD;

    T advancement_len = HURCHALLA_PRB_PARALLEL2_STARTING_LENGTH;
    T best_advancement = (expected_iterations >> 4);
    if (advancement_len < best_advancement)
        advancement_len = best_advancement;
    T pre_length = 2*advancement_len + 2;

    V b1 = mf.getUnityValue();
    b1 = mf.add(b1, b1);                     // sets b1 = mf.convertIn(2)
    V b2 = mf.add(b1, mf.getUnityValue());   // sets b2 = mf.convertIn(3)

    // negate c so that we can use fusedSquareSub inside the loop instead of
    // fusedSquareAdd (fusedSquareSub may be slightly more efficient).
    C negative_c = mf.negate(c);

    for (T i = 0; i < pre_length; ++i) {
        b1 = mf.fusedSquareSub(b1, negative_c);
        b2 = mf.fusedSquareSub(b2, negative_c);
    }

    expected_iterations = pre_length;

    V product1 = mf.getUnityValue();
    V product2 = mf.getUnityValue();
    while (true) {
        V a_fixed1 = b1;
        V a_fixed2 = b2;
        for (T i = 0; i < advancement_len; ++i) {
            b1 = mf.fusedSquareSub(b1, negative_c);
            b2 = mf.fusedSquareSub(b2, negative_c);
        }
        expected_iterations += advancement_len;

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunsafe-loop-optimizations"
#endif
        // Because num is not prime, there are at least two factors, each of
        // which is <= sqrt(ut_numeric_limits<T>::max()).  Since the worst case
        // length for the smallest hidden cycle that we will eventually uncover
        // is a length no larger than the smallest factor, the variable
        // advancement_len will never grow larger than double that factor.  So
        // since advancement_len <= 2*smallestfactor, and
        // 2*smallestfactor <= 2*sqrt(ut_numeric_limits<T>::max()),
        // we know advancement_len <= 2*sqrt(ut_numeric_limits<T>::max()).
        //
        // Our only problem would be if gcd_threshold is (significantly) larger
        // than 1 << (ut_numeric_limits<T>::digits - 1), which perhaps could
        // happen if T is some tiny integral type.  In such a scenario, the
        // addition of (i + gcd_threshold) in the loop below maybe could
        // overflow.  Let's use static_assert to cause a compile error in that
        // situation, so we have no worry about overflow:
        static_assert(gcd_threshold <
                      (static_cast<T>(1) << (ut_numeric_limits<T>::digits -1)));

        for (T i=0; i < advancement_len; i = static_cast<T>(i+gcd_threshold)) {

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic pop
#endif
            T gcd_loop_len = (gcd_threshold < static_cast<T>(advancement_len-i))
                           ? gcd_threshold : static_cast<T>(advancement_len-i);

            V absValDiff1, absValDiff2;
            T expected_iterations_tmp = expected_iterations;
            for (T j = 0; j < gcd_loop_len; ++j) {
                b1 = mf.fusedSquareSub(b1, negative_c);
                b2 = mf.fusedSquareSub(b2, negative_c);

                HPBC_INVARIANT2(mf.convertOut(product1) > 0);
                HPBC_INVARIANT2(mf.convertOut(product2) > 0);
                // modular unordered subtract isn't the same as absolute value
                // of a subtraction, but it works equally well for pollard rho.
                absValDiff1 = mf.unorderedSubtract(a_fixed1, b1);
                absValDiff2 = mf.unorderedSubtract(a_fixed2, b2);
                bool isZero1, isZero2;
                // do montgomery multiplication of product*absValDiff, and set
                // isZero to (mf.getCanonicalValue(result) == mf.getZeroValue())
                V result1 = mf.multiply(product1, absValDiff1, isZero1);
                V result2 = mf.multiply(product2, absValDiff2, isZero2);
#if 1
                if (isZero1) {
#else
                if (isZero1 || isZero2) {
#endif
                    // Since result == 0, we know that absValDiff == 0 -or-
                    // product and absValDiff together had all the factors of
                    // num (and likely more than just that), though neither
                    // could have contained all factors since they're both
                    // always mod num.  It's fairly obvious product could have a
                    // factor when absValDiff == 0 (they're uncorrelated) and
                    // the above shows product must have a factor when
                    // absValDiff != 0.  So we need to test product for a factor
                    // before checking for absValDiff == 0.

#if 0
                    // I don't think this is necessary.  If result2 is legit but
                    // we don't do this, then product2 will just get one less
                    // multiply applied than it could have had.
                    // I wouldn't expect that to make any noticable difference.
                    if (!isZero2)
                        product2 = result2;
#endif
                    break;
                }
#if 1
                product1 = result1;
                if (isZero2) {
                    break;
                }
#else
                product1 = result1;
#endif
                product2 = result2;
                ++expected_iterations_tmp;
            }
            expected_iterations = expected_iterations_tmp;
            bool isZero;
            V product_combined = mf.multiply(product1, product2, isZero);
            if (isZero)   // product1 and product2 together had all the factors
                product_combined = product1;  // we'll do fine with 1 factor!

            // The following is a more efficient way to compute
            // p = greatest_common_divisor(mf.convertOut(product_combined), num)
            T p = mf.gcd_with_modulus(product_combined, [](auto x, auto y)
                       { return ::hurchalla::greatest_common_divisor(x, y); } );
            // Since product_combined is in the range [1,num) and num is
            // required to be > 1, GCD will never return num or 0.  So we know
            // the GCD will be in the range [1, num).
            HPBC_ASSERT2(1<=p && p<num);
            if (p > 1)
                return p;
            if (mf.getCanonicalValue(absValDiff1) == mf.getZeroValue())
                return 0; // sequence1 cycled before we could find a factor
            if (mf.getCanonicalValue(absValDiff2) == mf.getZeroValue())
                return 0; // sequence2 cycled before we could find a factor
        }
        advancement_len = static_cast<T>(2 * advancement_len);
    }
  }
};

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif


}}  // end namespace

#endif
