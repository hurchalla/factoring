// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_BRENT_TRIAL_PARALLEL3_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_BRENT_TRIAL_PARALLEL3_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// See PollardRhoTrial.h for descriptions of the input parameters, template
// parameters, and the meaning of the return value.
// This functor's input and output is exactly the same.  Only the performance
// characteristics should be different from PollardRhoTrial.


#ifndef HURCHALLA_PRB_GCD_THRESHOLD
#  define HURCHALLA_PRB_PARALELL3_GCD_THRESHOLD 608
#else
#  define HURCHALLA_PRB_PARALELL3_GCD_THRESHOLD HURCHALLA_PRB_GCD_THRESHOLD
#endif
#ifndef HURCHALLA_PRB_STARTING_LENGTH
#  define HURCHALLA_PRB_PARALLEL3_STARTING_LENGTH 19
#else
#  define HURCHALLA_PRB_PARALLEL3_STARTING_LENGTH HURCHALLA_PRB_STARTING_LENGTH
#endif


#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4701)
#endif

// The following functor is a simple parallel (instruction level, not thread)
// variation on Pollard-Rho-Brent that advances two separate sequences.
// The instruction level parallelism allows most CPUs to advance two sequences
// with only slightly more cost than one sequence.

// M must be a Montgomery Form type; ordinarily you'd use a template class
// instantiation of hurchalla/montgomery_arithmetic/MontgomeryForm.h
template <class M>
struct PollardRhoBrentTrialParallelOpt {
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

    constexpr T gcd_threshold = HURCHALLA_PRB_PARALELL3_GCD_THRESHOLD;

    T advancement_len = HURCHALLA_PRB_PARALLEL3_STARTING_LENGTH;
    T best_advancement = static_cast<T>(expected_iterations >> 4);
    if (advancement_len < best_advancement)
        advancement_len = best_advancement;
    T pre_length = static_cast<T>(2*advancement_len + 2);

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

    V a_fixed1 = b1;
    for (T i = 0; i < advancement_len; ++i) {
        b1 = mf.fusedSquareSub(b1, negative_c);
        b2 = mf.fusedSquareSub(b2, negative_c);
    }

    expected_iterations = static_cast<T>(pre_length + advancement_len);

    V product = mf.getUnityValue();
    while (true) {

        V a_fixed2 = b2;

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
            V absValDiff;
            T expected_iterations_tmp = expected_iterations;
            { T j = 0; do {
                b1 = mf.fusedSquareSub(b1, negative_c);
                b2 = mf.fusedSquareSub(b2, negative_c);

                HPBC_INVARIANT2(mf.convertOut(product) > 0);
                // modular unordered subtract isn't the same as absolute value
                // of a subtraction, but it works equally well for pollard rho.
                absValDiff = mf.unorderedSubtract(a_fixed1, b1);
                bool isZero;
                // do montgomery multiplication of product*absValDiff, and set
                // isZero to (mf.getCanonicalValue(result) == mf.getZeroValue())
                V result = mf.multiply(product, absValDiff, isZero);
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
                ++expected_iterations_tmp;
            } while (++j < gcd_loop_len); }

            expected_iterations = expected_iterations_tmp;

            // The following is a more efficient way to compute 
            // p = greatest_common_divisor(mf.convertOut(product), num)
            T p = mf.gcd_with_modulus(product, [](auto x, auto y)
                       { return ::hurchalla::greatest_common_divisor(x, y); } );
            // Since product is in the range [1,num) and num is required to
            // be > 1, GCD will never return num or 0.  So we know the GCD will
            // be in the range [1, num).
            HPBC_ASSERT2(1<=p && p<num);
            if (p > 1)
                return p;
            if (mf.getCanonicalValue(absValDiff) == mf.getZeroValue())
                return 0; // the sequence cycled before we could find a factor
        }

        a_fixed1 = b1;

        T quarter_len = static_cast<T>(advancement_len / 2);
        advancement_len = static_cast<T>(2 * advancement_len);

        { T i = 0; do {
            b1 = mf.fusedSquareSub(b1, negative_c);
            b2 = mf.fusedSquareSub(b2, negative_c);
        } while (++i < quarter_len); }
        expected_iterations = static_cast<T>(expected_iterations + quarter_len);

        T three_quarter_len = static_cast<T>(advancement_len - quarter_len);

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunsafe-loop-optimizations"
#endif
        // for info on this loop, see comments in the gcd loop above.
        for (T i=0; i < three_quarter_len; i= static_cast<T>(i+gcd_threshold)) {

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic pop
#endif
            T gcd_loop_len=(gcd_threshold < static_cast<T>(three_quarter_len-i))
                          ? gcd_threshold : static_cast<T>(three_quarter_len-i);
            V absValDiff;
            T expected_iterations_tmp = expected_iterations;
            { T j = 0; do {
                b1 = mf.fusedSquareSub(b1, negative_c);
                b2 = mf.fusedSquareSub(b2, negative_c);

                HPBC_INVARIANT2(mf.convertOut(product) > 0);
                absValDiff = mf.unorderedSubtract(a_fixed2, b2);
                bool isZero;
                V result = mf.multiply(product, absValDiff, isZero);
                if (isZero)
                    break;
                product = result;
                ++expected_iterations_tmp;
            } while (++j < gcd_loop_len); }

            expected_iterations = expected_iterations_tmp;

            T p = mf.gcd_with_modulus(product, [](auto x, auto y)
                       { return ::hurchalla::greatest_common_divisor(x, y); } );
            HPBC_ASSERT2(1<=p && p<num);
            if (p > 1)
                return p;
            if (mf.getCanonicalValue(absValDiff) == mf.getZeroValue())
                return 0;
        }
/*
    T threshold;
    if (advancement_len <= gcd_threshold + gcd_threshold/4)
        threshold = advancement_len;
    else if (advancement_len <= 2*gcd_threshold + gcd_threshold/4 + gcd_threshold/8)   // 1.25 to 2.375
        threshold = 1 + advancement_len / 2;  // .6125 to 1.1875  * gcd_threshold
    else if (advancement_len <= 3*gcd_threshold + gcd_threshold/2)   // 2.375 to 3.5
        threshold = 1 + advancement_len / 3;  // .792 to 1.167  * gcd_threshold
    else if (advancement_len <= 4*gcd_threshold + gcd_threshold/2)   // 3.5 to 4.5
        threshold = 1 + advancement_len / 4;  // 0.875 to 1.125  * gcd_threshold
    else
        threshold = gcd_threshold;   // or  1 + advancement_len / (advancement_len / gcd_threshold);
*/
    }
  }
};

#if defined(_MSC_VER)
#  pragma warning(pop)
#endif


}}  // end namespace

#endif
