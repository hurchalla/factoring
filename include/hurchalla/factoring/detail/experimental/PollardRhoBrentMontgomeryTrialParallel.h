// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_BRENT_MONTGOMERY_TRIAL_PARALLEL_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_BRENT_MONTGOMERY_TRIAL_PARALLEL_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/factoring/detail/PollardRhoTrial.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"

namespace hurchalla { namespace detail {


// See PollardRhoTrial.h for descriptions of the input parameters, template
// parameters, and the meaning of the return value.
// This functor's input and output is exactly the same.  Only the performance
// characteristics should be different from PollardRhoTrial.


#ifndef HURCHALLA_PRBMTP_ONE_THIRD_GCD_THRESHOLD
#  define HURCHALLA_PRBMTP_ONE_THIRD_GCD_THRESHOLD 256
#endif


// This file is based on section 3(Brent's Improvement to Monte Carlo) of
// "Speeding the Pollard and Elliptic Curve Methods of Factorization", by
// Peter L. Montgomery.
// https://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866113-7/S0025-5718-1987-0866113-7.pdf


// The following functor is an adaptation of the function above, using a
// Montgomery Form type M, and Montgomery domain values and arithmetic.
// For type M, ordinarily you'll use a template class instantiation of
// hurchalla/montgomery_arithmetic/MontgomeryForm.h
template <class M>
struct PollardRhoBrentMontgomeryTrialParallel {
typename M::T_type operator()(const M& mf, typename M::CanonicalValue c)
{
    using T = typename M::T_type;
    using V = typename M::MontgomeryValue;
    using C = typename M::CanonicalValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");

    T num = mf.getModulus();
    HPBC_PRECONDITION2(num > 2);
    HPBC_PRECONDITION2(!is_prime_miller_rabin_integral(num));

    // If we used INITIAL_VALUE>2, we'd need to mod it by num (or we'd need to
    // require num > INITIAL_VALUE)
    //constexpr T INITIAL_VALUE = 2;    // note that b gets set to 2 below
    //HPBC_ASSERT2(INITIAL_VALUE < num);
    constexpr T PRE_CYCLE_SIZE = 48;

    constexpr T ONE_THIRD_INITIAL_CYCLE_SIZE = 8;
    T one_third_gcd_threshold = HURCHALLA_PRBMTP_ONE_THIRD_GCD_THRESHOLD;

    T one_third_cycle_size = ONE_THIRD_INITIAL_CYCLE_SIZE;
    V b = mf.getUnityValue();
    b = mf.add(b, b);   // sets b = mf.convertIn(2)
    V bz = mf.add(b, b);         // sets bz = mf.convertIn(4)
    // negate c so that we can use fmsub inside the loop instead of fmadd (fmsub
    // is very slightly more efficient).
    C negative_c = mf.negate(c);

    for (T i = 0; i < PRE_CYCLE_SIZE; ++i) {
        b = mf.fmsub(b, b, negative_c);
        bz = mf.fmsub(bz, bz, negative_c);
    }

    V product = mf.getUnityValue();
    V productz = mf.getUnityValue();
    while (true) {
        V t1 = b;
        // for negt2, we want -(b*b+c)
        b = mf.multiply(b, b);
        b = mf.subtract(negative_c, b);
        V negt2 = b;
        b = mf.fmsub(b, b, negative_c);

        V t1z = bz;
        bz = mf.multiply(bz, bz);
        bz = mf.subtract(negative_c, bz);
        V negt2z = bz;
        bz = mf.fmsub(bz, bz, negative_c);

        V t3 = b;
        V a1 = mf.subtract(t3, negt2);
        V a2 = mf.add(a1, t1);
        V tmpna3 = mf.fmadd(a1, t1, negative_c);
        V nega3 = mf.subtract(mf.multiply(negt2, t3), tmpna3);
        C nega4 = mf.getCanonicalValue(mf.multiply(a1, mf.add(negt2, nega3)));

        V t3z = bz;
        V a1z = mf.subtract(t3z, negt2z);
        V a2z = mf.add(a1z, t1z);
        V tmpna3z = mf.fmadd(a1z, t1z, negative_c);
        V nega3z = mf.subtract(mf.multiply(negt2z, t3z), tmpna3z);
        C nega4z = mf.getCanonicalValue(mf.multiply(a1z, mf.add(negt2z, nega3z)));

        for (T i = 0; i < 3*one_third_cycle_size; ++i) {
            b = mf.fmsub(b, b, negative_c);
            bz = mf.fmsub(bz, bz, negative_c);
        }
#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunsafe-loop-optimizations"
#endif
        // Because num is not prime, there are at least two factors, each of
        // which is <= sqrt(ut_numeric_limits<T>::max()).  Since the smallest
        // hidden cycle that we will eventually uncover has a length equal to
        // the smallest factor, cycle_size will never grow larger than double
        // that factor.  And since
        // cycle_size <= 2*smallestfactor <= 2*sqrt(ut_numeric_limits<T>::max())
        // we know cycle_size <= 2*sqrt(ut_numeric_limits<T>::max())
        //
        // Our only problem would be if one_third_gcd_threshold is
        // (significantly) larger than 1 << (ut_numeric_limits<T>::digits - 1),
        // which perhaps could happen if T is some tiny integral type.  In such
        // a scenario, the addition of (i + one_third_gcd_threshold) in the loop
        // below maybe could overflow.  Let's use static_assert to cause a
        // compile error in that situation, so we have no worry about overflow:
        static_assert(one_third_gcd_threshold <
                      (static_cast<T>(1) << (ut_numeric_limits<T>::digits -1)));

        for (T i = 0; i < one_third_cycle_size; i = static_cast<T>(i + one_third_gcd_threshold)) {

#if defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER)
#  pragma GCC diagnostic pop
#endif
            T gcd_loop_size = (one_third_gcd_threshold < static_cast<T>(one_third_cycle_size - i))
                               ? one_third_gcd_threshold : static_cast<T>(one_third_cycle_size - i);
            V gxi, gxiz;
            T j;
            for (j = 0; j < gcd_loop_size; ++j) {
                HPBC_INVARIANT2(mf.convertOut(product) > 0);
                HPBC_INVARIANT2(mf.convertOut(productz) > 0);

                V b2 = mf.fmsub(b, b, negative_c);
                V diffa2 = mf.subtract(b, a2);
                V diffna3 = mf.subtract(b2, nega3);
                gxi = mf.fmsub(diffa2, diffna3, nega4);

                V b2z = mf.fmsub(bz, bz, negative_c);
                V diffa2z = mf.subtract(bz, a2z);
                V diffna3z = mf.subtract(b2z, nega3z);
                gxiz = mf.fmsub(diffa2z, diffna3z, nega4z);

                bool isZero, isZeroz;
                V result = mf.multiplyIsZero(product, gxi, isZero);
                V resultz = mf.multiplyIsZero(productz, gxiz, isZeroz);
                if (isZero) {
                    // Since result == 0, we know that gxi == 0 -or-
                    // product and gxi together had all the factors of
                    // num (and likely more than just that), though neither
                    // could have contained all factors since they're both
                    // always mod num.  It's fairly obvious product could have a
                    // factor when gxi == 0 (they're uncorrelated) and
                    // the above shows product must have a factor when
                    // gxi != 0.  So we need to test product for a factor
                    // before checking for gxi == 0.
                    if (!isZeroz)
                        productz = resultz;
#if 1
                    // This might be a bit unnecessary... I'm not sure how important
                    // it is, but it should put us on similar ground to plain
                    // PollardRhoBrent.
                    if (mf.getCanonicalValue(gxi) == mf.getZeroValue()) {
                        // if gxi == 0, backtrack slightly
                        V absValDiff1 = mf.unorderedSubtract(b, t1);
                        V result2 = mf.multiplyIsZero(product, absValDiff1, isZero);
                        if (isZero)
                            break;
                        product = result2;
                        V absValDiff3 = mf.unorderedSubtract(b, t3);
                        V result3 = mf.multiplyIsZero(product, absValDiff3, isZero);
                        if (isZero)
                            break;
                        product = result3;
                        // we skip absvaldiff2 (we could have skippedabsvaldiff3
                        // instead - the order doesn't matter), because when gxi ==0
                        // this means that one of the absvaldiffs is zero, or
                        // the product of all three absvaldiffs is zero.  We don't
                        // bother to do anything special if all three absvaldiffs
                        // are zero and treat that simply as the sequences cycled
                        // before we could find a factor.  However when all factors
                        // are in the product of the three absvaldiffs, since there
                        // must be at least two non-trivial factors, this means we
                        // only need to look at two of the absvaldiffs to find at
                        // least one factor (which will get incorporated into
                        // product above when the absvaldiff != 0)
                    }
#endif
                    break;
                }
                product = result;

                if (isZeroz) {
#if 1
                    if (mf.getCanonicalValue(gxiz) == mf.getZeroValue()) {
                        V absValDiff1z = mf.unorderedSubtract(bz, t1z);
                        V result2z = mf.multiplyIsZero(productz, absValDiff1z, isZeroz);
                        if (isZeroz)
                            break;
                        productz = result2z;
                        V absValDiff3z = mf.unorderedSubtract(bz, t3z);
                        V result3z = mf.multiplyIsZero(productz, absValDiff3z, isZeroz);
                        if (isZeroz)
                            break;
                        productz = result3z;
                    }
#endif
                    break;
                }
                productz = resultz;

                b = mf.fmsub(b2, b2, negative_c);
                b = mf.fmsub(b, b, negative_c);

                bz = mf.fmsub(b2z, b2z, negative_c);
                bz = mf.fmsub(bz, bz, negative_c);
            }

            bool isZero3;
            V result3 = mf.multiplyIsZero(product, productz, isZero3);
            if (isZero3) {
                // we know neither product nor productz are 0, due to loop
                // invariant above.  So together they must have all factors,
                // though neither one has all of them by itself (that product
                // would have been zero if that were the case).  So both must
                // have at least one factor, and we can just pick one to use
                // for the gcd.
                result3 = product;   // or we could have used productz
            }

            // The following is a more efficient way to compute
            // p = greatest_common_divisor(mf.convertOut(product), num)
            T p = mf.template gcd_with_modulus<PrtGcdFunctor>(result3);
            // Since product is in the range [1,num) and num is required to
            // be > 1, GCD will never return num or 0.  So we know the GCD will
            // be in the range [1, num).
            HPBC_ASSERT2(1<=p && p<num);
            if (p > 1) 
                return p;
            if (mf.getCanonicalValue(gxi) == mf.getZeroValue()) 
                return 0; // the sequence cycled before we could find a factor
            if (mf.getCanonicalValue(gxiz) == mf.getZeroValue()) 
                return 0; // the sequence cycled before we could find a factor
        }
        one_third_cycle_size = static_cast<T>(2*one_third_cycle_size);
    }
}
};


}}  // end namespace

#endif
