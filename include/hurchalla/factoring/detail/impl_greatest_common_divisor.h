// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IMPL_GREATEST_COMMON_DIVISOR_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_GREATEST_COMMON_DIVISOR_H_INCLUDED


#include "hurchalla/util/count_trailing_zeros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"

namespace hurchalla { namespace detail {


// Note: we use a struct with static functions in order to disallow ADL
struct impl_greatest_common_divisor {

// The binary GCD algorithm is usually considerably faster than the Euclidean
// GCD algorithm.  However for native types T, some new CPUs have very fast
// dividers that potentially could make the Euclidean GCD implementation faster
// than the Binary GCD implementation.  You can predefine the macro 
// HURCHALLA_PREFER_EUCLIDEAN_GCD in such a case.  You will also need to make
// sure that you have predefined HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.

#if defined(HURCHALLA_PREFER_EUCLIDEAN_GCD) && \
        defined(HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE)
  // Euclidean GCD.  For details on the Euclidean GCD, see
  // https://en.wikipedia.org/wiki/Greatest_common_divisor
  // Note that it's more efficient to provide inputs a<=b (otherwise the first
  // loop of the algorithm effectively performs swap(a,b) via a % operation).
  template <typename T>
  static HURCHALLA_FORCE_INLINE
  typename std::enable_if<(ut_numeric_limits<T>::digits <=
                         HURCHALLA_TARGET_BIT_WIDTH), T>::type
  call(T a, T b)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(a > 0 || b > 0);

    while (a != 0) {
        T tmp = a;
        a = static_cast<T>(b % a);
        b = tmp;
    }
    HPBC_POSTCONDITION2(b > 0);
    return b;
  }
#endif


  // Binary/Stein GCD.
  // For Binary GCD info see https://en.wikipedia.org/wiki/Binary_GCD_algorithm
  // This function is adapted from the Rust iterative version there.
  template <typename T>
#if defined(HURCHALLA_PREFER_EUCLIDEAN_GCD) && \
        defined(HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE)
  static HURCHALLA_FORCE_INLINE
  typename std::enable_if<!(ut_numeric_limits<T>::digits <=
                            HURCHALLA_TARGET_BIT_WIDTH), T>::type
#else
  static HURCHALLA_FORCE_INLINE T
#endif
  call(T u, T v)
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // We use this precondition solely to ensure a nonzero return value
    HPBC_PRECONDITION2(u > 0 || v > 0);

    namespace hc = ::hurchalla;
    if (u == 0) {
        HPBC_POSTCONDITION2(v > 0);
        return v;
    }
    if (v != 0) {
        int i = hc::count_trailing_zeros(u);
        int j = hc::count_trailing_zeros(v);
        u = static_cast<T>(u >> i);
        v = static_cast<T>(v >> j);
           //int k = (i < j) ? i : j;
        int k = hc::conditional_select((i < j), i, j);

        while (true) {
            HPBC_ASSERT2(u % 2 == 1);
            HPBC_ASSERT2(v % 2 == 1);
            T tmp = u;
            T sub1 = static_cast<T>(v - tmp);
            T sub2 = static_cast<T>(tmp - v);
            if (tmp == v)
                break;
               // u = (tmp >= v) ? v : tmp;
            u = hc::conditional_select((tmp >= v), v, tmp);
            // set v to the absolute value of (v - tmp)
               // v = (tmp >= v) ? sub2 : sub1;
            v = hc::conditional_select((tmp >= v), sub2, sub1);
            HPBC_ASSERT2(u % 2 == 1);
            // In an earlier version of this function, the line below used
            // count_trailing_zeros(v) instead of count_trailing_zeros(sub1),
            // which had been more of a standard way to write the algorithm.
            // But as pointed out by https://gmplib.org/manual/Binary-GCD, "in
            // twos complement the number of low zero bits on u-v is the same
            // as v-u, so counting or testing can begin on u-v without waiting
            // for abs(u-v) to be determined."  Hence we are able to use sub1
            // for the argument, which saves a few cycles because it lets the
            // CPU exploit instruction level parallelism.
            j = hc::count_trailing_zeros(sub1);
            v = static_cast<T>(v >> j);
        }
        u = static_cast<T>(u << k);
    }
    HPBC_POSTCONDITION2(u > 0);
    return u;
  }

}; // end struct impl_greatest_common_divisor


}} // end namespace

#endif
