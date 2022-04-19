// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FACTORIZE_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_H_INCLUDED


#include "hurchalla/factoring/detail/impl_factorize.h"
#include "hurchalla/factoring/detail/PollardRhoIsPrime.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstddef>
#include <array>
#include <vector>
#include <numeric>
#include <functional>

namespace hurchalla {


// factorize() uses the Pollar-Rho factorization algorithm, with Brent's
// improvements.  See https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
// I have made novel improvements to the algorithm; both the code and
// algorithms are well-optimized.
// Prior to Pollard-Rho there is a small prime trial division stage.  Upon
// beginning Pollard-Rho, we test for primality before trying to extract each
// factor, by using the the deterministic Miller-Rabin algorithm - we usually
// speed up this algorithm by using one of the very small hash tables (~100
// bytes for example) from
// factoring/include/hurchalla/factoring/detail/miller_rabin_bases/
//
// This resulting factorization algorithm/function is likely to be the fastest
// method available for factoring arbitrary 64 bit numbers.
// [Again for arbitary inputs, Hart's One Line Factoring algorithm and/or
// Lehman's method have a good chance to be fastest for factoring non-large
// 32 bit numbers; ECM will likely be fastest for 128 bit and 256 bit numbers;
// and for larger numbers still, see Quadratic Sieve and GNFS.]

// ------------------------------------

// Returns a std::array that contains all factors of x, and writes the total
// number of factors to num_factors.  The array entries with index < num_factor
// are the factors.
// Note that ut_numeric_limits is a utility class that is used here to
// automatically get the bit width of your type T.  For example, if your type T
// is uint32_t, then this function returns std::array<uint32_t, 32>.  (The bit
// width has significance because it is impossible to have more factors than
// type T's bit width.)
// T can be any integral type <= 128 bits.
template <typename T>
std::array<T, ut_numeric_limits<T>::digits>
factorize(T x, int& num_factors)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::digits <= 128, "");
    HPBC_PRECONDITION(x >= 2);  // 0 and 1 do not have prime factorizations

    namespace hd = ::hurchalla::detail;
    // The max possible number of factors occurs when all factors equal 2,
    // so ut_numeric_limits<T>::digits is sufficient to hold all factors.
    std::array<T, ut_numeric_limits<T>::digits> arr =
                hd::impl_factorize::factorize_to_array(
                                       x, num_factors, hd::PollardRhoIsPrime());
    // After calling this function, a client should not index the returned
    // array with an index >= num_factors.  As a defensive measure, we'll set
    // all array entries at or beyond num_factors to 0 - this may help to make
    // an indexing error more obvious if a caller later makes this mistake.
    for (auto i = static_cast<std::size_t>(num_factors); i < arr.size(); ++i)
        arr[i] = 0;

    HPBC_POSTCONDITION(num_factors > 0);
    HPBC_POSTCONDITION(static_cast<std::size_t>(num_factors) <= arr.size());
    // all the factors multiplied together should == x
    HPBC_POSTCONDITION(x == std::accumulate(arr.begin(),
             arr.begin()+num_factors, static_cast<T>(1), std::multiplies<T>()));
    return arr;
}


// This version of factorize takes a std::vector, which it clears of any
// existing elements and then fills with the factors of x.  When this function
// returns, the size of the vector is the number of factors.
// Note that this version might be preferable to the array version of factorize
// if you wish to save stack space; vector is heap allocated but the array
// version needs memory on the stack for its returned array (for example with
// type T of uint64_t, the array needs 512 bytes on stack).
//
// T can be any integral type <= 128 bits.
template <typename T>
void factorize(T x, std::vector<T>& factors)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<T>::digits <= 128, "");
    HPBC_PRECONDITION(x >= 2);  // 0 and 1 do not have prime factorizations
    factors.clear();

    namespace hd = ::hurchalla::detail;
    hd::impl_factorize::factorize_to_vector(
                                           x, factors, hd::PollardRhoIsPrime());
    HPBC_POSTCONDITION(factors.size() > 0);
    // The max possible vector size needed for factors is when all of them are 2
    constexpr int max_num_factors = ut_numeric_limits<T>::digits;
    HPBC_POSTCONDITION(factors.size() <= max_num_factors);
    // all the factors multiplied together should == x
    HPBC_POSTCONDITION(x == std::accumulate(factors.begin(), factors.end(),
                                      static_cast<T>(1), std::multiplies<T>()));
}


}  // end namespace

#endif
