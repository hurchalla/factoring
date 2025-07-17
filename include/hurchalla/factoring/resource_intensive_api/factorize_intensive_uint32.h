// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FACTORIZE_INTENSIVE32_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_INTENSIVE32_H_INCLUDED


#if !defined(HURCHALLA_FACTORING_DISALLOW_INLINE_ASM) && \
        !defined(HURCHALLA_ALLOW_INLINE_ASM_ALL)
#  define HURCHALLA_ALLOW_INLINE_ASM_ALL
#endif

#include "hurchalla/factoring/resource_intensive_api/IsPrimeTable.h"
#include "hurchalla/factoring/detail/impl_factorize.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstddef>
#include <cstdint>
#include <array>
#include <numeric>
#include <functional>

namespace hurchalla {


// Returns an array that contains all factors of x, and writes the total number
// of factors to num_factors.  The array entries with index < num_factors are
// the factors.
// The argument 'expect_arbitrary_size_factors' does not affect the results, but
// it will optimize the factoring to be faster when you know or expect all the
// factors will be large (assuming you set it to false), or it will optimize the
// factoring for arbitrary size factors (if set to true, which is the default).
//
// This function is generally recommended only for the case where you have
// already created the IsPrimeTable object for other uses (presumably for
// primality testing), or for the unusual case where you need to maximize
// performance of factoring an extremely large collection of 32 bit integers.
// Prior to calling this function you will need to create the IsPrimeTable
// object 'ipt'.  'Ipt' has very high overhead - it needs ~256MB of memory, and
// it takes a few seconds to construct.  If you can disregard the considerable
// contruction time (and memory use) of ipt, this function will almost always be
// faster on a given system than factorize() for uint32_t values.
inline std::array<std::uint32_t, 32>
factorize_intensive_uint32(std::uint32_t x, unsigned int& num_factors,
                           const IsPrimeTable<std::uint32_t>& ipt,
                           bool expect_arbitrary_size_factors = true)
{
    using T = std::uint32_t;
    static_assert(ut_numeric_limits<T>::digits == 32);

    using namespace hurchalla::detail;
    // The max possible number of factors occurs when all factors equal 2,
    // so 32 is sufficient to hold all factors of a uint32_t number.
    std::array<T, 32> arr = impl_factorize::factorize_to_array(x, num_factors,
                      [&ipt](const auto& mf) { return ipt(mf.getModulus()); },
                      expect_arbitrary_size_factors);
    // After calling this function, a client should not index the returned
    // array with an index >= num_factors.  As a defensive measure, we'll set
    // all array entries at or beyond num_factors to 0 - this may help to make
    // an indexing error more obvious if a caller later makes this mistake.
    for (auto i = num_factors; i < arr.size(); ++i)
        arr[i] = 0;

    HPBC_POSTCONDITION(num_factors > 0);
    HPBC_POSTCONDITION(num_factors <= arr.size());
    // all the factors multiplied together should == x
    HPBC_POSTCONDITION(x == std::accumulate(arr.begin(),
             arr.begin()+num_factors, static_cast<T>(1), std::multiplies<T>()));
    return arr;
}


}  // end namespace

#endif
