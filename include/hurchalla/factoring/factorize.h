// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_FACTORIZE_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_H_INCLUDED


#include "hurchalla/factoring/detail/impl_factorize.h"
#include "hurchalla/factoring/detail/PollardRhoIsPrime.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstddef>
#include <array>
#include <vector>

namespace hurchalla {


// Returns an array that contains all factors of x, and writes the total number
// of factors to num_factors.  The array entries with index < num_factors are the
// factors.
template <typename T>
std::array<T, ut_numeric_limits<T>::digits>
factorize(T x, int& num_factors)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    using namespace hurchalla::detail;
    // The max possible number of factors occurs when all factors equal 2,
    // so ut_numeric_limits<T>::digits is sufficient to hold all factors.
    std::array<T, ut_numeric_limits<T>::digits> arr =
                   impl_factorize_to_array(x, num_factors, PollardRhoIsPrime());
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


// As a general recommendation, use factorize().  You may wish to use this
// function if your stack size is limited, since it uses (the heap allocated)
// std::vector.  Note that if you used factorize(), the returned std::array for
// type T of uint32_t would take 128 bytes, uint64_t would take 512 bytes, and
// __uint128_t would take 2kb.
//
// Returns a vector that contains all factors of x.
template <typename T>
std::vector<T> factorize_to_vector(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    using namespace hurchalla::detail;
    // The max possible vector size needed for factors is when all of them are 2
    constexpr int max_num_factors = ut_numeric_limits<T>::digits;
    std::vector<T> vec =
              impl_factorize_to_vector(x, max_num_factors, PollardRhoIsPrime());

    HPBC_POSTCONDITION(vec.size() > 0);
    HPBC_POSTCONDITION(vec.size() <= max_num_factors);
    // all the factors multiplied together should == x
    HPBC_POSTCONDITION(x == std::accumulate(vec.begin(), vec.end(),
                                      static_cast<T>(1), std::multiplies<T>()));
    return vec;
}


}  // end namespace

#endif
