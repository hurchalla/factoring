// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IMPL_FACTORIZE_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_FACTORIZE_H_INCLUDED


#include "hurchalla/factoring/detail/factorize_dispatch.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstddef>
#include <array>
#include <vector>

namespace hurchalla { namespace detail {


template <typename T, class PrimalityFunctor>
std::array<T, ut_numeric_limits<T>::digits>
impl_factorize_to_array(
                     T x, int& num_factors, const PrimalityFunctor& is_prime_mf)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    // The max possible number of factors occurs when all factors equal 2
    constexpr std::size_t array_size = ut_numeric_limits<T>::digits;
    std::array<T, array_size> arr;

    struct FactorArrayAdapter {
        using value_type = T;
        explicit FactorArrayAdapter(std::array<T, array_size>& a) :
                                                       ar(a), factor_count(0) {}
        void push_back(const value_type& val)
        {
            HPBC_ASSERT2(factor_count < array_size);
            ar[factor_count] = val;
            ++factor_count;
        }
        std::size_t size() { return factor_count; }
    private:
        std::array<T, array_size>& ar;
        std::size_t factor_count;
    };
    FactorArrayAdapter faa(arr);
    detail::factorize_dispatch(std::back_inserter(faa), x, is_prime_mf);
    num_factors = static_cast<int>(faa.size());

    HPBC_POSTCONDITION(num_factors > 0);
    HPBC_POSTCONDITION(static_cast<std::size_t>(num_factors) <= array_size);
    return arr;
}


template <typename T, class PrimalityFunctor>
std::vector<T> impl_factorize_to_vector(T x, int max_num_factors,
                                            const PrimalityFunctor& is_prime_mf)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    std::vector<T> vec;
    vec.reserve(max_num_factors);
    detail::factorize_dispatch(std::back_inserter(vec), x, is_prime_mf);

    HPBC_POSTCONDITION(vec.size() > 0);
    HPBC_POSTCONDITION(vec.size() <= max_num_factors);
    return vec;
}


}} // end namespace

#endif
