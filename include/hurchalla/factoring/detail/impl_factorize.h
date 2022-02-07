// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IMPL_FACTORIZE_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_FACTORIZE_H_INCLUDED


#include "hurchalla/factoring/detail/PrimeTrialDivisionWarren.h"
#include "hurchalla/factoring/detail/PrimeTrialDivisionMayer.h"
#include "hurchalla/factoring/detail/factorize_trialdivision.h"
#include "hurchalla/factoring/detail/factorize_pollard_rho.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstddef>
#include <array>
#include <vector>

namespace hurchalla { namespace detail {


namespace fz_detail {

#ifndef HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE
#  define HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionWarren
//#  define HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionMayer
#endif

#ifndef HURCHALLA_PR_TRIAL_DIVISION_SIZE
// FYI there are 54 primes below 256
// On quick benchmarks on Haswell, 135 worked well for PrimeTrialDivisionWarren
#  define HURCHALLA_PR_TRIAL_DIVISION_SIZE 135
#endif


template <class OutputIt, typename T, class PrimalityFunctor>
OutputIt
factorize_dispatch(OutputIt iter, T x, const PrimalityFunctor& is_prime_mf)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    static_assert(ut_numeric_limits<T>::digits % 2 == 0, "");
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations.
    // If x fails this precondition, we'll just write it to the destination -
    // technically the result is undefined with an unmet precondition.
    if (x < 2) {
        *iter++ = x;
        return iter;
    }
    constexpr T sqrtR = static_cast<T>(1)<<(ut_numeric_limits<T>::digits/2);

    T q, next_prime;
    iter = factorize_trialdivision<HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE,
                      HURCHALLA_PR_TRIAL_DIVISION_SIZE>(iter, q, next_prime, x);
    HPBC_ASSERT2(q >= 1);  // factorize_trialdivision() guarantees this
    if (q == 1)   // if factorize_trialdivision() completely factored x
        return iter;
    // factorize_trialdivision() guarantees that any factor of x that is less
    // than next_prime*next_prime must be prime.
    T threshold_always_prime = (next_prime < sqrtR) ?
          static_cast<T>(next_prime * next_prime) : ut_numeric_limits<T>::max();

    iter = factorize_pollard_rho(iter, q, is_prime_mf, threshold_always_prime);
    return iter;
}

} // end namespace fz_detail



template <typename T, class PrimalityFunctor>
std::array<T, ut_numeric_limits<T>::digits>
impl_factorize_to_array(T x, int& num_factors,
                        const PrimalityFunctor& is_prime_mf)
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
    fz_detail::factorize_dispatch(std::back_inserter(faa), x, is_prime_mf);
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
    fz_detail::factorize_dispatch(std::back_inserter(vec), x, is_prime_mf);

    HPBC_POSTCONDITION(vec.size() > 0);
    HPBC_POSTCONDITION(vec.size() <= max_num_factors);
    return vec;
}


}} // end namespace

#endif
