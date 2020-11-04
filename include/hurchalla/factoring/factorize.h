// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_FACTORIZE_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORIZE_H_INCLUDED


#include "hurchalla/factoring/detail/impl_factorize.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <iterator>
#include <limits>
#include <array>
#include <vector>

namespace hurchalla { namespace factoring {


template <typename T>
std::array<T, hurchalla::modular_arithmetic::ma_numeric_limits<T>::digits>
factorize(T x, int& num_factors)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");

    // The max possible number of factors occurs when all factors equal 2
    constexpr auto max_num_factors = ma::ma_numeric_limits<T>::digits;
    std::array<T, max_num_factors> arr;

    struct FactorArrayAdapter {
        using value_type = T;
        explicit FactorArrayAdapter(std::array<T, max_num_factors>& a) :
                                                       arr(a), num_factors(0) {}
        void push_back(const value_type& val)
        {
            HPBC_ASSERT2(num_factors < max_num_factors);
            arr[num_factors] = val;
            ++num_factors;
        }
        std::size_t size() { return num_factors; }
    private:
        std::array<T, max_num_factors>& arr;
        std::size_t num_factors;
    };
    FactorArrayAdapter faa(arr);
    impl_factorize(std::back_inserter(faa), x);
    num_factors = static_cast<decltype(num_factors)>(faa.size());

    HPBC_POSTCONDITION(num_factors > 0);
    HPBC_POSTCONDITION(num_factors <= max_num_factors);
    // all the factors multiplied together should == x
    HPBC_POSTCONDITION(x == std::accumulate(arr.begin(),
             arr.begin()+num_factors, static_cast<T>(1), std::multiplies<T>()));
    return arr;
}


// As a general recommendation, use factorize().  However, if your stack size is
// severely limited, the following function may be preferable since it uses (the
// heap allocated) std::vector.  Note that if you used factorize(), the returned
// std::array for type T of uint32_t would take 128 bytes, uint64_t would take
// 512 bytes, and __uint128_t would take 2kb.
template <typename T>
std::vector<T> factorize_to_vector(T x)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");

    // The max possible vector size needed for factors is when all of them are 2
    constexpr auto max_num_factors = ma::ma_numeric_limits<T>::digits;
    std::vector<T> vec;
    vec.reserve(max_num_factors);
    impl_factorize(std::back_inserter(vec), x);

    HPBC_POSTCONDITION(vec.size() > 0);
    HPBC_POSTCONDITION(vec.size() <= max_num_factors);
    // all the factors multiplied together should == x
    HPBC_POSTCONDITION(x == std::accumulate(vec.begin(), vec.end(),
                                      static_cast<T>(1), std::multiplies<T>()));
    return vec;
}


}}  // end namespace

#endif
