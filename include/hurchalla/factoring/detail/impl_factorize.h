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
#include "hurchalla/factoring/detail/FactorizeStage2.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <cstddef>
#include <cstdint>
#include <array>
#include <vector>

namespace hurchalla { namespace detail {


// Do *NOT* change the values for the macro in the code immediately below.  If
// you wish to use a different value for the macro (which is fine), please
// predefine the macro when compiling.
#ifndef HURCHALLA_FACTORING_ECM_THRESHOLD_BITS
#  define HURCHALLA_FACTORING_ECM_THRESHOLD_BITS 34
#endif


#ifndef HURCHALLA_TRIAL_DIVISION_TEMPLATE
#  define HURCHALLA_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionWarren
//#  define HURCHALLA_TRIAL_DIVISION_TEMPLATE PrimeTrialDivisionMayer
#endif

#ifndef HURCHALLA_TRIAL_DIVISION_SIZE
// FYI there are 54 primes below 256
// On quick benchmarks on Haswell, 139 worked well for PrimeTrialDivisionWarren
#  define HURCHALLA_TRIAL_DIVISION_SIZE 139
#endif



// Note: we use a struct with static functions in order to disallow ADL
struct impl_factorize {
private:

  template <int EcmMinBits, int MaxBitsX, class OutputIt,
            typename T, class PrimalityFunctor>
  static OutputIt
  dispatch(OutputIt iter, T x, const PrimalityFunctor& is_prime_functor,
           bool expect_arbitrary_size_factors)
  {
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    static_assert(ut_numeric_limits<T>::digits % 2 == 0);
    HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations.
    // If x fails this precondition, we'll just write it to the destination -
    // technically the result is undefined with an unmet precondition.
    if (x < 2) {
        *iter++ = x;
        return iter;
    }
    constexpr T sqrtR = static_cast<T>(1)<<(ut_numeric_limits<T>::digits/2);

    T q, next_prime;
    if (expect_arbitrary_size_factors) {
        iter = factorize_trialdivision::call<HURCHALLA_TRIAL_DIVISION_TEMPLATE,
                         HURCHALLA_TRIAL_DIVISION_SIZE>(iter, q, next_prime, x);
    } else {
        // Since we don't expect arbitrary size factors, we can assume it would
        // be a waste of time to look for small factors (via trial division).
        // However, to guarantee correctness we do need to attempt to extract
        // any and all occurences of the factor 2, so that we can later satisfy
        // Montgomery arithmetic's precondition that its modulus must be odd.
        q = x;
        while (q % 2 == 0) {
            *iter++ = 2;
            q = static_cast<T>(q / 2);
        }
        next_prime = 3;
    }

    HPBC_ASSERT2(q >= 1);  // factorize_trialdivision() guarantees this
    if (q == 1)   // if factorize_trialdivision() completely factored x
        return iter;
    // factorize_trialdivision() guarantees that any factor of q that is less
    // than next_prime*next_prime must be prime.
    T always_prime_limit = (next_prime < sqrtR) ?
          static_cast<T>(next_prime * next_prime) : ut_numeric_limits<T>::max();

    FactorizeStage2<EcmMinBits, MaxBitsX, T>
                          factorize_stage2(always_prime_limit,
                                           expect_arbitrary_size_factors);
    iter = factorize_stage2(iter, is_prime_functor, q);
    return iter;
  }


public:

  template <int EcmMinBits = HURCHALLA_FACTORING_ECM_THRESHOLD_BITS,
            typename T, class PrimalityFunctor>
  static std::array<T, ut_numeric_limits<T>::digits>
  factorize_to_array(T x, unsigned int& num_factors,
                     const PrimalityFunctor& is_prime_functor,
                     bool expect_arbitrary_size_factors)
  {
    static_assert(EcmMinBits > 0);
    static_assert(ut_numeric_limits<T>::is_integer);

    using U = typename extensible_make_unsigned<T>::type;

    // The max possible number of factors occurs when all factors equal 2
    constexpr std::size_t array_size = ut_numeric_limits<T>::digits;
    std::array<T, array_size> arr;

    struct FactorArrayAdapter {
        using value_type = U;
        explicit FactorArrayAdapter(std::array<T, array_size>& a) :
                                                       ar(a), factor_count(0) {}
        void push_back(const value_type& val)
        {
            HPBC_ASSERT2(factor_count < array_size);
            HPBC_ASSERT2(val <= static_cast<U>(ut_numeric_limits<T>::max()));
            ar[factor_count] = static_cast<T>(val);
            ++factor_count;
        }
        std::size_t size() { return factor_count; }
    private:
        std::array<T, array_size>& ar;
        std::size_t factor_count;
    };
    FactorArrayAdapter faa(arr);
    // MaxBitsX lets the called function know the compile-time limit to the
    // possible range of x, so that it can use an efficient Monty type.
    constexpr int MaxBitsX = ut_numeric_limits<T>::digits;
    dispatch<EcmMinBits, MaxBitsX>(std::back_inserter(faa), static_cast<U>(x),
                               is_prime_functor, expect_arbitrary_size_factors);
    num_factors = static_cast<unsigned int>(faa.size());

    HPBC_POSTCONDITION(num_factors > 0);
    HPBC_POSTCONDITION(num_factors <= array_size);
    return arr;
  }


  template <int EcmMinBits = HURCHALLA_FACTORING_ECM_THRESHOLD_BITS,
            typename T, class PrimalityFunctor>
  static void factorize_to_vector(T x, std::vector<T>& vec,
                                  const PrimalityFunctor& is_prime_functor,
                                  bool expect_arbitrary_size_factors)
  {
    static_assert(EcmMinBits > 0);
    static_assert(ut_numeric_limits<T>::is_integer);
    using U = typename extensible_make_unsigned<T>::type;

    // The max possible vector size needed for factors is when all of them are 2
    constexpr int max_num_factors = ut_numeric_limits<T>::digits;
    vec.reserve(max_num_factors);

    // convenient adapter to correctly cast from unsigned to (possibly)signed T
    struct FactorVectorAdapter {
        using value_type = U;
        explicit FactorVectorAdapter(std::vector<T>& a) : v(a) {}
        void push_back(const value_type& val)
        {
            HPBC_ASSERT2(val <= static_cast<U>(ut_numeric_limits<T>::max()));
            v.push_back(static_cast<T>(val));
        }
    private:
        std::vector<T>& v;
    };
    FactorVectorAdapter fva(vec);
    // MaxBitsX lets the called function know the compile-time limit to the
    // possible range of x, so that it can use an efficient Monty type.
    constexpr int MaxBitsX = ut_numeric_limits<T>::digits;
    dispatch<EcmMinBits, MaxBitsX>(std::back_inserter(fva), static_cast<U>(x),
                               is_prime_functor, expect_arbitrary_size_factors);
    HPBC_POSTCONDITION(vec.size() > 0);
    HPBC_POSTCONDITION(vec.size() <= max_num_factors);
  }

}; // end struct impl_factorize


}} // end namespace

#endif
