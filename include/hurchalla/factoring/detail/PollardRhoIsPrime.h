// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/detail/ImplIsPrimeIntensive.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include <array>
#include <vector>

namespace hurchalla { namespace detail {


struct PollardRhoIsPrime {

template <typename MontType>
typename std::enable_if<
      (ut_numeric_limits<typename MontType::T_type>::digits == 128), bool>::type
operator()(const MontType& mf) const
{
    return is_prime_miller_rabin(mf);
}

template <typename MontType>
typename std::enable_if<
       (ut_numeric_limits<typename MontType::T_type>::digits == 64), bool>::type
operator()(const MontType& mf) const
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(mf.getModulus() > 1);
#if 0
    if (mf.getModulus() < UINT64_C(350269456337)) {
        constexpr std::size_t TRIAL_SIZE = 2;
        return is_prime_miller_rabin64_3_350269456337<TRIAL_SIZE>(mf);
    }
#endif
#if 1
    if (mf.getModulus() < (static_cast<T>(1) << 32)) {
        constexpr std::size_t TOTAL_BASES = 2;
        constexpr std::size_t TRIAL_SIZE = 2;
        return MillerRabinMontgomery
                          <MontType, 32, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
    }
#else
    static const ImplIsPrimeIntensive<std::uint32_t, true> isprime_sieve32;
    if (mf.getModulus() < (static_cast<T>(1) << 32)) {
        return isprime_sieve32(static_cast<std::uint32_t>(mf.getModulus()));
    }
#endif
    constexpr std::size_t TOTAL_BASES = 3;
    constexpr std::size_t TRIAL_SIZE = 2;
    return MillerRabinMontgomery
                          <MontType, 64, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}

template <typename MontType>
typename std::enable_if<
       (ut_numeric_limits<typename MontType::T_type>::digits == 32), bool>::type
operator()(const MontType& mf) const
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(mf.getModulus() > 1);
#if 1
    constexpr std::size_t TOTAL_BASES = 2;
    constexpr std::size_t TRIAL_SIZE = 2;
    return MillerRabinMontgomery<MontType, ut_numeric_limits<T>::digits,
                                 TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
#else
    static const ImplIsPrimeIntensive<std::uint32_t, true> isprime_sieve32;
    return isprime_sieve32(static_cast<std::uint32_t>(mf.getModulus()));
#endif
}

template <typename MontType>
typename std::enable_if<
       (ut_numeric_limits<typename MontType::T_type>::digits < 32), bool>::type
operator()(const MontType& mf) const
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(mf.getModulus() > 1);
    constexpr std::size_t TOTAL_BASES = 1;
    constexpr std::size_t TRIAL_SIZE = 1;
    return MillerRabinMontgomery<MontType, ut_numeric_limits<T>::digits,
                                 TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}

};


}} // end namespace

#endif
