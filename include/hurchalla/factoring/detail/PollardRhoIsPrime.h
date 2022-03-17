// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_POLLARD_RHO_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_POLLARD_RHO_IS_PRIME_H_INCLUDED


// This is simply a support Functor that lets Pollard-Rho factoring determine
// primality, in this case by selecting and calling the most suitable
// Miller-Rabin primality test.  (note that there is another Pollard-Rho
// primality support functor within factorize_intensive_uint32.h that is a
// generic lambda using Sieve of Eratosthenes as its primality test)


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace detail {


struct PollardRhoIsPrime {

template <typename MontType>
typename std::enable_if<
   (ut_numeric_limits<typename MontType::IntegerType>::digits==128), bool>::type
operator()(const MontType& mf) const
{
    return is_prime_miller_rabin::call(mf);
}

// We could have made this file much simpler if we had just called
// is_prime_miller_rabin(mf) for everything rather than using custom code below.
// - However -
// we choose to optimize the TRIAL_SIZE we use, under an assumption that the
// modulus we will test will be slightly more likely to be prime than composite.
// That is how factorize_pollard_rho() works - when factoring a number, in total
// its recursions will always test exactly one more prime number (the number is
// either a factor or the original number) than composite number for primality.

template <typename MontType>
typename std::enable_if<
  (ut_numeric_limits<typename MontType::IntegerType>::digits == 64), bool>::type
operator()(const MontType& mf) const
{
    using T = typename MontType::IntegerType;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(mf.getModulus() > 1);
#if HURCHALLA_TARGET_BIT_WIDTH >= 64
    // factorize_pollard_rho() should be the only user of this functor, and that
    // template function will always use a native type when possible (e.g. on a
    // 32 bit system it will switch from using a 64 bit type to instead use a 32
    // bit type if the modulus fits in the 32 bit type).  Thus the following
    // conditional would always be false for a 64 bit target system, but the
    // compiler won't know this and would still generate code for the clause. 
    // We use the preprocessor ifdef above to to avoid generating code for a
    // conditional that we know would always be false.  Even if we are wrong
    // about this (because something changed in the future most likely), we'll
    // still be okay because this conditional is used to run a more efficient
    // version of the tests further below.  The tests further below will still
    // produce the correct result if they are run instead.
    if (mf.getModulus() < (static_cast<T>(1) << 32)) {
        constexpr std::size_t TOTAL_BASES = 2;
        constexpr std::size_t TRIAL_SIZE = 2;
        return MillerRabinMontgomery
                          <MontType, 32, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
    }
#endif
    if (mf.getModulus() < UINT64_C(350269456337)) {
        constexpr std::size_t TRIAL_SIZE = 3;
        return is_prime_miller_rabin_special::
                                         case_350269456337_64_3<TRIAL_SIZE>(mf);
    }
    constexpr std::size_t TOTAL_BASES = 5;
    constexpr std::size_t TRIAL_SIZE = 3;
    return MillerRabinMontgomery
                          <MontType, 64, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}

template <typename MontType>
typename std::enable_if<
  (ut_numeric_limits<typename MontType::IntegerType>::digits == 32), bool>::type
operator()(const MontType& mf) const
{
    using T = typename MontType::IntegerType;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(mf.getModulus() > 1);
    constexpr std::size_t TOTAL_BASES = 2;
    constexpr std::size_t TRIAL_SIZE = 2;
    return MillerRabinMontgomery<MontType, ut_numeric_limits<T>::digits,
                                 TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}

template <typename MontType>
typename std::enable_if<
   (ut_numeric_limits<typename MontType::IntegerType>::digits < 32), bool>::type
operator()(const MontType& mf) const
{
    using T = typename MontType::IntegerType;
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
