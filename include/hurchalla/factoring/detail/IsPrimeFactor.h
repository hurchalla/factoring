// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IS_PRIME_FACTOR_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_FACTOR_H_INCLUDED


// This is a support functor designed to help a factoring function determine
// primality.  It accomplishes this by selecting and calling the most suitable
// Miller-Rabin primality test.  (note that there is another factoring primality
// support functor within factorize_intensive_uint32.h that is a generic lambda
// using the Sieve of Eratosthenes for its primality test)
//    We could have made this file extremely simple by just delegating to
// is_prime_miller_rabin::call(mf) for everything..
//    - However -
// we choose to optimize the TRIAL_SIZE we use, under an assumption that the
// modulus we will test will be slightly more likely to be prime than composite.
// That is what occurs in FactorizeStage2() - when factoring a number, in total
// its recursions will always test exactly one more prime number than composite
// number for primality.


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <type_traits>

namespace hurchalla { namespace detail {


struct IsPrimeFactor {

template <typename MontType>
typename std::enable_if<
   (ut_numeric_limits<typename MontType::IntegerType>::digits==128), bool>::type
operator()(const MontType& mf) const
{
    return is_prime_miller_rabin::call_mont(mf);
}


template <typename MontType>
typename std::enable_if<
  (ut_numeric_limits<typename MontType::IntegerType>::digits == 64), bool>::type
operator()(const MontType& mf) const
{
    using T = typename MontType::IntegerType;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(mf.getModulus() > 1);

    if constexpr (HURCHALLA_TARGET_BIT_WIDTH >= 64) {
        // We expect the factorization functions to be the only callers of this
        // functor, and they will always use a native type when possible.  E.g.
        // on a 32 bit system they will switch from using a 64 bit type to using
        // a 32 bit type if the modulus fits in the 32 bit type, and so when
        // they call this functor with a 64 bit type we already know that the
        // modulus does not fit in a 32 bit type.  Thus the conditional below
        // would always be false for a 32 bit target system, but the compiler
        // has no way to know this and so it would still generate code for the
        // condition's clause.  This is why we use the "constexpr if" above,
        // which will prevent the compiler from generating the useless code.
        //    Even if we are wrong about this (because something changed in the
        // future most likely), we'll still be okay because this conditional is
        // used to run a more efficient version of the tests further below.  The
        // tests further below will still produce the correct result if they are
        // run instead.
        static_assert(ut_numeric_limits<T>::digits == 64);
        if (mf.getModulus() < (static_cast<T>(1) << 32)) {
            constexpr std::size_t TOTAL_BASES = 2;
            constexpr std::size_t TRIAL_SIZE = 2;
            return MillerRabinMontgomery
                              <MontType, 32, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
        }
    }

    if (mf.getModulus() < (static_cast<std::uint64_t>(1) << 44)) {
        constexpr std::size_t TOTAL_BASES = 3;
        constexpr std::size_t TRIAL_SIZE = 3;
        return MillerRabinMontgomery
                          <MontType, 44, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
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
