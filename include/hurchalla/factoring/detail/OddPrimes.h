// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_ODDPRIMES_H_INCLUDED
#define HURCHALLA_FACTORING_ODDPRIMES_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_bruteforce.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/programming_by_contract.h"
#include <array>
#include <type_traits>
#include <cstdint>
#include <cstddef>

namespace hurchalla { namespace detail {


// Note: we use a struct with static functions in order to disallow ADL
struct OddPrimes {
private:
  // Returns a uint64_t std::array of the first N=SIZE odd primes.
  // Intended for use at compile-time, to assist get_array().
  template <int SIZE>
  static constexpr std::array<std::uint64_t, static_cast<std::size_t>(SIZE)>
  get_array64()
  {
    static_assert(SIZE > 0);
    using T = std::uint64_t;
    std::array<T, static_cast<std::size_t>(SIZE)> oddprimes64{};
    std::size_t i = 0;
    for (T x=3; i < static_cast<std::size_t>(SIZE); x=static_cast<T>(x+2)) {
        if (is_prime_bruteforce::call(x))
            oddprimes64[i++] = x;
        // Use a compile-time assert that the next loop iteration will not
        // overflow x.  Note: it should be impossible for any computationally
        // feasible use of this function to overflow x, since T is uint64_t.
        HPBC_CONSTEXPR_ASSERT(x != ut_numeric_limits<T>::max());
    }
    return oddprimes64;
  }
  // Intended for use at compile-time, to assist get_next_prime().
  template <typename T, T oddprime>
  static constexpr T get_next_prime_internal()
  {
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    static_assert(oddprime % 2 == 1);
    static_assert(oddprime != ut_numeric_limits<T>::max());
    T x = static_cast<T>(oddprime + 2);
    while (!is_prime_bruteforce::call(x)) {
        HPBC_CONSTEXPR_ASSERT(x != ut_numeric_limits<T>::max());
        x = static_cast<T>(x + 2);
    }
    return x;
  }

public:

  // Intended for use at compile-time, to initialize a constexpr variable.
  // Returns a std::array of the first N=SIZE odd primes.
  template <int SIZE>
  static constexpr auto get_array()
  {
    static_assert(SIZE > 0);
    constexpr std::array<std::uint64_t, static_cast<std::size_t>(SIZE)>
                                              oddprimes64 = get_array64<SIZE>();
    constexpr std::uint64_t lastprime = oddprimes64[SIZE - 1];
    // sanity check to be sure oddprimes64 is compile time initialized
    static_assert(oddprimes64[0] != 0);

    // get the smallest type that we can use to store all the primes
#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wtype-limits"
#endif
    using T = typename std::conditional<
                 (lastprime <= ut_numeric_limits<std::uint8_t>::max()),
                 std::uint8_t,
                 typename std::conditional<
                     (lastprime <= ut_numeric_limits<std::uint16_t>::max()),
                     std::uint16_t,
                     typename std::conditional<
                         (lastprime <= ut_numeric_limits<std::uint32_t>::max()),
                         std::uint32_t,
                         std::uint64_t
                     >::type
                 >::type
              >::type;
#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

    std::array<T, static_cast<std::size_t>(SIZE)> oddprimes{};
    for (std::size_t i = 0; i < static_cast<std::size_t>(SIZE); ++i)
        oddprimes[i] = static_cast<T>(oddprimes64[i]);
    return oddprimes;
  }


  // Intended for use at compile-time, to initialize a constexpr variable.
  // Returns the first prime larger than oddprime.
  template <typename T, T oddprime>
  static constexpr auto get_next_prime()
  {
    static_assert(ut_numeric_limits<T>::is_integer);
    static_assert(!ut_numeric_limits<T>::is_signed);
    static_assert(oddprime % 2 == 1);
    using T2 = typename sized_uint<ut_numeric_limits<T>::digits * 2>::type;
    constexpr T2 next = get_next_prime_internal<T2,static_cast<T2>(oddprime)>();
    // sanity check: ensure 'next' is compile time initialized
    static_assert(next != 0);
    using U = typename std::conditional<
                            (next <= ut_numeric_limits<T>::max()),
                             T, T2>::type;
    return static_cast<U>(next);
  }


  // Intended for use at compile-time, to initialize a constexpr variable.
  // Returns the square of a constant without overflow.
  template <typename T, T number>
  static constexpr auto get_constant_squared()
  {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // get a type that we can use to square the number without overflow.
    using T2 = typename std::conditional<
             (number < (static_cast<T>(1) << (ut_numeric_limits<T>::digits/2))),
             T,
             typename sized_uint<ut_numeric_limits<T>::digits * 2>::type
          >::type;
    using P2 = typename safely_promote_unsigned<T2>::type;
    constexpr P2 result = static_cast<P2>(number) * static_cast<P2>(number);
    // sanity check: ensure result is compile time initialized
    static_assert(result != 0 || number == 0);
    return static_cast<T2>(result);
  }
};


}}  // end namespace

#endif
