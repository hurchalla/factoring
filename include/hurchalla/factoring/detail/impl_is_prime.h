// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/is_prime_wheel210.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>
#include <array>
#include <vector>
#include <cstdint>
#include <cstddef>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace factoring {


template <typename T>
bool impl_is_prime(T x)
{
    HPBC_PRECONDITION2(x >= 0);
    namespace mont = hurchalla::montgomery_arithmetic;
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");
    using std::uint64_t; using std::uint32_t; using std::uint16_t;
    using std::size_t;

    // First try small trial divisions to find easy factors or easy primality.
    // If primality still unknown, use miller-rabin to prove prime or composite.

    // TODO: empirically find a decent value for max_trial_factor.  It will need
    // to fit in a type T variable.
    constexpr T max_trial_factor = 255;
    bool isSuccessful;
    bool isPrime = is_prime_wheel210(x, &isSuccessful, max_trial_factor);
    if (isSuccessful)
        return isPrime;

    // At this point, we didn't find any factors and we couldn't detect whether
    // x is prime.  We'll fall back to determining primality via miller-rabin.

    HPBC_ASSERT2(x % 2 == 1);
    if (x > ma::ma_numeric_limits<uint64_t>::max())
        return is_prime_miller_rabin(mont::MontgomeryForm<T>(x));
    else if (x > ma::ma_numeric_limits<uint32_t>::max())
        return is_prime_miller_rabin(mont::MontgomeryForm<uint64_t>(
                                                     static_cast<uint64_t>(x)));
    else if (x > ma::ma_numeric_limits<uint16_t>::max())
        return is_prime_miller_rabin(mont::MontgomeryForm<uint32_t>(
                                                     static_cast<uint32_t>(x)));
    else
        return is_prime_miller_rabin(mont::MontgomeryForm<T>(x));
}


}}  // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
