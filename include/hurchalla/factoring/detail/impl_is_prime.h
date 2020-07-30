// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/small_trial_division256.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/factoring/detail/FactorsContainerAdapter.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/modular_arithmetic/traits/extensible_make_unsigned.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>
#include <limits>
#include <array>
#include <vector>
#include <cstdint>

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace factoring {


template <typename T>
bool impl_is_prime(T x)
{
    HPBC_PRECONDITION2(x >= 0);
    namespace ma = hurchalla::montgomery_arithmetic;

    if (x < 2)  // handle the noncomposite yet nonprime numbers
        return false;

    // First try small trial divisions to find easy factors or easy primality.
    // If primality still unknown, use miller-rabin to prove prime or composite.

    {  // trial divisions
        using U= typename modular_arithmetic::extensible_make_unsigned<T>::type;

        // The max possible array length for factors is if all factors are 2
        constexpr size_t max_factors = std::numeric_limits<U>::digits;
#if 1 // use the stack (std::array).  Note: a 128 bit type T would take 2kb
        std::array<U, max_factors> factors;
#else // use the heap (std::vector)
        std::vector<U> factors;
#endif
        FactorsContainerAdapter<decltype(factors)> fca(factors);
        fca.reserve(max_factors);
        small_trial_division256(fca, static_cast<U>(x));
        if (fca.size() > 0)  //size>0 shows x is noncomposite, or we got factors
            return (fca.size() == 1);  // size 1 indicates x is noncomposite.

        // At this point we know fca.size() == 0, so we didn't find any factors
        // and couldn't detect if x is composite.  We'll fallback to determining
        // primality via the slower miller-rabin test below.
        // Note: we also know x >= 65536, since a postcondition of
        // small_trial_division256() is that if the number passed to it is
        // <= 65535, then the final fca.size() >= 1.  Conversely, if fca.size()
        // == 0, then the number passed in had to be > 65535.
    }

    using std::uint32_t; using std::uint64_t;

    // miller-rabin
    HPBC_ASSERT2(x % 2 == 1);
    if (x > std::numeric_limits<uint64_t>::max()) {
        return is_prime_miller_rabin(ma::MontgomeryForm<T>(x));
    } else if (x > std::numeric_limits<uint32_t>::max()) {
        if (std::is_same<T, int64_t>::value)
            return is_prime_miller_rabin(ma::MontgomeryForm<T>(x));
        else
            return is_prime_miller_rabin(ma::MontgomeryForm<uint64_t>(
                                                     static_cast<uint64_t>(x)));
    } else {
        if (std::numeric_limits<T>::digits < 32)
            return is_prime_miller_rabin(ma::MontgomeryForm<T>(x));
        else
            return is_prime_miller_rabin(ma::MontgomeryForm<uint32_t>(
                                                     static_cast<uint32_t>(x)));
    }
}


}}  // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
