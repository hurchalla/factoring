
#ifndef HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_IS_PRIME_H_INCLUDED


#include "hurchalla/factoring/detail/small_trial_division.h"
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/modular_arithmetic/traits/extensible_make_unsigned.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <type_traits>
#include <limits>

namespace hurchalla { namespace factoring {


template <typename T>
bool impl_is_prime(T x)
{
    HPBC_PRECONDITION2(x >= 0);

    // First try small trial divisions to find easy factors or easy primality.
    // If primality still unknown, use miller-rabin to prove prime or composite.

    {  // trial divisions
        if (x < 2)
            return false;
        namespace ma = hurchalla::modular_arithmetic;
        using U = typename ma::extensible_make_unsigned<T>::type;
        // The longest possible array of factors results when all factors are 2.
        // std::numeric_limits<U>::digits  gives this length.
        constexpr int factors_len = std::numeric_limits<U>::digits;
        auto factors_up = std::unique_ptr<U[]>(new U[factors_len]);
        T* factors = factors_up.get();
        // U factors[factors_len];   // alternative, on stack.
                                     // Static would be no good- not thread safe
        U u = static_cast<U>(x);
        int num_factors_found = small_trial_division(factors, factors_len, u);
        HPBC_ASSERT2(num_factors_found >= 0);
        x = static_cast<T>(u);
        if (num_factors_found > 0)
            return false;
        else if (x == 1)  // small_trial_division says there are no more factors
            return true;
    }

    namespace ma = hurchalla::montgomery_arithmetic;

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

#endif
