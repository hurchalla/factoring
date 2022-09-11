// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include <type_traits>
#include <stdbool.h>
#include <cstdint>
#include <iostream>

#include "microecm_c.h"
#include "get_single_factor.h"
#include "hurchalla/factoring/is_prime.h"



template <typename T>
uint64_t test_factoring(T min, T max)
{
    uint64_t num_errors = 0;
    if (min >= hurchalla::ut_numeric_limits<T>::max())
        return 0;
    if (min % 2 == 0)
        min++;
    if (max <= 0)
        return 0;
    if (max % 2 == 0)
        max--;
    assert(min % 2 == 1 && max % 2 == 1);

    for (T x = min; x < max; x=x+2) {
        if (x < 2 || hurchalla::is_prime(x))
            continue;
#ifdef EXPECT_ARBITRARY_SIZE_FACTORS
        bool expect_arbitrary_factors = true;
#else
        bool expect_arbitrary_factors = false;
#endif
#ifdef USE_ECM_C_INTERFACE
        static_assert(std::is_same<T, uint64_t>::value);
        uint64_t random_val = 0;
        T ecm_result = getfactor_uecm(x, expect_arbitrary_factors, &random_val);
#else
        T ecm_result = hurchalla::get_single_factor_ecm(x,
                                                      expect_arbitrary_factors);
#endif
        if (ecm_result == 0 || x % ecm_result != 0 || ecm_result == x) {
            if constexpr (hurchalla::ut_numeric_limits<T>::digits <= 64) {
                num_errors++;
                std::cout << "Error: ECM, on value " << x << "\n";
            }
            else {
                num_errors++;
                std::cout << "Error: 128 bit ECM\n";
            }
        }
    }
#if 0   // #if 0  lets us skip Pollard-Rho testing
    for (T x = min; x < max; x=x+2) {
        if (x < 2 || hurchalla::is_prime(x))
            continue;
        T pr_result = hurchalla::get_single_factor_pollard_rho(x);
        if (pr_result == 0 || x % pr_result != 0 || pr_result == x) {
            if constexpr (hurchalla::ut_numeric_limits<T>::digits <= 64) {
                num_errors++;
                std::cout << "Error: pr on " << x << "\n";
            }
            else {
                num_errors++;
                std::cout << "Error: 128 bit pr\n";
            }
        }
    }
#endif
    return num_errors;
}

// Simple preliminary testing, for experimental verification

int main()
{
// test uint64_t
    std::cout << "\n";
    {
        uint64_t num_errors = 0;
        using T = uint64_t;
        T min = 0; T max = 100000;
        num_errors += test_factoring(min, max);

        max = hurchalla::ut_numeric_limits<T>::max();
        min = max - 100000;
        num_errors += test_factoring(min, max);

        std::cout << "uint64_t testing complete with " << num_errors << " errors\n";
    }

// test __uint128_t, if possible
#ifndef USE_ECM_C_INTERFACE   // the ECM C interface only supports uint64_t
#ifndef _MSC_VER   // MSVC doesn't support __uint128_t
    {
        uint64_t num_errors = 0;
        using T = __uint128_t;
        T min = 0; T max = 100000;
        num_errors += test_factoring(min, max);

        // ECM could do 128bit true max, but it's computationally infeasible for
        // pollard rho; so we use a much smaller max.
        max = static_cast<T>(1) << 70;
        min = max - 100000;
        num_errors += test_factoring(min, max);

        std::cout << "__uint128_t testing complete with " << num_errors << " errors\n";
    }
#endif
#endif

    std::cout << "\n";
}

