// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/factoring/detail/experimental/get_single_factor.h"
#include "hurchalla/factoring/is_prime.h"
#include <cstdint>
#include <iostream>


template <typename T>
void test_factoring(T min, T max)
{
    for (T x = min; x < max; ++x) {
        if (x < 2 || hurchalla::is_prime(x))
            continue;
        T ecm_result = hurchalla::get_single_factor_ecm(x);
        T pr_result = hurchalla::get_single_factor_pollard_rho(x);
        if (ecm_result == 0 || x % ecm_result != 0 || ecm_result == x) {
            if constexpr (hurchalla::ut_numeric_limits<T>::digits <= 64)
                std::cout << "Error: ecm on " << x << "\n";
            else
                std::cout << "Error: 128 bit ecm\n";
        }
        if (pr_result == 0 || x % pr_result != 0 || pr_result == x) {
            if constexpr (hurchalla::ut_numeric_limits<T>::digits <= 64)
                std::cout << "Error: pr on " << x << "\n";
            else
                std::cout << "Error: 128 bit pr\n";
        }
    }
}

// Simple preliminary testing, for experimental verification

int main()
{
    {
        using T = uint64_t;
        T min = 0; T max = 100000;
        test_factoring(min, max);
        max = hurchalla::ut_numeric_limits<T>::max();
        min = max - 100000;
        test_factoring(min, max);
        std::cout << "uint64_t done\n";
    }
    {
        using T = __uint128_t;
        T min = 0; T max = 100000;
        test_factoring(min, max);
        // ECM could do 128bit true max, but it's computationally infeasible for
        // pollard rho; so we use a much smaller max.
        max = static_cast<T>(1) << 70;
        min = max - 100000;
        test_factoring(min, max);
        std::cout << "__uint128_t done\n";
    }
}

