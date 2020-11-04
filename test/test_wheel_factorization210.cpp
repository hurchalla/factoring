// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#include "factorize_bruteforce.h"
#include "hurchalla/factoring/detail/wheel_factorization210.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <cstddef>
#include <iterator>
#include <vector>
#include <array>
#include <algorithm>
#include <numeric>
#include <functional>

namespace {


#if 0
// A comparative performance test - we shouldn't normally run this.
// 6/7/20 haswell results: for uint32_t, template version ran about 2.5x slower.
//            changing both to uint64_t, template version ran about 5.5x slower.
//            (These results were roughly as expected)

    TEST(HurchallaFactoringIsPrime, perf_uint32_small_trial_div) {
        namespace hf = hurchalla::factoring;
        using T = std::uint32_t;
        T max = std::numeric_limits<T>::max();
        T dummy = 0;
        constexpr int factors_len = std::numeric_limits<T>::digits;
        T factors[factors_len];
        for (T x = max - 65336*32; x < max; ++x) {
            T y = x;
            dummy += (T)hf::small_trial_division(factors, factors_len, y);
        }
        EXPECT_TRUE(dummy != 0);
    }
    TEST(HurchallaFactoringIsPrime, perf_uint32_small_trial_div_template) {
        namespace hf = hurchalla::factoring;
        using T = std::uint32_t;
        T max = std::numeric_limits<T>::max();
        T dummy = 0;
        constexpr int factors_len = std::numeric_limits<T>::digits;
        T factors[factors_len];
        for (T x = max - 65336*32; x < max; ++x) {
            T y = x;
            dummy += (T)hf::small_trial_division<T>(factors, factors_len, y);
        }
        EXPECT_TRUE(dummy != 0);
    }
#endif



    TEST(HurchallaFactoringWheelFactorization210, exhaustive_uint16_t) {
        namespace hf = hurchalla::factoring;
        using T = std::uint16_t;
        for (T x = std::numeric_limits<T>::max(); x >= 2; --x) {
            std::vector<T> answer = hf::factorize_bruteforce(x);
            std::vector<T> factors;
            T q;
            hf::wheel_factorization210(std::back_inserter(factors), q, x);
            std::sort(answer.begin(), answer.end());
            std::sort(factors.begin(), factors.end());
            SCOPED_TRACE(testing::Message() << "x == " << x);
            EXPECT_TRUE(factors == answer);
        }
    }


    template <typename T>
    void test_factor_wheel210(T x, const std::vector<T>& answer)
    {
        namespace hf = hurchalla::factoring;
        namespace ma = hurchalla::modular_arithmetic;
        T q;

        // first test using std::vector
        std::vector<T> vec;
        hf::wheel_factorization210(std::back_inserter(vec), q, x);
        EXPECT_TRUE(q == 1);
        // at this time, I haven't made a guarantee for wheel_factorization210()
        // that the destination range will be sorted, so we'll sort it here.
        std::sort(vec.begin(), vec.end());
        EXPECT_TRUE(vec == answer);

        // second test using std::array
        // the max possible number of factors occurs when all factors equal 2
        constexpr auto max_num_factors = ma::ma_numeric_limits<T>::digits;
        std::array<T, max_num_factors> arr;
        struct FactorArrayAdapter {
            using value_type = T;
            explicit FactorArrayAdapter(std::array<T, max_num_factors>& a) :
                                                       arr(a), num_factors(0) {}
            void push_back(const value_type& val)
            {
                HPBC_ASSERT2(num_factors < max_num_factors);
                arr[num_factors] = val;
                ++num_factors;
            }
            std::size_t size() { return num_factors; }
        private:
            std::array<T, max_num_factors>& arr;
            std::size_t num_factors;
        };
        FactorArrayAdapter faa(arr);

        hf::wheel_factorization210(std::back_inserter(faa), q, x);
        auto num_factors = faa.size();
        EXPECT_TRUE(q == 1);
        EXPECT_TRUE(num_factors == answer.size());
        std::sort(arr.begin(), arr.begin()+num_factors);
        EXPECT_TRUE(std::equal(arr.begin(), arr.begin()+num_factors,
                                                               answer.begin()));
    }


    TEST(HurchallaFactoringWheelFactorization210, basic_tests) {
        using U = std::uint64_t;
        std::vector<U> answer = { 2, 3, 5, 13, 17 };
        // multiply all the factors in answer to get the number to factorize.
        U x = std::accumulate(answer.begin(), answer.end(), static_cast<U>(1),
                                                          std::multiplies<U>());
        SCOPED_TRACE(testing::Message() << "x == " << x);
        test_factor_wheel210<U>(x, answer);
    }

#ifdef __SIZEOF_INT128__
    TEST(HurchallaFactoringWheelFactorization210, basic_tests_128bit) {
        using U = __uint128_t;
        std::vector<U> answer = { 2, 3, 5, 13, 17 };
        // multiply all the factors in answer to get the number to factorize.
        U x = std::accumulate(answer.begin(), answer.end(), static_cast<U>(1),
                                                          std::multiplies<U>());
        test_factor_wheel210<U>(x, answer);
    }
#endif


} // end namespace
