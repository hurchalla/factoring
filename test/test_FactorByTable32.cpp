// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#include "hurchalla/factoring/resource_intensive_api/FactorByTable32.h"
#include "hurchalla/factoring/resource_intensive_api/IsPrimeIntensive.h"
#include "hurchalla/factoring/resource_intensive_api/factorize_intensive_uint32.h"
#include "hurchalla/factoring/factorize.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"

#include "gtest/gtest.h"

#include <cstdint>
#include <limits>
#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <filesystem>
#include <string>
#include <optional>
#include <chrono>
#include <iostream>
#include <random>
#include <unordered_set>
#include <fstream>

namespace {

using namespace hurchalla;
using namespace std::chrono;
using dsec = duration<double>;


template <typename T>
T calculate_x(const std::vector<T>& answer)
{
    return std::accumulate(answer.begin(), answer.end(), static_cast<T>(1),
                                                      std::multiplies<T>());
}

template <typename T, int BITLEN, bool FAVOR_SMALL_SIZE>
void test_factorize(const std::vector<T>& answer,
                    const FactorByTable<BITLEN, FAVOR_SMALL_SIZE>& factorTable)
{
    // multiply all the factors in answer to get the number to factorize.
    T x = calculate_x(answer);
    HPBC_PRECONDITION2(x < (static_cast<uint64_t>(1) << BITLEN));

    int num_factors;
    auto arr = factorTable(x, num_factors);
    EXPECT_TRUE(num_factors == static_cast<int>(answer.size()));
    // at this time, I haven't made a guarantee for FactorByTable32()
    // that the destination range will be sorted, so we'll sort it here.
    std::sort(arr.begin(), arr.begin()+num_factors);
    EXPECT_TRUE(std::equal(arr.begin(), arr.begin()+num_factors,
                                                           answer.begin()));
}

template <int BITLEN, bool FAVOR_SMALL_SIZE>
void test_all_valid_inputs(const FactorByTable<BITLEN, FAVOR_SMALL_SIZE>& factorTable)
{
    IsPrimeIntensive<uint32_t, true> is_prime;
    constexpr uint32_t maxvalid = static_cast<uint32_t>(
                                      (static_cast<uint64_t>(1) << BITLEN) - 1);
    for (uint32_t x = maxvalid; x >= 2; --x)
    {
        int numfactors;
        auto arr = factorTable(x, numfactors);
        uint32_t product = 1;
        for (int i=0; i<numfactors; ++i) {
            EXPECT_TRUE(arr[i] > 1);
            EXPECT_TRUE(is_prime(arr[i]));
            product = product * arr[i];
        }
        EXPECT_TRUE(product == x);
    }
}


template <int BITLEN, bool FAVOR_SMALL_SIZE>
dsec quick_bench(const FactorByTable<BITLEN, FAVOR_SMALL_SIZE>& factorTable,
                 uint32_t min, uint32_t max, uint32_t samplesize)
{
    HPBC_PRECONDITION2(min >= 2);
    HPBC_PRECONDITION2(max > min);
    // Next line to avoid perf problems, by ensuring at worst no more
    // than half of the random numbers we generate will already exist
    // in our table (of numbers to reject).
    // Assert the samplesize is less than half the range size.
    HPBC_PRECONDITION2(samplesize < (max - min + 1)/2);

    auto t0 = steady_clock::now();
    std::vector<uint32_t> randomvec;
    {
        // static so that we start the program with the same seed, but don't
        // keep resetting to that same seed every time we call this function
        static std::mt19937 rng;
        std::uniform_int_distribution<uint32_t> dist(min, max);
        using hashtable = std::unordered_set<uint32_t>;
        hashtable ht;
        for (uint32_t i = 0; i < samplesize; ++i) {
            uint32_t val;
            while (true) {
                val = dist(rng);
                while (val >= min && (val % 2 == 0))
                    --val;
                if (val < min)
                    continue;
                std::pair<hashtable::iterator, bool> pair = ht.insert(val);
                bool insert_succeeded = pair.second;
                if (insert_succeeded)  // if val wasn't already in ht
                    break;
            }
            randomvec.push_back(val);
        }
    }
    auto t1 = steady_clock::now();
//    std::cout << "randomvec creation time " << dsec(t1-t0).count() << "\n";

    bool impossible_happened = false;

    t0 = steady_clock::now();
    for (uint32_t x : randomvec) {
       if (x == 0) {
            impossible_happened = true;
            break;
        }
    }
    t1 = steady_clock::now();
//    std::cout << "vector access time " << dsec(t1-t0).count() << "\n";


// can enable this section to also bench the other factorization functions
#if 0
    t0 = steady_clock::now();
    for (uint32_t x : randomvec) {
        int num_factors;
        auto arr = factorize(x, num_factors);
        // We need to prevent the compiler from completely removing
        // the factorize calls due to the array never being used.
        // So we'll check arr[0] (which is never 0) just so it's used.
        if (arr[0] == 0) {
            impossible_happened = true;
            break;
        }
    }
    t1 = steady_clock::now();
    std::cout << "factorize time " << dsec(t1-t0).count() << "\n";

    IsPrimeIntensive<uint32_t, true> is_prime;
    t0 = steady_clock::now();
    for (uint32_t x : randomvec) {
        int num_factors;
        auto arr = factorize_intensive_uint32(x, num_factors, is_prime);
        if (arr[0] == 0) {
            impossible_happened = true;
            break;
        }
    }
    t1 = steady_clock::now();
    std::cout << "factorize_intensive_uint32 time " << dsec(t1-t0).count() << "\n";
#endif


    t0 = steady_clock::now();
    for (uint32_t x : randomvec) {
        int num_factors;
        auto arr = factorTable(x, num_factors);
        if (arr[0] == 0) {
            impossible_happened = true;
            break;
        }
    }
    t1 = steady_clock::now();
//    std::cout << "   FactorByTable time " << dsec(t1-t0).count() << "\n";

    if (impossible_happened)
        std::cout << "impossible\n";

    return dsec(t1-t0);
}



bool file_exists(const std::string& filepath)
{
    std::ifstream fs(filepath);
    return fs.is_open();  // could also check nothing went wrong via fs.good()
}


template <int BITLEN, bool FAVOR_SMALL_SIZE>
void basic_tests_bit_limited(bool test_all_inputs, bool run_benchmark_32bit)
{
    using Path = std::filesystem::path;
    Path path = std::filesystem::current_path();
    std::string pathstr = path.string();
#ifdef _WIN32
    std::string separator("\\");
#else
    std::string separator("/");
#endif
    std::string size_desciption = (FAVOR_SMALL_SIZE) ? "_smaller_" : "_bigger_";
    std::string filepath = pathstr + separator + "factor_table" +
                              size_desciption + std::to_string(BITLEN) + ".bin";

    // we use std::optional as a workaround to allow us to construct a
    // FactorByTable within an assertion, and then have the factortable
    // still be in scope after the assertion.
    std::optional<FactorByTable<BITLEN, FAVOR_SMALL_SIZE>> opt;
    if (!file_exists(filepath)) {
        std::cout << "  Constructing a factor table... this may take a few minutes...\n";
        ASSERT_NO_THROW(opt.emplace());
        ASSERT_NO_THROW(opt->writeTableToFile(filepath.c_str()));
    } else {
        ASSERT_NO_THROW(opt.emplace(filepath.c_str(), false));
    }

    FactorByTable<BITLEN, FAVOR_SMALL_SIZE> factorTable(std::move(*opt));

    using U = std::uint32_t;

    // basic tests
    static_assert(BITLEN >= 13);
    std::vector<U> answer1a = { 2, 3, 5, 13, 17 };
    SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer1a));
    test_factorize<U>(answer1a, factorTable);

    static_assert(BITLEN >= 16);
    std::vector<U> answer1b = { 241, 251 };
    SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer1b));
    test_factorize<U>(answer1b, factorTable);

    // hard semiprimes
    if constexpr (BITLEN >= 32) {
        U twoPow16 = static_cast<U>(1) << 16;
        // use largest primes < (1<<16):
        // (1<<16) minus { 15, 17, 39, 57, 87, 89, 99, 113, 117, 123 }
        std::vector<U> answer2 = { twoPow16 - 17, twoPow16 - 15 };
        SCOPED_TRACE(testing::Message() << "x == " << calculate_x(answer2));
        test_factorize(answer2, factorTable);
    }

    // test_all_valid_inputs() takes a long time to run with a 32 bit table
    // (~30 mins), so ordinarily we would want to have test_all_inputs == false.
    // Even with a 24 bit table, it takes minutes in debug mode (fyi in
    // release mode 24bit is very quick though)
    if (test_all_inputs)
        test_all_valid_inputs(factorTable);

    if constexpr (BITLEN == 32) {
        if (run_benchmark_32bit) {
            dsec::rep besttime = 0;
            for (int i=0; i<10; ++i) {
                uint32_t max = std::numeric_limits<uint32_t>::max();
                uint32_t min = max/2;
                uint32_t samplesize = 4000000;
                dsec tmp = quick_bench(factorTable, min, max, samplesize);
                if (besttime == 0 || tmp.count() < besttime)
                    besttime = tmp.count();
            }
            std::cout << "best bench time " << besttime << "\n";
        }
    }
}



// 24 bit tables are very fast to create, so these tests can be part of our
// routine unit testing.  In contrast, creating (the default) 32 bit factor
// tables from scratch takes too long for everyday unit testing.
TEST(HurchallaFactorByTable32, basic_tests_24bit_limited_smaller) {
    constexpr int BITLIMIT = 24;
    constexpr bool FAVOR_SMALL_SIZE = true;
#if ((defined(_MSC_VER) && !defined(_DEBUG)) || defined(__OPTIMIZE__))
    bool test_all_inputs = false;
#else
    bool test_all_inputs = false;
#endif
    bool run_benchmark_32bit = false;
    basic_tests_bit_limited<BITLIMIT, FAVOR_SMALL_SIZE>(test_all_inputs, run_benchmark_32bit);
}
TEST(HurchallaFactorByTable32, basic_tests_24bit_limited_bigger) {
    constexpr int BITLIMIT = 24;
    constexpr bool FAVOR_SMALL_SIZE = false;
#if ((defined(_MSC_VER) && !defined(_DEBUG)) || defined(__OPTIMIZE__))
    bool test_all_inputs = false;
#else
    bool test_all_inputs = false;
#endif
    bool run_benchmark_32bit = false;
    basic_tests_bit_limited<BITLIMIT, FAVOR_SMALL_SIZE>(test_all_inputs, run_benchmark_32bit);
}


// This section is disabled for automatic everyday unit testing, since we
// generally assume the factor table doesn't already exist on disk for us to
// read in, and creating 32 bit factor tables from scratch for these tests
// typically takes a few minutes for each test.
#if 0

// the following tests take a while without optimization, so if we're building
// without any optimizations we usually want to skip them
#if ((defined(_MSC_VER) && !defined(_DEBUG)) || defined(__OPTIMIZE__))
TEST(HurchallaFactorByTable32, basic_tests_table_32bit_smaller) {
    constexpr int BITLIMIT = 32;
    constexpr bool FAVOR_SMALL_SIZE = true;
    bool test_all_inputs = false;
    bool run_benchmark_32bit = true;
    basic_tests_bit_limited<BITLIMIT, FAVOR_SMALL_SIZE>(test_all_inputs, run_benchmark_32bit);
}
TEST(HurchallaFactorByTable32, basic_tests_table_32bit_bigger) {
    constexpr int BITLIMIT = 32;
    constexpr bool FAVOR_SMALL_SIZE = false;
    bool test_all_inputs = false;
    bool run_benchmark_32bit = true;
    basic_tests_bit_limited<BITLIMIT, FAVOR_SMALL_SIZE>(test_all_inputs, run_benchmark_32bit);
}
#endif

#endif

 
} // end namespace
