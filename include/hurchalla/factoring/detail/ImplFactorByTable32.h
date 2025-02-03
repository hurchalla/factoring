// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_IMPL_FACTOR_BY_TABLE32_H_INCLUDED
#define HURCHALLA_FACTORING_IMPL_FACTOR_BY_TABLE32_H_INCLUDED

#include "hurchalla/factoring/factorize.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/BitpackedUintVector.h"
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <exception>
#include <type_traits>
#include <limits>

namespace hurchalla { namespace detail {


// we use DUMMY to (hopefully) prevent the enclosed static array from using any
// exe space or memory when this class isn't used.
template <typename DUMMY = void>
struct PrimesUnder65536
{
    static_assert(std::is_same<DUMMY, void>::value);

    // Ordered array of the 6542 primes that are less than 2^16.
    static constexpr std::uint16_t prime[] = {
        #include "primes_under_65536.include"
    };
    static constexpr int num_primes = 6542;
    static_assert(sizeof(prime)/sizeof(prime[0]) == num_primes);
};



template <int INPUT_BIT_LIMIT, bool FAVOR_SMALL_SIZE>
class ImplFactorByTable32
{
private:
    // Instantiating with FAVOR_SMALL_SIZE = true saves ~15% memory/disk space,
    // at a cost of the factor table's entries not being a power of 2 in their
    // bitsize (which takes longer to get/set).  Less significantly, it also
    // costs an additional memory access (to the buffer in PrimesUnder65536)
    // when factoring.
    // When FAVOR_SMALL_SIZE == true, the expected memory/disk space use (if
    // wheel_divisors[] contains all primes up to 13) is ~1.4GB.  When
    // FAVOR_SMALL_SIZE == false, expected use is ~1.6GB.
    // Initial speed testing suggests factoring is a little over 10% faster on
    // x86 when choosing FAVOR_SMALL_SIZE = false, rather than true.

    // wheel_divisors must contain the number 2.  In theory it could contain all
    // primes up to 19 to save memory/disk space, but using only up to 13 is
    // recommended:
    // using all primes up to 19 would require ~20MB of exe space for our static
    // constexpr array wheel_reindex; indeed this would be a favorable tradeoff
    // of increased wheel_reindex size for reduced factor table size, but most
    // compilers are unlikely to be able to create such a large wheel_reindex at
    // (constexpr) compile time.  In practice, using all primes up to 17 (which
    // needs ~1MB of exe space for wheel_reindex) is a likely limit, and even
    // that can be too large for many compilers.
    // Using all primes up to 13 is the more practical option (it requires ~30k
    // for wheel_reindex).  Or you can very safely use fewer primes if you want
    // or need; using all primes up to 11, or all primes up to 7, is fine.
    // If you really want to use all primes up to 19, you will likely need to
    // create all the entries for wheel_reindex in a separate program that
    // writes the entries to file, and then you'll need to #include that file
    // where you define wheel_reindex in this class.  See the definition of
    // PrimesUnder65536 for a #include example.
    static constexpr std::uint8_t wheel_divisors[] = { 2, 3, 5, 7, 11, 13 };

    // On some systems, using a few extra prime divisors may or may not improve
    // factoring speed of arbitrary numbers (they often contain small factors),
    // but you would need to benchmark to know.  Regardless, using extra prime
    // divisors will almost always be unnecessary overhead if you know you will
    // be factoring 'hard' numbers, since the nature of 'hard to factor numbers'
    // is that they have no small factors.
    // If we expect extra divisors to slow down factoring, we can use an empty
    // extra_divisors array.  An empty array is my default recommendation.
    //   It's illegal for a C style array to have zero elements, so we use
    // std::array; unfortunately we have to explicitly specify the size.
#if 1
    static constexpr std::array<std::uint8_t, 0> extra_divisors{};
#elif 1
    static constexpr std::array<std::uint8_t, 3> extra_divisors{ 17, 19, 23 };
#else
    static constexpr std::array<std::uint8_t, 9> extra_divisors{ 17, 19, 23, 29, 31, 37, 41, 43, 47 };
#endif


// ----------- Usually you should change nothing past this point. --------------


    static constexpr std::size_t size_wheel_divs = sizeof(wheel_divisors) / sizeof(wheel_divisors[0]);

    // We work around an apparent bug in clang v16 to v18 (and perhaps beyond)
    // where wheel_divisors[0] can't be accessed at compile time - clang states
    // in a compilation output note that "read of element of array without known
    // bound is not allowed in a constant expression".  Obviously the bound for
    // wheel_divisors could be determined during compilation, and so this seems
    // to be a clang bug.  We use copy_array() defined here to work around this
    // issue.
    template <typename T, unsigned int N>
    static constexpr std::array<std::remove_cv_t<T>, N> copy_array(T (&arr)[N])
    {
        std::array<std::remove_cv_t<T>, N> tmp{};
        for (unsigned int i=0; i<N; ++i)
            tmp[i] = arr[i];
        return tmp;
    }
    static constexpr std::array<std::uint8_t, size_wheel_divs>
                                        wheel_divs = copy_array(wheel_divisors);


    // multiple places in this file rely on 2 being a wheel divisor
    static_assert(wheel_divs[0] == 2);
    // wheel_divs must be prime, and it saves no memory space to go past 19
    static_assert(size_wheel_divs > 0);
    static_assert(wheel_divs[0] == 2);
    static_assert(size_wheel_divs <= 1 || wheel_divs[1] == 3);
    static_assert(size_wheel_divs <= 2 || wheel_divs[2] == 5);
    static_assert(size_wheel_divs <= 3 || wheel_divs[3] == 7);
    static_assert(size_wheel_divs <= 4 || wheel_divs[4] == 11);
    static_assert(size_wheel_divs <= 5 || wheel_divs[5] == 13);
    static_assert(size_wheel_divs <= 6 || wheel_divs[6] == 17);
    static_assert(size_wheel_divs <= 7 || wheel_divs[7] == 19);
    static_assert(size_wheel_divs < 8); // there's little/no benefit past 19


    // use DUMMY to work around an undefined class error when initializing the
    // constexpr variable 'wheel_size' below, likely due to the class being
    // incomplete at constexpr init time.
    template <class DUMMY = void>
    static constexpr std::uint64_t get_wheel_size()
    {
        static_assert(size_wheel_divs > 0);
        uint64_t n = 1;
        for (std::size_t i=0; i<size_wheel_divs; ++i)
            n *= wheel_divs[i];
        return n;
    }
    static constexpr std::uint64_t wheel_size = get_wheel_size();
    static_assert(wheel_size % 2 == 0);
    static constexpr std::uint64_t half_wheel_size = wheel_size/2;
 

    template <class DUMMY = void>
    static constexpr uint64_t get_num_spokes()
    {
        // We should be able to get the count via Euler's totient function,
        // given the knowledge all divisors are prime (we know all are prime
        // because of the static_asserts after defining wheel_divs).
        static_assert(size_wheel_divs > 0);
        static_assert(wheel_divs[0] == 2);
        uint64_t count = 1;
        for (std::size_t i=0; i<size_wheel_divs; ++i)
            count = count * (wheel_divs[i] - 1);

        return count;
    }
    static constexpr uint64_t num_spokes = get_num_spokes();


#ifdef __GNUC__
#  pragma GCC diagnostic push
// gcc issues a pointless warning for the std::conditionals below that a
// comparison is always false, which we disable here
#  pragma GCC diagnostic ignored "-Wtype-limits"
#endif

    // WI_T should be the smallest type that's able to hold any wheel index.
    using WI_T = typename std::conditional<
         (wheel_size <= ut_numeric_limits<std::uint8_t>::max()),
         std::uint8_t,
         typename std::conditional<
             (wheel_size <= ut_numeric_limits<std::uint16_t>::max()),
             std::uint16_t,
             typename std::conditional<
                 (wheel_size <= ut_numeric_limits<std::uint32_t>::max()),
                 std::uint32_t,
                 std::uint64_t
             >::type
         >::type
      >::type;
    static_assert(wheel_size <= std::numeric_limits<WI_T>::max());


    // SI_T should be the smallest type that's able to hold any spoke index.
    using SI_T = typename std::conditional<
         (num_spokes <= ut_numeric_limits<std::uint8_t>::max()),
         std::uint8_t,
         typename std::conditional<
             (num_spokes <= ut_numeric_limits<std::uint16_t>::max()),
             std::uint16_t,
             typename std::conditional<
                 (num_spokes <= ut_numeric_limits<std::uint32_t>::max()),
                 std::uint32_t,
                 std::uint64_t
             >::type
         >::type
      >::type;
    static_assert(num_spokes <= std::numeric_limits<SI_T>::max());

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif


    // the returned array's indices represent odd numbers.
    template <class DUMMY = void>
    static constexpr std::array<SI_T, half_wheel_size>
    get_wheel_reindex()
    {
        static_assert(wheel_divs[0] == 2);
        // Create sieve to test coprimality to the wheel_divs.
        // Indices represent odd numbers.
        // We know all wheel_divs after wheel_divs[0] are odd, because
        // of the static_asserts after defining wheel_divs.

        std::array<SI_T, half_wheel_size> odd_coprime_sieve{};
        for (WI_T i=0; i<half_wheel_size; ++i)
            odd_coprime_sieve[i] = 1;

        // we skip j=0 since implicitly wheel_divs[0] is already handled.
        for (std::size_t j=1; j<size_wheel_divs; ++j) {
            std::uint64_t wd = wheel_divs[j];
            std::uint64_t two_wd = wd + wd;
            odd_coprime_sieve[static_cast<std::size_t>(wd/2)] = 0;
            for (std::uint64_t i=wd*wd; i<wheel_size; i+=two_wd) {
                HPBC_CONSTEXPR_ASSERT(i/2 < half_wheel_size);
                static_assert((wheel_size-1)/2 <=
                              std::numeric_limits<std::size_t>::max());
                odd_coprime_sieve[static_cast<std::size_t>(i/2)] = 0;
            }
        }

        // repurpose the odd_coprime_sieve to be used for reindexing
        SI_T count = 0;
        for (WI_T i=0; i<half_wheel_size; ++i) {
            if (odd_coprime_sieve[i] != 0)
                odd_coprime_sieve[i] = count++;
        }
        HPBC_CONSTEXPR_ASSERT(count == num_spokes);
        return odd_coprime_sieve;
    }
    static constexpr std::array<SI_T, half_wheel_size> wheel_reindex = get_wheel_reindex();


    template <class DUMMY = void>
    static constexpr std::array<WI_T, num_spokes> get_spokes()
    {
        std::array<WI_T, num_spokes> tmp_spokes{};
        // wheel_reindex[0] will be 0, but it's a reindex nonetheless.
        // (keep in mind wheel_reindex's indices represent odd numbers)
        static_assert(wheel_reindex[0] == 0);
        // Thus we set tmp_spokes[0] manually, and start count at 1 to include it.
        tmp_spokes[0] = static_cast<WI_T>(2*0 + 1);
        SI_T count = 1;
        for (WI_T i=1; i<half_wheel_size; ++i) {
            if (wheel_reindex[i] != 0)
                tmp_spokes[count++] = static_cast<WI_T>(2*i + 1);
        }
        HPBC_CONSTEXPR_ASSERT(count == num_spokes);
        return tmp_spokes;
    }
    static constexpr std::array<WI_T, num_spokes> spokes = get_spokes();


// This next function should be correct, but it is disabled for now because
// this class does not use or need it, at the moment.
#if 0
    template <class DUMMY = void>
    static constexpr std::uint32_t get_number_from_table_index(std::uint32_t index)
    {
        std::uint32_t quotient = index/num_spokes;
        std::uint32_t remainder = index - quotient*num_spokes;
        HPBC_CONSTEXPR_ASSERT(remainder < num_spokes);
        std::uint64_t n = static_cast<std::uint64_t>(quotient)*wheel_size + spokes[remainder];

        // although this tests n, it's really a test that index isn't too large.
        HPBC_CONSTEXPR_PRECONDITION(n <= std::numeric_limits<std::uint32_t>::max());

        std::uint32_t n32 = static_cast<std::uint32_t>(n);
        return n32;
    }
#endif

    template <class DUMMY = void>
    static constexpr std::uint32_t get_table_index_from_number(std::uint32_t n)
    {
        // precondition: n is not divisible by any of the wheel_divs
        HPBC_CONSTEXPR_PRECONDITION([&]{
            static_assert(size_wheel_divs > 0);
            bool is_coprime = true;
            for (std::size_t i=0; i<size_wheel_divs; ++i)
                if (n % wheel_divs[i] == 0)
                    is_coprime = false;
            return is_coprime;
        }());

        std::uint32_t quotient = static_cast<std::uint32_t>(n/wheel_size);
        std::uint32_t remainder = static_cast<std::uint32_t>(n - quotient*wheel_size);
        HPBC_CONSTEXPR_ASSERT(remainder < wheel_size);
        static_assert(num_spokes < wheel_size);
        std::uint32_t index = static_cast<std::uint32_t>(
                              quotient*num_spokes + wheel_reindex[remainder/2]);
        return index;
    }


    template <class DUMMY = void>
    static constexpr uint32_t get_num_table_elements()
    {
        uint32_t maxval = std::numeric_limits<uint32_t>::max();
        static_assert(0 < INPUT_BIT_LIMIT && INPUT_BIT_LIMIT <= 32);
        if constexpr (INPUT_BIT_LIMIT < 32)
            maxval = (static_cast<uint32_t>(1) << INPUT_BIT_LIMIT) - 1;
        bool finished = false;
        while (!finished) {
            size_t i=0;
            for (; i<size_wheel_divs; ++i) {
                if (maxval % wheel_divs[i] == 0) {
                    --maxval;
                    break;
                }
            }
            finished = (i >= size_wheel_divs);
        }
        uint32_t index_of_maxval = get_table_index_from_number(maxval);
        HPBC_CONSTEXPR_ASSERT(index_of_maxval < std::numeric_limits<uint32_t>::max());
        return index_of_maxval + 1;
    }
    static constexpr uint32_t num_table_elements = get_num_table_elements();


    // For FAVOR_SMALL_SIZE == true:
    // We will encode factor values with a reduced number of bits by realizing
    // that all factors are prime.  Thus we will encode by using the particular
    // index into PrimesUnder65536::prime array that corresponds to a prime we
    // need for a factor.
    // 13 bits is enough to cover the indices of all primes < 2^16, since there
    // are 6542 of those primes, which is less than 2^13.
    // We use one extra bit to store whether the quotient (the number divided by
    // the factor stored in the table entry) is prime, which is 14 bits total.
    //
    // For FAVOR_SMALL_SIZE == false:
    // When we're not encoding with prime indices, we naively need 2^16 numbers
    // to represent at least one of the factors of a uint32_t input. But we know
    // we will never need to look up any numbers divisible by a wheel_divisor,
    // so if we wanted we could encode by using get_table_index_from_number().
    // We'd need a little under 2^14 numbers in that case.  However we take a
    // much simpler approach and only filter out all even numbers (since
    // wheel_divs[0] == 2), which means we only need 2^15 numbers.  When
    // we add an extra bit to store whether the quotient (see above) is prime,
    // that gets us to 16, and 16 is an efficient table entry size.
    //
    static constexpr unsigned int
                              TABLE_ENTRY_BITLEN = (FAVOR_SMALL_SIZE) ? 14 : 16;


    using TableType = BitpackedUintVector<uint16_t, TABLE_ENTRY_BITLEN>;


    static TableType makePopulatedTable()
    {
        using size_type = typename TableType::size_type;
        // ensure the initializer for table was valid
        static_assert(sizeof(size_type) >= sizeof(num_table_elements));

        TableType tmp_table(num_table_elements);

        // set up a temporary primes_reindex (it's only needed for code of
        // FAVOR_SMALL_SIZE == true).
        std::vector<std::uint16_t> primes_reindex(65536, 0);
        if constexpr (FAVOR_SMALL_SIZE) {
            for (std::uint16_t i=0; i < PrimesUnder65536<>::num_primes; ++i)
                primes_reindex[PrimesUnder65536<>::prime[i]] = i;
            // Note!!  the prime 2 reindexes to 0, and we rely on this below.
            HPBC_ASSERT2(primes_reindex[2] == 0);

            // Since the factorize() call will factor out any 2 via the wheel,
            // we know 2 will never (normally) be an entry that's used in the
            // factor table.  We'll repurpose 2 to mean that there are no factors
            // (i.e. the number is prime).  2 reindexes to the value 0.
            static_assert(wheel_divs[0] == 2);
        }

        for (uint32_t i = 0;; ++i) {
            for (uint32_t j = 0; j < num_spokes; ++j) {
                uint64_t n = static_cast<uint64_t>(i)*wheel_size + spokes[j];
                static_assert(0 < INPUT_BIT_LIMIT && INPUT_BIT_LIMIT <= 32);
                if (n > (static_cast<uint64_t>(1) << INPUT_BIT_LIMIT))
                    goto end;
                uint32_t n32 = static_cast<uint32_t>(n);
                uint32_t index = static_cast<uint32_t>(i*num_spokes + j);
                HPBC_ASSERT2(index == get_table_index_from_number(n32));

                uint32_t encoded;
                if (n32 < 2) {
                    // factoring n32 in these cases would be undefined - just
                    // use 0, which treats n32 as if it is prime.
                    encoded = 0;
                }
                else {
                    unsigned int num_factors;
                    auto array = hurchalla::factorize(n32, num_factors);
                    HPBC_ASSERT2(num_factors > 0);  // factorize() guarantees
                    if (num_factors == 1) {
                        // n32 is prime, which we indicate with 0
                        encoded = 0;
                    } else {
                        HPBC_ASSERT2(num_factors >= 2);
                        uint32_t largest_storable_factor = 0;
                        for (unsigned int k = 0; k < num_factors; ++k) {
                            if (array[k] > largest_storable_factor &&
                                        array[k] < (static_cast<uint32_t>(1) << 16))
                                largest_storable_factor = array[k];
                        }
                        HPBC_ASSERT2(largest_storable_factor != 0);

                        HPBC_ASSERT2(largest_storable_factor <=
                                     std::numeric_limits<uint16_t>::max());

                        // wheel_divs[0] == 2, so any number that survives
                        // the wheel's trial division would have already had all
                        // factors == 2 removed.
                        static_assert(wheel_divs[0] == 2);
                        if constexpr (FAVOR_SMALL_SIZE) {
                            // By the logic above, the factor 2 shouldn't be
                            // needed. (and it would be unacceptable here as a
                            // factor because primes_reindex[2] == 0, and we use
                            // 0 to indicate n32 is prime.)
                            HPBC_ASSERT2(primes_reindex[2] == 0);
                            HPBC_ASSERT2(largest_storable_factor != 2);
                            // Encode the factor in less than 16 bits, by using
                            // the index into PrimesUnder65536::prime.
                            encoded = primes_reindex[largest_storable_factor];
                        } else {
                            HPBC_ASSERT2(largest_storable_factor != 2);
                            // all other factors are prime, and thus odd.
                            HPBC_ASSERT2(largest_storable_factor % 2 == 1);
                            // since all our factors are odd, we can encode
                            // div 2 (and later decode using 2*encoded + 1).
                            encoded = largest_storable_factor / 2;
                            // since 0 and 1 can't be factors, encoded isn't 0.
                            HPBC_ASSERT2(encoded != 0); //we reserve 0 to mean prime
                        }

                        // quotient is n32 / decoded_factor.
                        // We'll use a bit to store whether quotient is prime
                        bool quotient_is_prime = (num_factors == 2);
                        HPBC_ASSERT2(encoded <
                          (static_cast<uint32_t>(1) << (TABLE_ENTRY_BITLEN-1)));
                        encoded = encoded << 1;
                        uint32_t lowbit = (quotient_is_prime) ? 1 : 0;
                        encoded = encoded | lowbit;

                        HPBC_ASSERT2(encoded != 0);
                    }
                }
                // write at index the (encoded) prime factor of n32
                HPBC_ASSERT2(
                    encoded < (static_cast<uint32_t>(1) << TABLE_ENTRY_BITLEN));
                static_assert(sizeof(size_type) >= sizeof(index));
                tmp_table.setAt(static_cast<size_type>(index),
                                static_cast<uint16_t>(encoded));
            }
        }
    end:
        return tmp_table;
    }


    static bool can_open_file(const char* filepath)
    {
        std::ifstream fs(filepath);
        return fs.is_open();
    }


private:
    TableType table;

public:
    ImplFactorByTable32(const ImplFactorByTable32&) = delete;
    ImplFactorByTable32(ImplFactorByTable32&& other) :
                table(std::move(other.table)) {}

    ImplFactorByTable32() : table(makePopulatedTable()) {}

    // can throw from a file open failure, a read failure, or mismatch in values
    // read vs values expected
    ImplFactorByTable32(const char* table_filepath, bool createTableIfCantOpen = true) :
            table((can_open_file(table_filepath) || !createTableIfCantOpen)
                  ? deserialize(table_filepath) : makePopulatedTable())
    {}


    // can throw from a file open failure, write failure
    void writeTableToFile(const char* table_filepath) const
    {
        serialize(table_filepath);
    }


    std::array<std::uint32_t, 32> 
    factorize(std::uint32_t x, unsigned int& num_factors) const
    {
        HPBC_PRECONDITION2(x >= 2);  // 0 and 1 do not have prime factorizations.
        static_assert(0 < INPUT_BIT_LIMIT && INPUT_BIT_LIMIT <= 32);
        HPBC_PRECONDITION2(x < (static_cast<uint64_t>(1) << INPUT_BIT_LIMIT));

        std::array<std::uint32_t, 32> factors{};
        num_factors = 0;
        if (x < 2)
            return factors;

        uint32_t q = x;
        hurchalla::Unroll<size_wheel_divs>::call([&](size_t i){
            while (q % wheel_divs[i] == 0) {
                factors[num_factors++] = wheel_divs[i];
                q = static_cast<uint32_t>(q / wheel_divs[i]);
            }
        });
        if constexpr (extra_divisors.size() != 0) {
            hurchalla::Unroll<extra_divisors.size()>::call([&](size_t i){
                while (q % extra_divisors[i] == 0) {
                    factors[num_factors++] = extra_divisors[i];
                    q = static_cast<uint32_t>(q / extra_divisors[i]);
                }
            });
        }

        HPBC_ASSERT2(q >= 1);
        while(q != 1) {
#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif
            if (HPBC_ASSERT2_MACRO_IS_ACTIVE) {
                for (unsigned int i=0; i<size_wheel_divs; ++i) {
                    HPBC_ASSERT2(q % wheel_divs[i] != 0);
                }
            }
#if defined(_MSC_VER)
#  pragma warning(pop)
#endif
            uint32_t index = get_table_index_from_number(q);

            using size_type = typename TableType::size_type;
            static_assert(sizeof(size_type) >= sizeof(index));
            uint16_t encoded = table.getAt(index);

            // quotient is q / decoded_factor; i.e. the value
            // of q in the next iteration of this loop.
            bool quotient_is_prime = encoded & 1;
            encoded = encoded >> 1;

            uint32_t tmp;
            if constexpr (FAVOR_SMALL_SIZE) {
                HPBC_ASSERT2(encoded < PrimesUnder65536<>::num_primes);
                tmp = PrimesUnder65536<>::prime[encoded];
            } else {
                tmp = static_cast<uint32_t>(2*encoded + 1);
            }

            // encoded == 0 indicates q is prime
            uint32_t qfactor = (encoded == 0) ? q : tmp;

            factors[num_factors++] = qfactor;
            HPBC_ASSERT2(q % qfactor == 0);
            q = q / qfactor;
            HPBC_ASSERT2(q >= 1);

            if (quotient_is_prime) {
                factors[num_factors++] = q;
                return factors;
            }
        }
        return factors;
    }


    struct FileFormatError : public std::runtime_error {
        FileFormatError(const char* what_arg) : std::runtime_error(what_arg) {}
    };


private:

    // can throw from a file open failure, write failure, or vector resize failure
    void serialize(const char* table_filepath) const
    {
        std::ofstream ofs(table_filepath, std::ios::binary);
        if (ofs.fail())
            throw std::ios_base::failure("couldn't open file");
        ofs.exceptions(std::ofstream::failbit | std::ofstream::badbit |
                       std::ofstream::eofbit);

        uint32_t format = table.getFormatID();
        std::array<unsigned char, 4> format_array;
        for (unsigned int i = 0; i<format_array.size(); ++i)
            format_array[i] = static_cast<unsigned char>(format >> 8*i);

        auto count = table.size();
        HPBC_ASSERT2(count <= std::numeric_limits<uint32_t>::max());
        uint32_t count32 = static_cast<uint32_t>(count);
        static_assert(std::numeric_limits<unsigned char>::digits == 8);
        std::array<unsigned char, 4> count_array;
        for (unsigned int i = 0; i<count_array.size(); ++i)
            count_array[i] = static_cast<unsigned char>(count32 >> 8*i);

        std::size_t datasize = table.dataSizeBytes();
        static_assert(wheel_divs[0] == 2);  // should be enough to guarantee the next assert
        HPBC_ASSERT2(datasize <= std::numeric_limits<uint32_t>::max());
        uint32_t datasize32 = static_cast<uint32_t>(datasize);
        std::array<unsigned char, 4> datasize_array;
        for (unsigned int i = 0; i<datasize_array.size(); ++i)
            datasize_array[i] = static_cast<unsigned char>(datasize32 >> 8*i);

        ofs.write(reinterpret_cast<char*>(format_array.data()), format_array.size());
        ofs.write(reinterpret_cast<char*>(count_array.data()), count_array.size());
        ofs.write(reinterpret_cast<char*>(datasize_array.data()), datasize_array.size());

        // Annoyingly std::streamsize is signed and std::size_t is unsigned.
        // So potentially the table could have a data buffer with size
        // that (quite expectedly) fits in size_t, but that is too large to
        // write via std::streamsize.  Use static_asserts to prevent this.
        static_assert(sizeof(typename TableType::size_type) >= sizeof(num_table_elements));
        HPBC_ASSERT2(count32 == num_table_elements);
        constexpr std::size_t expected_datasize =
                             TableType::dataSizeBytes(num_table_elements);
        static_assert(expected_datasize
                       <= std::numeric_limits<std::streamsize>::max());
        HPBC_ASSERT2(datasize == expected_datasize);  // might as well check too

        ofs.write(reinterpret_cast<const char*>(table.data()),
                  static_cast<std::streamsize>(datasize));
        ofs.close();
    }


    // can throw from a file open failure, a read failure, or vector resize failure
    TableType deserialize(const char* table_filepath)
    {
        std::ifstream ifs(table_filepath, std::ios::binary);
        if (ifs.fail())
            throw std::ios_base::failure("couldn't open file");
        ifs.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                       std::ifstream::eofbit);

        std::array<uint8_t, 4> format_array;
        std::array<uint8_t, 4> count_array;
        std::array<uint8_t, 4> datasize_array;
        ifs.read(reinterpret_cast<char*>(format_array.data()), 4);
        ifs.read(reinterpret_cast<char*>(count_array.data()), 4);
        ifs.read(reinterpret_cast<char*>(datasize_array.data()), 4);

        uint32_t format32 = 0;
        uint32_t count32 = 0;
        uint32_t datasize32 = 0;
        for (unsigned int i = 0; i<format_array.size(); ++i)
            format32 += static_cast<decltype(format32)>(format_array[i]) << 8*i;
        for (unsigned int i = 0; i<count_array.size(); ++i)
            count32 += static_cast<decltype(count32)>(count_array[i]) << 8*i;
        for (unsigned int i = 0; i<datasize_array.size(); ++i)
            datasize32 += static_cast<decltype(datasize32)>(datasize_array[i]) << 8*i;

        static_assert(sizeof(std::size_t) >= sizeof(uint32_t));
        std::size_t datasize = static_cast<std::size_t>(datasize32);
        std::unique_ptr<unsigned char[]> data = std::make_unique<unsigned char[]>(datasize);

        // Annoyingly std::streamsize is signed and std::size_t is unsigned.
        // So potentially the table could have a data buffer with size
        // that (quite expectedly) fits in size_t, but that is too large to
        // read via std::streamsize.  Use static_asserts to prevent this.
        static_assert(sizeof(typename TableType::size_type) >=
                      sizeof(num_table_elements));
        constexpr std::size_t expected_datasize =
                             TableType::dataSizeBytes(num_table_elements);
        static_assert(expected_datasize
                       <= std::numeric_limits<std::streamsize>::max());
        if (datasize != expected_datasize || count32 != num_table_elements ||
                    format32 != TableType::getFormatID())
            throw FileFormatError("mismatch in values read vs values expected");

        ifs.read(reinterpret_cast<char*>(data.get()),
                 static_cast<std::streamsize>(datasize));
        ifs.close();

        static_assert(sizeof(typename TableType::size_type) >= sizeof(count32));
        return TableType(std::move(data), datasize, count32);
    }
};


}} // end namespace

#endif
