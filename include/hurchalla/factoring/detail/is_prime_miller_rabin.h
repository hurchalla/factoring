// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_IS_PRIME_MILLER_RABIN_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_MILLER_RABIN_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinProbabilisticBases128.h"

#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_2.h"
#ifdef HURCHALLA_CHOOSE_ALTERNATE_MILLER_RABIN_BASES64_3
# include "hurchalla/factoring/detail/miller_rabin_bases/alternative_tables/MillerRabinBases64_3_alt.h"
#else
# include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_3.h"
#endif
#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_4.h"
#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_5.h"
#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_6.h"
#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases64_7.h"

#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases62_2.h"
#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases62_3.h"

#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases32_1.h"
#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases32_2.h"
#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases32_3.h"

#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases16_1.h"
#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases16_2.h"

#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>
#include <array>
#include <cstdint>
#include <cstddef>

#ifdef __GNUC__
#  pragma GCC diagnostic push
// old versions of gcc and clang give unnecessary warnings about single braced
// initialization lists with std::array (newer versions fixed this).
#  pragma GCC diagnostic ignored "-Wmissing-braces"
#endif

namespace hurchalla { namespace detail {


// Using a large TRIAL_SIZE when possible will often significantly improve
// performance, due to increased opportunity for instruction level parallelism:
// the executing instructions can more efficiently use the CPU's pipelined and
// superscalar execution units.  For example on Intel Haswell, using TRIAL_SIZE
// 3 doubles performance compared to 3 trials using TRIAL_SIZE 1.  You could
// see performance increases up through TRIAL_SIZE 5 and beyond (though it may
// be hard to find an application that needs very large trial sizes).  As
// always, you will need to measure to know the performance impact it will have
// on your machine.  The downside of increasing the TRIAL_SIZE is that it will
// increase the code size and instruction cache usage (to some degree it may
// decrease temporal and spatial code locality in the i-cache).


// FYI: For GCC and TRIAL_SIZE of 3 and 4, this function is faster if we make an
// (enable_if) overload of the function with all the TRIAL_SIZE loops manually
// unrolled.  Manual unroll makes either no difference or is slightly slower for
// clang - I haven't measured icc or msvc.
template <std::size_t TRIAL_SIZE, typename MontType>
HURCHALLA_FLATTEN
//typename std::enable_if<(TRIAL_SIZE != 1), bool>::type
bool mr_trial(const MontType& mf,
              const std::array<typename MontType::T_type, TRIAL_SIZE>& bases,
              typename MontType::T_type d,
              int r)
{
    using T = typename MontType::T_type;
    using V = typename MontType::MontgomeryValue;
    using C = typename MontType::CanonicalValue;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(TRIAL_SIZE > 0, "");

    auto zero = mf.getZeroValue();
    auto unity = mf.getUnityValue();
    auto negativeOne = mf.getNegativeOneValue();

    std::array<V, TRIAL_SIZE> mv_base;
    HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i)
        mv_base[i] = mf.convertIn(bases[i]);

    std::array<bool, TRIAL_SIZE> isProbPrime;
    HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i)
        isProbPrime[i] = (mf.getCanonicalValue(mv_base[i]) == zero);

    std::array<V, TRIAL_SIZE> result = mf.pow(mv_base, d);

    std::array<C, TRIAL_SIZE> canonical;
    HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i)
        canonical[i] = mf.getCanonicalValue(result[i]);

#ifdef HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS
// Usually we don't want to allow even numbers, since they can't be used with
// any true Montgomery Form.  Only MontgomeryStandardMathWrapper will work with
// evens, and it only works because that class wraps standard modular arithmetic
// in the montgomery interface, rather than actually using montgomery math.
// There's also rarely a benefit to allowing even numbers, since they are so
// trivial to test for primality.
    HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i)
        isProbPrime[i] |= (canonical[i] == unity);
    if (r == 0) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i) {
            if (!isProbPrime[i])
                return false;
        }
        return true;
    }
    HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i)
        isProbPrime[i] |= (canonical[i] == negativeOne);
#else
    HPBC_PRECONDITION2(r > 0);
    HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i) {
        isProbPrime[i] = (canonical[i] == unity) | (canonical[i] == negativeOne)
                    | (isProbPrime[i]);
    }
#endif

    for (int j=1; j<r; ++j) {
        HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i)
            result[i] = mf.multiply(result[i], result[i]);
        HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i)
            isProbPrime[i] |= (mf.getCanonicalValue(result[i]) == negativeOne);
    }

    HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t i=0; i<TRIAL_SIZE; ++i) {
        if (!isProbPrime[i])
            return false;
    }
    return true;
}


// This function overload is a very slight win when testing numbers that are all
// composite, but loses by a larger margin when testing primes or a mixture of
// composites and primes.  I measured with gcc and clang (no icc or msvc).
// I don't expect there to normally be a reason to strongly favor composites,
// and so this function is disabled by default via #if 0.  If you enable this
// function, you'll need to also uncomment the enable_if line for the first
// mr_trial() function overload above.
#if 0
template <std::size_t TRIAL_SIZE, typename MontType>
HURCHALLA_FLATTEN
typename std::enable_if<(TRIAL_SIZE == 1), bool>::type
mr_trial(const MontType& mf,
         const std::array<typename MontType::T_type, TRIAL_SIZE>& bases,
         typename MontType::T_type d,
         int r)
{
    static_assert(TRIAL_SIZE == 1, "");
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");

    auto zero = mf.getZeroValue();
    auto unity = mf.getUnityValue();
    auto negativeOne = mf.getNegativeOneValue();

    // note: mf.convertIn(bases[0]) effectively behaves as if it performs
    // bases[0] % num, prior to converting to montgomery form.
    auto mv_base = mf.convertIn(bases[0]);
    if (mf.getCanonicalValue(mv_base) == zero)
        return true;
// **
    auto result = mf.pow(mv_base, d);
// **
    auto canonicalResult = mf.getCanonicalValue(result);

#ifdef HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS
    if (canonicalResult == unity)
        return true;
    if (r == 0)
        return false;
    if (canonicalResult == negativeOne)
        return true;
#else
    HPBC_PRECONDITION2(r > 0);
    if (canonicalResult == unity || canonicalResult == negativeOne)
        return true;
#endif

    for (int j=1; j<r; ++j) {
        result = mf.multiply(result, result);
        if (mf.getCanonicalValue(result) == negativeOne)
            return true;
    }
    return false;
}
#endif



// Miller-Rabin first step:
// write num-1 as 2^r * d by factoring powers of 2 from num-1
template <typename T>
void extract_powers_of_two_from_num_minus_one(T num, T& d, int& r)
{
    HPBC_PRECONDITION2(num >= 2);
    d = static_cast<T>(num - 1);
    HPBC_ASSERT2(d > 0);
    r = 0;
// Usually we don't want to allow even numbers, since they can't be used with
// any true Montgomery Form.  Only MontgomeryStandardMathWrapper will work with
// evens, and it only works because that class wraps standard modular arithmetic
// in the montgomery interface, rather than using actual montgomery math.
// There's also rarely any benefit to allowing even numbers, since they are so
// trivial to test for primality.
#ifdef HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS
    while (d % 2 == 0) {
        ++r;
        d = static_cast<T>(d/2);
    }
#else
    HPBC_ASSERT2(num % 2 == 1);
    HPBC_ASSERT2(d % 2 == 0);
    do {
        ++r;
        d = static_cast<T>(d/2);
    } while (d % 2 == 0);
    HPBC_ASSERT2(r > 0);
#endif
}



template <std::size_t TRIAL_SIZE, std::size_t TOTAL_BASES,
          typename B, typename MontType>
HURCHALLA_FORCE_INLINE
typename std::enable_if<(TOTAL_BASES % TRIAL_SIZE != 0), bool>::type
miller_rabin_trials(
                const MontType& mf,
                const std::array<B,TOTAL_BASES>& bases
                )
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<B>::is_integer, "");
    static_assert(ut_numeric_limits<T>::max() >=
                  ut_numeric_limits<B>::max(), "");
    static_assert(TRIAL_SIZE > 0, "");
    static_assert(TOTAL_BASES > 0, "");

    static constexpr std::size_t REMAINDER = TOTAL_BASES % TRIAL_SIZE;
    static_assert(REMAINDER > 0, "");

    // convertIn() inside mr_trial() effectively performs bases[i] % num, prior
    // to using any base.  While setting the bases % num is essential, it would
    // be redundant work if we did it here (or anywhere in addition to mr_trial)
    T num = mf.getModulus();
    T d;  int r;
    extract_powers_of_two_from_num_minus_one(num, d, r);
    {
        std::array<T, REMAINDER> bases_chunk;
        HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t j=0; j < REMAINDER; ++j)
            bases_chunk[j] = bases[j];
        bool isProbablyPrime = mr_trial(mf, bases_chunk, d, r);
        if (!isProbablyPrime)
            return false;
    }
    for (std::size_t i=REMAINDER; i < TOTAL_BASES; i += TRIAL_SIZE) {
        std::array<T, TRIAL_SIZE> bases_chunk;
        HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t j=0; j < TRIAL_SIZE; ++j)
            bases_chunk[j] = bases[i + j];
        bool isProbablyPrime = mr_trial(mf, bases_chunk, d, r);
        if (!isProbablyPrime)
            return false;
    }
    return true;
}

template <std::size_t TRIAL_SIZE, std::size_t TOTAL_BASES,
          typename B, typename MontType>
HURCHALLA_FORCE_INLINE
typename std::enable_if<(TOTAL_BASES % TRIAL_SIZE == 0), bool>::type
miller_rabin_trials(
                const MontType& mf,
                const std::array<B,TOTAL_BASES>& bases
                )
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(ut_numeric_limits<B>::is_integer, "");
    static_assert(ut_numeric_limits<T>::max() >=
                  ut_numeric_limits<B>::max(), "");
    static_assert(TRIAL_SIZE > 0, "");
    static_assert(TOTAL_BASES > 0, "");
    static_assert(TOTAL_BASES % TRIAL_SIZE == 0, "");

    // convertIn() inside mr_trial() effectively performs bases[i] % num, prior
    // to using any base.  While setting the bases % num is essential, it would
    // be redundant work if we did it here (or anywhere in addition to mr_trial)
    T num = mf.getModulus();
    T d;  int r;
    extract_powers_of_two_from_num_minus_one(num, d, r);
    for (std::size_t i=0; i < TOTAL_BASES; i += TRIAL_SIZE) {
        std::array<T, TRIAL_SIZE> bases_chunk;
        HURCHALLA_REQUEST_UNROLL_LOOP for (std::size_t j=0; j < TRIAL_SIZE; ++j)
            bases_chunk[j] = bases[i + j];
        bool isProbablyPrime = mr_trial(mf, bases_chunk, d, r);
        if (!isProbablyPrime)
            return false;
    }
    return true;
}




// -----------------------------------------------------------------------------
// We will assume that the caller uses an optimally performing MontType,
// considering the particular value of mf's modulus (a small value can allow
// a small MontType).  By placing any work of finding and instantiating an ideal
// MontType onto the caller, we don't have to do any of that work here.
// The is_prime functions will still work properly with a sub-optimal MontType;
// we just don't concern ourselves here with trying to find a better MontType.
// Note that none of the primality testing functions in this file explicitly
// reduce any bases by the tested number: the miller rabin algorithm is
// implemented mostly in mr_trial(), which has a call to mf.convertIn(), and
// convertin() implicitly reduces each base mod num while converting the base
// into MontgomeryForm.  There's no need to ever explicitly reduce a base when
// convertin() effectively does it for us.
// -----------------------------------------------------------------------------


// -----------------
// Miller-rabin primality testing functions with fixed-width integer type
// preconditions on the modulus size.
// -----------------
// Template parameters:
// LOG2_MODULUS_LIMIT can be 128, 64, 32, or 16, and creates a precondition of
// modulus < (1 << LOG2_MODULUS_LIMIT).  Typically the type T satisfies this at
// compile time (e.g. a type T of uint64_t with a LOG2_MODULUS_LIMIT of 64), 
// but if a type T bit-width is larger than LOG2_MODULUS_LIMIT, then the caller
// needs to ensure at run-time that the precondition will be satisfied.
// Primality testing is generally faster with smaller LOG2_MODULUS_LIMIT.
// The TOTAL_BASES valid choices depend upon the choice of LOG2_MODULUS_LIMIT.
// Here are the valid combinations of LOG2_MODULUS_LIMIT and TOTAL_BASES, along
// with the resulting internal hash table size.  (Roughly speaking, using fewer
// TOTAL_BASES is faster, but it trades more memory/cache use by a hash table
// for possible speed gain.)
//   16 bit LOG2_MODULUS_LIMIT, 1 base TOTAL_BASES - 8 byte hash table
//   16 bit, 2 bases - 0 bytes (no hash table used)
//
//   32 bit, 1 base - 512 byte hash table
//   32 bit, 2 bases - 16 bytes
//   32 bit, 3 bases - 0 bytes (no hash table used)
//
//   64 bit, 2 bases - 448 KB hash table
//   64 bit, 3 bases - 25.5 KB
//   64 bit, 4 bases - 2.75 KB
//   64 bit, 5 bases - 320 bytes
//   64 bit, 6 bases - 40 bytes
//   64 bit, 7 bases - 0 bytes (no hash table used)
//
//   128 bit, 128 bases - 0 bytes (no hash table used, probabilistic testing)

// Note: when HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS is defined, this test
// works for both even and odd moduli (though only the MontyWrappedStandardMath
// MontyType supports an even modulus).

// helper constexpr function
template <int LOG2_MODULUS_LIMIT>
constexpr int get_MRM_POW2_LIMIT()
{
    int val=1;
    for (; val<LOG2_MODULUS_LIMIT; val*=2)
        ;
    return val;
}
// Primary template.  LOG2_MODULUS_LIMIT of 128 has a partial specialization
template <typename MontType, int LOG2_MODULUS_LIMIT, std::size_t TRIAL_SIZE,
          std::size_t TOTAL_BASES>
struct MillerRabinMontgomery {
  static bool is_prime(const MontType& mf)
  {
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // We restrict the allowed types of T to only those types for which this
    // function would be an efficient choice:
    static_assert((ut_numeric_limits<T>::digits/2 < LOG2_MODULUS_LIMIT &&
                   LOG2_MODULUS_LIMIT <= ut_numeric_limits<T>::digits) ||
              (LOG2_MODULUS_LIMIT < ut_numeric_limits<T>::digits &&
               ut_numeric_limits<T>::digits <= HURCHALLA_TARGET_BIT_WIDTH), "");
    T num = mf.getModulus();

    constexpr int POW2_LIMIT = get_MRM_POW2_LIMIT<LOG2_MODULUS_LIMIT>();
    static_assert(POW2_LIMIT >= LOG2_MODULUS_LIMIT, "");
    using U = typename sized_uint<POW2_LIMIT>::type;
    HPBC_PRECONDITION2(1 < num &&
                    num <= (static_cast<U>(1) << (LOG2_MODULUS_LIMIT-1)) - 1 +
                           (static_cast<U>(1) << (LOG2_MODULUS_LIMIT-1)));
    const auto bases = MillerRabinBases<LOG2_MODULUS_LIMIT, TOTAL_BASES>::
                                                       get(static_cast<U>(num));
    return miller_rabin_trials<TRIAL_SIZE>(mf, bases);
  }
};

// Partial specialization for 128 bit numbers.
// This version requires input < 2^128, and uses 128 bases (no hash tables).
template <typename MontType, std::size_t TRIAL_SIZE>
struct MillerRabinMontgomery<MontType, 128, TRIAL_SIZE, 128> {
  static bool is_prime(const MontType& mf)
  {
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");

    // Other algorithms for primality testing are likely to be far more
    // suitable for large 128 bit numbers than Miller-Rabin, but we can still
    // use Miller-Rabin to do the job.
    //
    // This particular test is unique amongst all others in this file, because
    // it is a probabilistic test.  By nature, the test can be designed for an
    // arbitrarily small chance of failure, and so this test has been tailored
    // (using a huge number of bases: 128) to have an almost inconceivably small
    // chance of ever returning a wrong result.  For details, see
    // MillerRabinProbabilisticBases128.h.
    // The huge number of bases makes it a very slow test when a number is
    // prime, compared to all the other tests in this file.  [But like all
    // miller-rabin tests, when given a composite number, it is on average very
    // fast since it usually can detect compositeness on the first trial.]

    static_assert(ut_numeric_limits<T>::digits == 128, "");
    // I believe it is *only* a good idea to use this specialization for 128 bit
    // types T.  Any smaller type T would be forced to use far more bases than
    // necessary (given the smaller range of possible moduli), and any larger
    // type would require a slower MontyType than 128 bit T would.  [For what
    // it's worth, I expect this specialization should work fine for non-128
    // bit types.]
    //   If for some reason you need this function specialization but your
    // MontyType T is larger than 128 bit, then construct a different MontyType
    // with 128 bit T, and use it to call this function.  Of course your modulus
    // value must be less than 2^128 in order to do this, regardless of your
    // original T size.  I expect it would be much faster to construct the
    // MontyType object and then call this function, than it would be to use
    // this function with a non-128 bit T MontyType.
    //   If your MontyType T is less than 128 bit, then use one of the
    // MillerRabinMontgomery structs that is smaller than this enclosing struct,
    // since it will provide far better performance.

    T num = mf.getModulus();
    HPBC_PRECONDITION2(1 < num);
    const auto& bases = MillerRabinProbabilisticBases128<>::bases;
    return miller_rabin_trials<TRIAL_SIZE>(mf, bases);
  }
};




// ------------------------
// Functions with special preconditions that limit the allowed range of modulus
// (it's inconvenient but can increase speed):
// ------------------------

// Note: when HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS is defined, these tests
// work for both even and odd moduli (though only the MontyWrappedStandardMath
// MontyType supports an even modulus).

// This function requires modulus < 360018361, and uses 2 bases (no hash table).
template <std::size_t TRIAL_SIZE, typename MontType>
bool is_prime_miller_rabin32_2_360018361(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // Technically any 32 bit or larger type should work, but we allow only the
    // types for which this function would be an efficient choice.
    static_assert(ut_numeric_limits<T>::digits == 32 ||
              (ut_numeric_limits<T>::digits > 32 &&
               ut_numeric_limits<T>::digits <= HURCHALLA_TARGET_BIT_WIDTH), "");
    T num = mf.getModulus();
    HPBC_PRECONDITION2(1 < num && num <= UINT32_C(360018361));
    // Dana Jacobsen, Wojciech Izykowski, and Marcin Panasiuk discovered these
    // bases; see https://miller-rabin.appspot.com/
    // I verified they are correct for all num < 360018361.
    const std::array<std::uint32_t, 2> bases =
                                    { UINT32_C(1143370), UINT32_C(2350307676) };
    return miller_rabin_trials<TRIAL_SIZE>(mf, bases);
}

// This function requires modulus < 1050535501, and uses 2 bases (no hash table)
template <std::size_t TRIAL_SIZE, typename MontType>
bool is_prime_miller_rabin64_2_1050535501(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // Technically any 64 bit or larger type should work, but we allow only the
    // types for which this function would be an efficient choice.
    static_assert(ut_numeric_limits<T>::digits == 64 ||
              (ut_numeric_limits<T>::digits > 64 &&
               ut_numeric_limits<T>::digits <= HURCHALLA_TARGET_BIT_WIDTH), "");
    T num = mf.getModulus();
    HPBC_PRECONDITION2(1 < num && num <= UINT32_C(1050535501));
    // Wojciech Izykowski and Marcin Panasiuk discovered these bases; see
    // https://miller-rabin.appspot.com/
    // I verified they are correct for all num < 1050535501.
    const std::array<std::uint64_t, 2> bases =
                         { UINT64_C(336781006125), UINT64_C(9639812373923155) };
    return miller_rabin_trials<TRIAL_SIZE>(mf, bases);
}

// This function requires modulus<350269456337, and uses 3 bases (no hash table)
template <std::size_t TRIAL_SIZE, typename MontType>
bool is_prime_miller_rabin64_3_350269456337(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // Technically any 64 bit or larger type should work, but we allow only the
    // types for which this function would be an efficient choice.
    static_assert(ut_numeric_limits<T>::digits == 64 ||
              (ut_numeric_limits<T>::digits > 64 &&
               ut_numeric_limits<T>::digits <= HURCHALLA_TARGET_BIT_WIDTH), "");
    T num = mf.getModulus();
    HPBC_PRECONDITION2(1 < num && num <= UINT64_C(350269456337));
    // Steve Worley discovered these bases; see https://miller-rabin.appspot.com
    // I verified they are correct for all num < 350269456337.
    const std::array<std::uint64_t, 3> bases = { UINT64_C(4230279247111683200),
               UINT64_C(14694767155120705706), UINT64_C(16641139526367750375) };
    return miller_rabin_trials<TRIAL_SIZE>(mf, bases);
}

// This function requires modulus < 55245642489451; uses 4 bases (no hash table)
template <std::size_t TRIAL_SIZE, typename MontType>
bool is_prime_miller_rabin64_4_55245642489451(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // Technically any 64 bit or larger type should work, but we allow only the
    // types for which this function would be an efficient choice.
    static_assert(ut_numeric_limits<T>::digits == 64 ||
              (ut_numeric_limits<T>::digits > 64 &&
               ut_numeric_limits<T>::digits <= HURCHALLA_TARGET_BIT_WIDTH), "");
    T num = mf.getModulus();
    HPBC_PRECONDITION2(1 < num && num <= UINT64_C(55245642489451));
    // I verified these bases are correct for all num < 55245642489451, assuming
    // Feitsma's database of Fermat pseudoprimes base 2 is valid.  Steve Worley
    // discovered the bases - see https://miller-rabin.appspot.com/
    const std::array<std::uint64_t, 4> bases =
              { UINT64_C(2), UINT64_C(141889084524735),
                UINT64_C(1199124725622454117), UINT64_C(11096072698276303650) };
    return miller_rabin_trials<TRIAL_SIZE>(mf, bases);
}

// This function requires modulus < 7999252175582851; has 5 bases(no hash table)
template <std::size_t TRIAL_SIZE, typename MontType>
bool is_prime_miller_rabin64_5_7999252175582851(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // Technically any 64 bit or larger type should work, but we allow only the
    // types for which this function would be an efficient choice.
    static_assert(ut_numeric_limits<T>::digits == 64 ||
              (ut_numeric_limits<T>::digits > 64 &&
               ut_numeric_limits<T>::digits <= HURCHALLA_TARGET_BIT_WIDTH), "");
    T num = mf.getModulus();
    HPBC_PRECONDITION2(1 < num && num <= UINT64_C(7999252175582851));
    // I verified these bases are correct for all num < 7999252175582851, if
    // Feitsma's database of Fermat pseudoprimes base 2 is valid.  Steve Worley
    // discovered the bases - see https://miller-rabin.appspot.com/
    const std::array<std::uint64_t, 5> bases =
           { UINT64_C(2), UINT64_C(4130806001517), UINT64_C(149795463772692060),
             UINT64_C(186635894390467037), UINT64_C(3967304179347715805) };
    return miller_rabin_trials<TRIAL_SIZE>(mf, bases);
}

// This function requires modulus < 585226005592931977; uses 6 bases (no hash)
template <std::size_t TRIAL_SIZE, typename MontType>
bool is_prime_miller_rabin64_6_585226005592931977(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // Technically any 64 bit or larger type should work, but we allow only the
    // types for which this function would be an efficient choice.
    static_assert(ut_numeric_limits<T>::digits == 64 ||
              (ut_numeric_limits<T>::digits > 64 &&
               ut_numeric_limits<T>::digits <= HURCHALLA_TARGET_BIT_WIDTH), "");
    T num = mf.getModulus();
    HPBC_PRECONDITION2(1 < num && num <= UINT64_C(585226005592931977));
    // I verified these bases are correct for all num < 585226005592931977, if
    // Feitsma's database of Fermat pseudoprimes base 2 is valid.  Steve Worley
    // discovered the bases - see https://miller-rabin.appspot.com/
    const std::array<std::uint64_t, 6> bases =
                { UINT64_C(2), UINT64_C(123635709730000),
                  UINT64_C(9233062284813009), UINT64_C(43835965440333360),
                  UINT64_C(761179012939631437), UINT64_C(1263739024124850375) };
    return miller_rabin_trials<TRIAL_SIZE>(mf, bases);
}




// ----------------------------
// The default functions for miller-rabin primality testing:
// ----------------------------

// Implementation Notes: there are three general principles guiding the choices
// in the default functions below.
// First, we avoid optimizing speed via choosing to use a small number of bases,
// whenever this would require a large hash table to be used.  We can't assume
// large memory and cache usage would be justified for an average caller.
// Nevertheless, we do use this optimization for the defaults if the associated
// hash table is tiny (<= 320 bytes).
// Second, in order to minimize the compiled machine code size, we try to avoid
// causing additional function template instantiations when it's reasonable.
// For example, the 64bit is_prime_miller_rabin template function below has an
// optimization to call is_prime_miller_rabin64_3_350269456337 which fits nicely
// with this guideline because the called function needs instantiations of
// mr_trial() with TRIAL_SIZE 1 and 2 - those exact same instantiations are
// also needed by the template function's MillerRabinMontgomery::is_prime call.
// Thus it causes very little increase in machine code size.
// Third, we normally choose an odd number of TOTAL_BASES and a TRIAL_SIZE of
// 2, which ensures that the first miller-rabin trial has a trial size of 1 and
// all the rest have a trial size of 2.  The first trial will almost always be
// able to detect a composite number, so we maximize its speed at doing this by
// having it check only one base.  If it does not detect a number to be
// composite, it is usually because the number is a prime, and so we assume it's
// unlikely that any of the further trials would detect composite-ness (which
// would let the test end early without running all trials).  Instead we assume
// all the remaining trials will probably run, and so for all remaining trials
// we check more than one base per trial, via setting TRIAL_SIZE to 2.  This
// increases instruction level parallelism in the code, which a pipelined and/or
// superscalar CPU will use to execute more instructions per cycle.  It would
// not be unusual for the functions that use a TRIAL_SIZE 2 to be 50% or more
// faster than they would be with a TRIAL_SIZE 1, so long as the total number of
// bases they process remains the same.  Since we are assuming that all bases
// will need to be processed if the first trial doesn't prove composite-ness, a
// TRIAL_SIZE of 2 provides a significant performance increase in the cases
// where the tested number is prime.  The downside is that functions that use a
// TRIAL_SIZE 2 compile to approximately double the size (in machine code) that
// they would with a TRIAL_SIZE 1.  The involved functions are not huge, so it
// can be justified (probably) for the defaults, but nevertheless it increases
// pressure on the instruction cache.  Increased use of i-cache always has the
// potential for a net negative impact on program performance.  This is why we
// generally avoid a TRIAL_SIZE of 3 or 4 in the defaults, despite the fact that
// when testing primes, we might perhaps expect a doubling in performance
// compared to TRIAL_SIZE 1, due to better instruction level parallelism.  The
// downside of a 3-4x increase in code size (for the involved functions) seems
// too great a cost for defaults.
//
// Just as a FYI reference, on Intel Haswell, performing a single miller-rabin
// trial with TRIAL_SIZE 2 (processing 2 bases) takes roughly 1.2x longer than a
// single trial with TRIAL_SIZE 1 (processing 1 base).  Performing a TRIAL_SIZE
// 3 trial takes about 1.5x longer than a TRIAL_SIZE 1 trial.  A TRIAL_SIZE 4
// trial takes about 1.85x longer than a TRIAL_SIZE 1 trial.  So for Intel
// Haswell, when processing the same number of total bases, TRIAL_SIZE 3 is
// about double the speed of TRIAL_SIZE 1, and TRIAL_SIZE 4 is faster yet.
// (This mostly ignores the more difficult to measure negative effects from
// increased code size and increased instruction cache use, though.)

template <typename MontType>
typename std::enable_if<
      (ut_numeric_limits<typename MontType::T_type>::digits == 128), bool>::type
is_prime_miller_rabin(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // 128 bit miller-rabin with 128 bases is going to be slow no matter what,
    // but a trial size of 4 will usually significantly improve performance
    // over trial size 1, due to more efficient use of the CPU's pipelined
    // and/or superscalar execution units.
    // We typically avoid a TRIAL_SIZE > 2 due to the machine code size increase
    // it causes, but this function processes so many bases that any negative
    // effect on instruction cache should be more than made up for by the
    // repeated speed gain from processing more bases per trial.
    constexpr std::size_t TOTAL_BASES = 128;
    constexpr std::size_t TRIAL_SIZE = 4;
    return MillerRabinMontgomery
                         <MontType, 128, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}

template <typename MontType>
typename std::enable_if<
       (ut_numeric_limits<typename MontType::T_type>::digits == 64), bool>::type
is_prime_miller_rabin(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    if (mf.getModulus() < UINT64_C(350269456337)) {
        // We can use a 3 base test when the modulus is small for a uint64_t.
        // It's faster than the 5 base test below, and needs no extra static
        // memory (it's not hashed) and has essentially no effect on code size.
        constexpr std::size_t TRIAL_SIZE = 2;
        return is_prime_miller_rabin64_3_350269456337<TRIAL_SIZE>(mf);
    }
    // 5 base (hashed) miller-rabin with a trial size 2 should be a good 64 bit
    // default.  It will start with a trial size of 1 (because 5%2 == 1), which
    // achieves max speed for composite numbers unless they are pseudoprimes to
    // the first base.  If testing the first base doesn't prove a number is
    // composite, then the tests will switch to trial size 2 for the remaining 4
    // bases, which should be a significant speed-up over trial size 1 (I
    // measure 1.67x speed-up on Intel Haswell) due to more efficient use of CPU
    // pipelined/superscalar execution units.
    // The obvious 64 bit alternative default is 7 base non-hashed miller-rabin.
    // We choose to use hashed 5 base because it is faster for primes to be able
    // to avoid testing 2 extra bases, and the hashed version's downside of
    // needing a 320 byte hash table in static memory (and accessing it) should
    // usually be a worthwhile trade-off for the speed improvement.
    // You can switch the default to 7 base miller-rabin by pre-defining the
    // macro HURCHALLA_DEFAULT_TO_UNHASHED_MILLER_RABIN
#ifdef HURCHALLA_DEFAULT_TO_UNHASHED_MILLER_RABIN
    constexpr std::size_t TOTAL_BASES = 7;
#else
    constexpr std::size_t TOTAL_BASES = 5;
#endif
    constexpr std::size_t TRIAL_SIZE = 2;
    return MillerRabinMontgomery
                          <MontType, 64, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}

template <typename MontType>
typename std::enable_if<
       (ut_numeric_limits<typename MontType::T_type>::digits == 32), bool>::type
is_prime_miller_rabin(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // 2 base (hashed) miller-rabin with a trial size 1 should be a good 32 bit
    // default, since it's faster than 3 base (non-hashed) miller-rabin yet
    // keeps essentially the same code size and static memory usage.
    // Some compilers may use 16 bytes of static memory for the hash table
    // rather than loading the values on the fly (which uses no static memory).
    // If this occurs and it is unacceptable for you, you can switch the default
    // to 3 base miller-rabin by pre-defining the macro
    // HURCHALLA_DEFAULT_TO_UNHASHED_MILLER_RABIN
#ifdef HURCHALLA_DEFAULT_TO_UNHASHED_MILLER_RABIN
    constexpr std::size_t TOTAL_BASES = 3;
#else
    constexpr std::size_t TOTAL_BASES = 2;
#endif
    constexpr std::size_t TRIAL_SIZE = 1;
    return MillerRabinMontgomery
                          <MontType, 32, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}

// It's questionable whether using miller-rabin is a good idea for primality
// testing uint16_t values.  For any intensive repeated primality testing,
// initializing a bit vector using sieve of eratosthenes will provide a small
// (8KB) lookup table that should be faster at determining primality.  For a one
// time primality test you could just trial divide by all primes < 256, and
// even if it's not quite as fast as this function, the computational load will
// be tiny.  On the plus side, this function may be the most convenient option
// if you are already using the miller-rabin tests in this file.
template <typename MontType>
typename std::enable_if<
       (ut_numeric_limits<typename MontType::T_type>::digits == 16), bool>::type
is_prime_miller_rabin(const MontType& mf)
{
    using T = typename MontType::T_type;
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!ut_numeric_limits<T>::is_signed, "");
    // 1 base (hashed) miller-rabin with a trial size 1 should be a good 16 bit
    // default, since it's faster than 2 base (non-hashed) miller-rabin yet
    // keeps essentially the same code size and static memory usage.
    // Some compilers may use 8 bytes of static memory for the hash table rather
    // than loading the values on the fly (which uses no static memory).
    // If this occurs and it is unacceptable for you, you can switch the default
    // to 2 base miller-rabin by pre-defining the macro
    // HURCHALLA_DEFAULT_TO_UNHASHED_MILLER_RABIN
#ifdef HURCHALLA_DEFAULT_TO_UNHASHED_MILLER_RABIN
    constexpr std::size_t TOTAL_BASES = 2;
#else
    constexpr std::size_t TOTAL_BASES = 1;
#endif
    constexpr std::size_t TRIAL_SIZE = 1;
    return MillerRabinMontgomery
                          <MontType, 16, TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}


// Integral argument versions:

// I'm assuming a caller would usually have already used trial division to find
// any factors<256 for a number x.  Thus a caller would (usually) have no need
// to determine primality for any x < 65536.  By this reasoning there is
// probably little benefit in an is_prime_miller_rabin_integral() function for
// uint16_t.  We start instead by enable_if'ing for types T <= 32 bit, and
// casting T to uint32_t to minimize further function instantations.

template <typename T>
typename std::enable_if<(ut_numeric_limits<T>::digits <= 32), bool>::type
is_prime_miller_rabin_integral(T x)
{
    HPBC_PRECONDITION2(x % 2 == 1);
    static_assert(ut_numeric_limits<T>::is_integer, "");
#ifndef HURCHALLA_TARGET_BIT_WIDTH
#   error "HURCHALLA_TARGET_BIT_WIDTH must be defined"
#endif
#if HURCHALLA_TARGET_BIT_WIDTH >= 64
    return is_prime_miller_rabin(
          MontgomeryQuarter<std::uint64_t>(static_cast<std::uint64_t>(x)));
#else
    using U = std::uint32_t;
    constexpr U Rdiv4 = static_cast<U>(1) << 30;
    return (x < Rdiv4) ?
            is_prime_miller_rabin(MontgomeryQuarter<U>(static_cast<U>(x))) :
            is_prime_miller_rabin(MontgomeryFull<U>(static_cast<U>(x)));
#endif
}

template <typename T>
typename std::enable_if<(32 < ut_numeric_limits<T>::digits) &&
                        (ut_numeric_limits<T>::digits <= 64),
                    bool>::type
is_prime_miller_rabin_integral(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION2(x % 2 == 1);
    if (x > ut_numeric_limits<std::uint32_t>::max()) {
        using U = std::uint64_t;
        constexpr U Rdiv4 = static_cast<U>(1) << 62;
        return (x < Rdiv4) ?
                is_prime_miller_rabin(MontgomeryQuarter<U>(static_cast<U>(x))) :
                is_prime_miller_rabin(MontgomeryFull<U>(static_cast<U>(x)));
    }
    else
        return is_prime_miller_rabin_integral(static_cast<std::uint32_t>(x));
}

template <typename T>
typename std::enable_if<(64 < ut_numeric_limits<T>::digits), bool>::type
is_prime_miller_rabin_integral(T x)
{
    static_assert(ut_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION2(x % 2 == 1);
    if (x > ut_numeric_limits<std::uint64_t>::max()) {
        using U = typename extensible_make_unsigned<T>::type;
        constexpr U Rdiv4 =
                        static_cast<U>(1) << (ut_numeric_limits<U>::digits - 2);
        return (x < Rdiv4) ?
                is_prime_miller_rabin(MontgomeryQuarter<U>(static_cast<U>(x))) :
                is_prime_miller_rabin(MontgomeryFull<U>(static_cast<U>(x)));
    }
    else
        return is_prime_miller_rabin_integral(static_cast<std::uint64_t>(x));
}


}}  // end namespace

#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif


#endif
