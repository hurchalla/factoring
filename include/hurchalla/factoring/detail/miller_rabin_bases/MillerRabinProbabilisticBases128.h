// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_PROBABILISTIC_BASES128_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_PROBABILISTIC_BASES128_H_INCLUDED


#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// These 128 bases let you conduct a probabilistic primality test for any 128
// bit unsigned int number, via miller-rabin primality testing.  A test using
// these bases should work equally well for both even and odd numbers.  Note
// however that montgomery arithmetic (which is one way to implement the miller-
// rabin test) always requires odds.

// Probability Details:
// Using Miller-Rabin for large 128 bit numbers, there is no known set of bases
// to use for a deterministic test.  Instead we depend upon providing a huge
// number of bases, so that there will be an incredibly low probability that
// miller-rabin testing will ever declare that a composite number is prime after
// completing all trials.  According to
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Miller%E2%80%93Rabin_test
// if k is the number of bases, then the probability is 4^(-k) that all trials
// will declare that a particular composite number is prime (it will never
// declare a prime is composite).  So with 64 bases, if we tested every number
// from 1 to 2^128, we might expect very roughly about 1 false declaration of
// prime.  With an additional 64 bases (for a total of 128), the probability
// that this single expected false declaration would still remain would be very
// roughly 1/2^128, or 1 chance in 340282366920938463463374607431768211456.  See
// the Wikipedia link's discussion of the conditional probability and Bayes'
// theorem for some info on why I wouldn't be surprised if the true probability
// is 100x more likely.  Yet even assuming it's 100000 times more likely, it
// would have odds of 1 in a decillion (roughly one over the number of grains of
// sand on earth, squared) that any number exists between 1 and 2^128 for which
// it would deliver a wrong answer.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY=void>
struct MillerRabinProbabilisticBases128 {
    static_assert(std::is_same<DUMMY, void>::value, "");
    // We use all primes for the bases, simply following the example of
    // http://oeis.org/A014233 (I don't know why they use all primes).
    static constexpr std::array<std::uint16_t, 128> bases = {
                                             2, 3, 5, 7, 11, 13, 17, 19, 23,
        29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
        103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
        179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251,
        257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337,
        347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
        431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
        509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
        607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683,
        691, 701, 709, 719 };
};
template <typename DUMMY> constexpr
std::array<std::uint16_t, 128> MillerRabinProbabilisticBases128<DUMMY>::bases;


}}  // end namespace

#endif
