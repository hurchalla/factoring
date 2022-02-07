// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla


// normally we wouldn't define this, since there's almost never a reason to
// test any even number for primality, but here we want to verify that the
// is_prime functions correctly identify even numbers (other than 2) as
// composite.
#define HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS 1


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <vector>
#include <cstddef>


// for easier testing, we'll wrap the is_prime_miller_rabin internals
template <int LOG2_MODULUS_LIMIT, std::size_t TRIAL_SIZE,
          std::size_t TOTAL_BASES, typename MontType>
bool is_prime_mr(const MontType& mf)
{
    namespace hc = hurchalla::detail;
    return hc::MillerRabinMontgomery<MontType, LOG2_MODULUS_LIMIT,
                                TRIAL_SIZE, TOTAL_BASES>::is_prime(mf);
}


// Initialize a bit vector that marks all prime indexes as true up to size,
// and marks all non-prime indexes as false.
std::vector<bool> init_primes(std::uint64_t size)
{
    using std::uint64_t;
    std::vector<bool> primes(size, true);
    primes[0] = false;
    primes[1] = false;
    const auto pb = primes.begin();
    const auto pe = primes.end();
    for (uint64_t j=2*2; j<size; j+=2)
        primes[j] = false;
    for (uint64_t i=3; i*i<size; i = std::find(pb+(i+1),pe,true) - pb)
        for (uint64_t j=i*i; j<size; j+=(2*i))
            primes[j] = false;
    return primes;
}


int main()
{
    namespace hc = hurchalla::detail;

    std::cout << "started\n";

    constexpr uint64_t LIMIT_U16 = static_cast<uint64_t>(1) << 16;
    auto primevec = init_primes(LIMIT_U16);

    std::cout << "primes are initialized\n";


// this bunch all succeeded.
    // -----Test 2 base uint16_t miller-rabin-----
    for (std::uint64_t i = 3; i < LIMIT_U16; i += 2) {
        hurchalla::MontgomeryForm<std::uint16_t>
                                              mf(static_cast<std::uint16_t>(i));
        if (is_prime_mr<16, 2, 2>(mf) != primevec[i]) {
            std::cout << "failure 1 on " << i << "\n";
            return 1;
        }
        if (is_prime_mr<16, 1, 2>(mf) != primevec[i]) {
            std::cout << "failure 2 on " << i << "\n";
            return 1;
        }
    }
    std::cout << "2 base MR for odds passed verification\n";
      // verify evens other than 2 are not id'd as prime
    for (std::uint64_t i = 4; i < LIMIT_U16; i += 2) {
        // MontgomeryStandardMathWrapper is the only MontgomeryForm
        // that allows an even modulus.  It can do this since it doesn't use
        // any montgomery math- it just wraps standard math in a monty interface
        hurchalla::MontgomeryStandardMathWrapper<std::uint16_t>
                                              mf(static_cast<std::uint16_t>(i));
        if (is_prime_mr<16, 2, 2>(mf)) {
            std::cout << "failure 3 on " << i << "\n";
            return 1;
        }
        if (is_prime_mr<16, 1, 2>(mf)) {
            std::cout << "failure 4 on " << i << "\n";
            return 1;
        }
    }
    std::cout << "2 base MR for evens passed verification\n";
    {  // verify 2 gets id'd as prime
        std::uint64_t i = 2;
        hurchalla::MontgomeryStandardMathWrapper<std::uint16_t>
                                              mf(static_cast<std::uint16_t>(i));
        if (!is_prime_mr<16, 2, 2>(mf)) {
            std::cout << "failure 5 on " << i << "\n";
            return 1;
        }
        if (!is_prime_mr<16, 1, 2>(mf)) {
            std::cout << "failure 6 on " << i << "\n";
            return 1;
        }
    }
    std::cout << "2 base MR for the number 2 passed verification\n";


// this bunch all succeeded.
    // -----Test 1 base uint16_t miller-rabin-----
    for (std::uint64_t i = 3; i < LIMIT_U16; i += 2) {
        hurchalla::MontgomeryForm<std::uint16_t>
                                              mf(static_cast<std::uint16_t>(i));
        if (is_prime_mr<16, 1, 1>(mf) != primevec[i]) {
            std::cout << "failure 7 on " << i << "\n";
            return 1;
        }
    }
    std::cout << "1 base MR for odds passed verification\n";
      // verify evens other than 2 are not id'd as prime
    for (std::uint64_t i = 4; i < LIMIT_U16; i += 2) {
        // MontgomeryStandardMathWrapper is the only MontgomeryForm
        // that allows an even modulus.  It can do this since it doesn't use
        // any montgomery math- it just wraps standard math in a monty interface
        hurchalla::MontgomeryStandardMathWrapper<std::uint16_t>
                                              mf(static_cast<std::uint16_t>(i));
        if (is_prime_mr<16, 1, 1>(mf)) {
            std::cout << "failure 8 on " << i << "\n";
            return 1;
        }
    }
    std::cout << "1 base MR for evens passed verification\n";
    {  // verify 2 gets id'd as prime
        std::uint64_t i = 2;
        hurchalla::MontgomeryStandardMathWrapper<std::uint16_t>
                                              mf(static_cast<std::uint16_t>(i));
        if (!is_prime_mr<16, 1, 1>(mf)) {
            std::cout << "failure 9 on " << i << "\n";
            return 1;
        }
    }
    std::cout << "1 base MR for the number 2 passed verification\n";


    std::cout << "Complete - all tests succeeded.\n";
    return 0;
}
