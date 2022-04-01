// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


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
#include <fstream>
#include <cassert>



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



int read_psps(std::string filename, std::vector<std::uint64_t>& psps)
{
    std::ifstream in_file(filename);
    if(!in_file)
    {  
        std::cout << "Error: file open failed for " << filename << "\n";
        return 1;
    }

    uint64_t psp;
    while (in_file >> psp) {
        psps.push_back(psp);
    }
    if (in_file.bad()) {
        std::cout << "Error: I/O error while reading input file\n";
        return 4;
    }
    else if (in_file.eof()) {
        // file read was successful
        return 0;
    }
    else if (in_file.fail()) {
        std::cout << "Error: non-integer data found in input file\n";
        return 5;
    }
    assert(false);  // we should never reach here.
    return 6;
}



int main(int argc, char* argv[])
{
    // arguments: program name, base2 pseudoprimes file
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << "  SOURCE_FILE_PSEUDOPRIMES_BASE2_2POW64\n";
        return 1;
    }

    using namespace hurchalla::detail;

    using IPMRS = is_prime_miller_rabin_special;

    // we will (hopefully) verify that no pseudoprimes get identified as prime.

    std::cout << "started\n";

    constexpr uint64_t LIMIT_U32 = static_cast<uint64_t>(1) << 32;
    auto primevec = init_primes(LIMIT_U32);

    std::vector<std::uint64_t> psps;
    int status = read_psps(argv[1], psps);
    if (status != 0)
        return 1;

    std::cout << "primes and psps are initialized\n";


    // ----Test uint64_t miller-rabin, modulus limit < 273919523041, 3 base----

    // This may take a half day to run (I modifed it for multithreading when I
    // ran it, but it's single threaded here for simplicity).
    std::cout << "Very long (~half day) test started...\n";

    for (std::uint64_t i = 3; i < UINT64_C(273919523041); i += 2) {
        using MF = hurchalla::MontgomeryQuarter<std::uint64_t>;
        MF mf(i);
        // By theory, miller-rabin is always correct when it claims a number is
        // composite, so to save time we'll only verify the test result when it
        // claims a number is prime (i.e. we will verify that there are no
        // pseudoprimes).

// we use #if 0 to test only one TRIAL_SIZE, in order to save time on this very
// long running test.
#if 0
        if (IPMRS::case_273919523041_64_3<3>(mf) == true) {
            if (MillerRabinMontgomery<MF, 64, 3, 3>::is_prime(mf) == false) {
                std::cout << "failure 13b on " << i << "\n";
                return 13;
            }
        }
        if (IPMRS::case_273919523041_64_3<1>(mf) == true) {
            if (MillerRabinMontgomery<MF, 64, 3, 3>::is_prime(mf) == false) {
                std::cout << "failure 15b on " << i << "\n";
                return 15;
            }
        }
#endif
        if (IPMRS::case_273919523041_64_3<2>(mf) == true) {
            if (MillerRabinMontgomery<MF, 64, 3, 3>::is_prime(mf) == false) {
                std::cout << "failure 14b on " << i << "\n";
                return 14;
            }
        }
    }
    std::cout << "case_273919523041_64_3 for odds is ok\n";
#if 0
// this would take a very long time, and we know from theory that this function
// will always show that any even number > 2 is composite, because the function
// is implemented with at least one even base.
      // verify evens other than 2 are not id'd as prime
    for (std::uint64_t i = 4; i < UINT64_C(273919523041); i += 2) {
        // MontgomeryStandardMathWrapper is the only MontgomeryForm
        // that allows an even modulus.  It can do this since it doesn't use
        // any montgomery math- it just wraps standard math in a monty interface
        hurchalla::MontgomeryStandardMathWrapper<std::uint64_t> mf(i);
        if (IPMRS::case_273919523041_64_3<3>(mf) != false) {
            std::cout << "failure 16b on " << i << "\n";
            return 16;
        }
        if (IPMRS::case_273919523041_64_3<2>(mf) != false) {
            std::cout << "failure 17b on " << i << "\n";
            return 17;
        }
        if (IPMRS::case_273919523041_64_3<1>(mf) != false) {
            std::cout << "failure 18b on " << i << "\n";
            return 18;
        }
    }
    std::cout << "case_273919523041_64_3 for evens ok\n";
#endif
    {  // verify 2 gets id'd as prime
        std::uint64_t i = 2;
        hurchalla::MontgomeryStandardMathWrapper<std::uint64_t> mf(i);
        if (IPMRS::case_273919523041_64_3<3>(mf) != true) {
            std::cout << "failure 19b on " << i << "\n";
            return 19;
        }
        if (IPMRS::case_273919523041_64_3<2>(mf) != true) {
            std::cout << "failure 20b on " << i << "\n";
            return 20;
        }
        if (IPMRS::case_273919523041_64_3<1>(mf) != true) {
            std::cout << "failure 21b on " << i << "\n";
            return 21;
        }
    }
    std::cout << "case_273919523041_64_3 for 2 is ok\n";


    // Note: we don't verify  case_3317044064679887385961981_128_13().
    // It covers far too large a range of numbers for us to feasibly check.  We
    // just assume that the finders of the bases wrote a correct proof that
    // the bases are valid for checking primality of any number less than
    // 3317044064679887385961981.


    std::cout << "Complete - all tests succeeded.\n";
    return 0;
}
