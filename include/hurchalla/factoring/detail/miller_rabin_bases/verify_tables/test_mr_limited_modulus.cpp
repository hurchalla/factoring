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

// this bunch all succeeded.
    // -----Test uint32_t miller-rabin, modulus limit < 360018361, 2 base-----
    for (std::uint32_t i = 3; i < 360018361; i += 2) {
        hurchalla::MontgomeryQuarter<std::uint32_t> mf(i);
        if (IPMRS::case_360018361_32_2<2>(mf) != primevec[i]) {
            std::cout << "failure 1 on " << i << "\n";
            return 1;
        }
        if (IPMRS::case_360018361_32_2<1>(mf) != primevec[i]) {
            std::cout << "failure 2 on " << i << "\n";
            return 2;
        }
    }
    std::cout << "case_360018361_32_2 for odds is ok\n";
      // verify evens other than 2 are not id'd as prime
    for (std::uint32_t i = 4; i < 360018361; i += 2) {
        // MontgomeryStandardMathWrapper is the only MontgomeryForm
        // that allows an even modulus.  It can do this since it doesn't use
        // any montgomery math- it just wraps standard math in a monty interface
        hurchalla::MontgomeryStandardMathWrapper<std::uint32_t> mf(i);
        if (IPMRS::case_360018361_32_2<2>(mf) != false) {
            std::cout << "failure 3 on " << i << "\n";
            return 3;
        }
        if (IPMRS::case_360018361_32_2<1>(mf) != false) {
            std::cout << "failure 4 on " << i << "\n";
            return 4;
        }
    }
    std::cout << "case_360018361_32_2 for evens is ok\n";
    {  // verify 2 gets id'd as prime
        std::uint32_t i = 2;
        hurchalla::MontgomeryStandardMathWrapper<std::uint32_t> mf(i);
        if (IPMRS::case_360018361_32_2<2>(mf) != true) {
            std::cout << "failure 5 on " << i << "\n";
            return 5;
        }
        if (IPMRS::case_360018361_32_2<1>(mf) != true) {
            std::cout << "failure 6 on " << i << "\n";
            return 6;
        }
    }
    std::cout << "case_360018361_32_2 for 2 is ok\n";


// this bunch all succeeded.
    // -----Test uint64_t miller-rabin, modulus limit < 1050535501, 2 base-----
    for (std::uint64_t i = 3; i < 1050535501; i += 2) {
        hurchalla::MontgomeryQuarter<std::uint64_t> mf(i);
        if (IPMRS::case_1050535501_64_2<2>(mf) != primevec[i]) {
            std::cout << "failure 7 on " << i << "\n";
            return 7;
        }
        if (IPMRS::case_1050535501_64_2<1>(mf) != primevec[i]) {
            std::cout << "failure 8 on " << i << "\n";
            return 8;
        }
    }
    std::cout << "case_1050535501_64_2 for odds is ok\n";
      // verify evens other than 2 are not id'd as prime
    for (std::uint64_t i = 4; i < 1050535501; i += 2) {
        // MontgomeryStandardMathWrapper is the only MontgomeryForm
        // that allows an even modulus.  It can do this since it doesn't use
        // any montgomery math- it just wraps standard math in a monty interface
        hurchalla::MontgomeryStandardMathWrapper<std::uint64_t> mf(i);
        if (IPMRS::case_1050535501_64_2<2>(mf) != false) {
            std::cout << "failure 9 on " << i << "\n";
            return 9;
        }
        if (IPMRS::case_1050535501_64_2<1>(mf) != false) {
            std::cout << "failure 10 on " << i << "\n";
            return 10;
        }
    }
    std::cout << "case_1050535501_64_2 for evens is ok\n";
    {  // verify 2 gets id'd as prime
        std::uint64_t i = 2;
        hurchalla::MontgomeryStandardMathWrapper<std::uint64_t> mf(i);
        if (IPMRS::case_1050535501_64_2<2>(mf) != true) {
            std::cout << "failure 11 on " << i << "\n";
            return 11;
        }
        if (IPMRS::case_1050535501_64_2<1>(mf) != true) {
            std::cout << "failure 12 on " << i << "\n";
            return 12;
        }
    }
    std::cout << "case_1050535501_64_2 for 2 is ok\n";


// this bunch all succeeded.
    // ----Test uint64_t miller-rabin, modulus limit < 350269456337, 3 base----

    // This may take a half day to run (I modifed it for multithreading when I
    // ran it, but it's single threaded here for simplicity).
    std::cout << "Very long (~half day) test started...\n";

    for (std::uint64_t i = 3; i < UINT64_C(350269456337); i += 2) {
        using MF = hurchalla::MontgomeryQuarter<std::uint64_t>;
        MF mf(i);
        // By theory, miller-rabin is always correct when it claims a number is
        // composite, so to save time we'll only verify the test result when it
        // claims a number is prime (i.e. we will verify that there are no
        // pseudoprimes).

// we use #if 0 to test only one TRIAL_SIZE, in order to save time on this very
// long running test.
#if 0
        if (IPMRS::case_350269456337_64_3<3>(mf) == true) {
            if (MillerRabinMontgomery<MF, 64, 3, 3>::is_prime(mf) == false) {
                std::cout << "failure 13 on " << i << "\n";
                return 13;
            }
        }
        if (IPMRS::case_350269456337_64_3<1>(mf) == true) {
            if (MillerRabinMontgomery<MF, 64, 3, 3>::is_prime(mf) == false) {
                std::cout << "failure 15 on " << i << "\n";
                return 15;
            }
        }
#endif
        if (IPMRS::case_350269456337_64_3<2>(mf) == true) {
            if (MillerRabinMontgomery<MF, 64, 3, 3>::is_prime(mf) == false) {
                std::cout << "failure 14 on " << i << "\n";
                return 14;
            }
        }
    }
    std::cout << "case_350269456337_64_3 for odds is ok\n";
#if 0
// this would take a very long time, and we know from theory that this function
// will always show that any even number > 2 is composite, because the function
// is implemented with at least one even base.
      // verify evens other than 2 are not id'd as prime
    for (std::uint64_t i = 4; i < UINT64_C(350269456337); i += 2) {
        // MontgomeryStandardMathWrapper is the only MontgomeryForm
        // that allows an even modulus.  It can do this since it doesn't use
        // any montgomery math- it just wraps standard math in a monty interface
        hurchalla::MontgomeryStandardMathWrapper<std::uint64_t> mf(i);
        if (IPMRS::case_350269456337_64_3<3>(mf) != false) {
            std::cout << "failure 16 on " << i << "\n";
            return 16;
        }
        if (IPMRS::case_350269456337_64_3<2>(mf) != false) {
            std::cout << "failure 17 on " << i << "\n";
            return 17;
        }
        if (IPMRS::case_350269456337_64_3<1>(mf) != false) {
            std::cout << "failure 18 on " << i << "\n";
            return 18;
        }
    }
    std::cout << "case_350269456337_64_3 for evens ok\n";
#endif
    {  // verify 2 gets id'd as prime
        std::uint64_t i = 2;
        hurchalla::MontgomeryStandardMathWrapper<std::uint64_t> mf(i);
        if (IPMRS::case_350269456337_64_3<3>(mf) != true) {
            std::cout << "failure 19 on " << i << "\n";
            return 19;
        }
        if (IPMRS::case_350269456337_64_3<2>(mf) != true) {
            std::cout << "failure 20 on " << i << "\n";
            return 20;
        }
        if (IPMRS::case_350269456337_64_3<1>(mf) != true) {
            std::cout << "failure 21 on " << i << "\n";
            return 21;
        }
    }
    std::cout << "case_350269456337_64_3 for 2 is ok\n";


    // ----Test uint64_t miller-rabin, modulus limit < 273919523041, 3 base----

    // This may take a half day to run (I modifed it for multithreading when I
    // ran it, but it's single threaded here for simplicity).
    std::cout << "Second very long (~half day) test started...\n";

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



// this bunch all succeeded.
    // ---Test uint64_t miller-rabin, modulus limit < 47636622961201, 4 base--
    // The implementation of case_47636622961201_64_4 has a base of 2, so we
    // only need to test the strong pseudoprimes to base 2.
    for (auto psp : psps) {
        if (psp >= UINT64_C(47636622961201))
            continue;
        hurchalla::MontgomeryQuarter<std::uint64_t> mf(psp);

        if (IPMRS::case_47636622961201_64_4<4>(mf)) {
            std::cout << "failure 22 on " << psp << "\n";
            return 22;
        }
        if (IPMRS::case_47636622961201_64_4<3>(mf)) {
            std::cout << "failure 23 on " << psp << "\n";
            return 23;
        }
        if (IPMRS::case_47636622961201_64_4<2>(mf)) {
            std::cout << "failure 24 on " << psp << "\n";
            return 24;
        }
        if (IPMRS::case_47636622961201_64_4<1>(mf)) {
            std::cout << "failure 25 on " << psp << "\n";
            return 25;
        }
    }
    std::cout << "case_47636622961201_64_4 is ok\n";


    // ---Test uint64_t miller-rabin, modulus limit < 55245642489451, 4 base--
    // The implementation of case_55245642489451_64_4 has a base of 2, so we
    // only need to test the strong pseudoprimes to base 2.
    for (auto psp : psps) {
        if (psp >= UINT64_C(55245642489451))
            continue;
        hurchalla::MontgomeryQuarter<std::uint64_t> mf(psp);

        if (IPMRS::case_55245642489451_64_4<4>(mf)) {
            std::cout << "failure 22 on " << psp << "\n";
            return 22;
        }
        if (IPMRS::case_55245642489451_64_4<3>(mf)) {
            std::cout << "failure 23 on " << psp << "\n";
            return 23;
        }
        if (IPMRS::case_55245642489451_64_4<2>(mf)) {
            std::cout << "failure 24 on " << psp << "\n";
            return 24;
        }
        if (IPMRS::case_55245642489451_64_4<1>(mf)) {
            std::cout << "failure 25 on " << psp << "\n";
            return 25;
        }
    }
    std::cout << "case_55245642489451_64_4 is ok\n";


// this bunch all succeeded.
    // ---Test uint64_t miller-rabin, modulus limit < 7999252175582851, 5 base--
    // The implementation of case_7999252175582851_64_5 has
    // a base of 2, so we only need to test the strong pseudoprimes to base 2.
    for (auto psp : psps) {
        if (psp >= UINT64_C(7999252175582851))
            continue;
        hurchalla::MontgomeryQuarter<std::uint64_t> mf(psp);

        if (IPMRS::case_7999252175582851_64_5<5>(mf)) {
            std::cout << "failure 26 on " << psp << "\n";
            return 26;
        }
        if (IPMRS::case_7999252175582851_64_5<4>(mf)) {
            std::cout << "failure 27 on " << psp << "\n";
            return 27;
        }
        if (IPMRS::case_7999252175582851_64_5<3>(mf)) {
            std::cout << "failure 28 on " << psp << "\n";
            return 28;
        }
        if (IPMRS::case_7999252175582851_64_5<2>(mf)) {
            std::cout << "failure 29 on " << psp << "\n";
            return 29;
        }
        if (IPMRS::case_7999252175582851_64_5<1>(mf)) {
            std::cout << "failure 30 on " << psp << "\n";
            return 30;
        }
    }
    std::cout << "case_7999252175582851_64_5 is ok\n";


// this bunch all succeeded.
    // ---Test uint64_t miller-rabin, modulus limit < 585226005592931977, 6 base
    // The implementation of case_585226005592931977_64_6 has a base 2, so we
    // only need to test the strong pseudoprimes to base 2.
    for (auto psp : psps) {
        if (psp >= UINT64_C(585226005592931977))
            continue;
        hurchalla::MontgomeryQuarter<std::uint64_t> mf(psp);

        if (IPMRS::case_585226005592931977_64_6<6>(mf)) {
            std::cout << "failure 31 on " << psp << "\n";
            return 31;
        }
        if (IPMRS::case_585226005592931977_64_6<5>(mf)) {
            std::cout << "failure 32 on " << psp << "\n";
            return 32;
        }
        if (IPMRS::case_585226005592931977_64_6<4>(mf)) {
            std::cout << "failure 33 on " << psp << "\n";
            return 33;
        }
        if (IPMRS::case_585226005592931977_64_6<3>(mf)) {
            std::cout << "failure 34 on " << psp << "\n";
            return 34;
        }
        if (IPMRS::case_585226005592931977_64_6<2>(mf)) {
            std::cout << "failure 35 on " << psp << "\n";
            return 35;
        }
        if (IPMRS::case_585226005592931977_64_6<1>(mf)) {
            std::cout << "failure 36 on " << psp << "\n";
            return 36;
        }
    }
    std::cout << "case_585226005592931977_64_6 is ok\n";


    std::cout << "Complete - all tests succeeded.\n";
    return 0;
}
