// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */


#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/montgomery_arithmetic/MontgomeryForm.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include <cstdint>
#include <fstream>
#include <iostream>
#include <cassert>
#include <string>
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
    // arguments: program name, pseudoprimes file
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << "  SOURCE_FILE_PSEUDOPRIMES_2POW64\n";
        return 1;
    }

    // we will (hopefully) verify that no pseudoprimes get identified as prime.

    std::cout << "started\n";

    namespace hc = hurchalla::detail;

    std::vector<std::uint64_t> psps;
    int status = read_psps(argv[1], psps);
    if (status != 0)
        return 1;

    std::cout << "all psps have been read\n";

#if 1
    // test the bases that cover the full range of uint64_t
    for (auto psp : psps) {
        hurchalla::MontgomeryForm<decltype(psp)> mf(psp);

// this bunch all succeeded.
        if (is_prime_mr<64, 7, 7>(mf)) {
            std::cout << "failure 1 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 6, 7>(mf)) {
            std::cout << "failure 2 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 5, 7>(mf)) {
            std::cout << "failure 3 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 4, 7>(mf)) {
            std::cout << "failure 4 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 3, 7>(mf)) {
            std::cout << "failure 5 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 7>(mf)) {
            std::cout << "failure 6 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 7>(mf)) {
            std::cout << "failure 7 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded (for both of the ifdef'd tables).
        if (is_prime_mr<64, 6, 6>(mf)) {
            std::cout << "failure 8 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 5, 6>(mf)) {
            std::cout << "failure 9 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 4, 6>(mf)) {
            std::cout << "failure 10 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 3, 6>(mf)) {
            std::cout << "failure 11 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 6>(mf)) {
            std::cout << "failure 12 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 6>(mf)) {
            std::cout << "failure 13 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<64, 5, 5>(mf)) {
            std::cout << "failure 14 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 4, 5>(mf)) {
            std::cout << "failure 15 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 3, 5>(mf)) {
            std::cout << "failure 16 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 5>(mf)) {
            std::cout << "failure 17 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 5>(mf)) {
            std::cout << "failure 18 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<64, 4, 4>(mf)) {
            std::cout << "failure 19 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 3, 4>(mf)) {
            std::cout << "failure 20 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 4>(mf)) {
            std::cout << "failure 21 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 4>(mf)) {
            std::cout << "failure 22 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.  Both for the normal 25KB 2 hash version,
//   and for the ../alternative_tables 32KB 1 hash version
        if (is_prime_mr<64, 3, 3>(mf)) {
            std::cout << "failure 23 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 3>(mf)) {
            std::cout << "failure 24 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 3>(mf)) {
            std::cout << "failure 25 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<64, 2, 2>(mf)) {
            std::cout << "failure 26 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 2>(mf)) {
            std::cout << "failure 27 on " << psp << "\n";
            return 1;
        }
    }
#endif


#if 1
    // test the bases that allow a max of 63bit modulus

    std::uint64_t Rdiv2 = (static_cast<std::uint64_t>(1) << 63);

    for (auto psp : psps) {
        if (psp >= Rdiv2)
            continue;
        hurchalla::MontgomeryForm<decltype(psp)> mf(psp);

// this bunch all succeeded.
        if (is_prime_mr<63, 6, 6>(mf)) {
            std::cout << "failure 28a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 5, 6>(mf)) {
            std::cout << "failure 29a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 4, 6>(mf)) {
            std::cout << "failure 30a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 3, 6>(mf)) {
            std::cout << "failure 31a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 2, 6>(mf)) {
            std::cout << "failure 32a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 1, 6>(mf)) {
            std::cout << "failure 33a on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<63, 5, 5>(mf)) {
            std::cout << "failure 34 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 4, 5>(mf)) {
            std::cout << "failure 35a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 3, 5>(mf)) {
            std::cout << "failure 36a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 2, 5>(mf)) {
            std::cout << "failure 37a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 1, 5>(mf)) {
            std::cout << "failure 38a on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<63, 4, 4>(mf)) {
            std::cout << "failure 39a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 3, 4>(mf)) {
            std::cout << "failure 40a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 2, 4>(mf)) {
            std::cout << "failure 41a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 1, 4>(mf)) {
            std::cout << "failure 42a on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<63, 3, 3>(mf)) {
            std::cout << "failure 43a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 2, 3>(mf)) {
            std::cout << "failure 44a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 1, 3>(mf)) {
            std::cout << "failure 45a on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<63, 2, 2>(mf)) {
            std::cout << "failure 46a on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<63, 1, 2>(mf)) {
            std::cout << "failure 47a on " << psp << "\n";
            return 1;
        }
    }
#endif


#if 1
    // test the bases that allow a max of 62bit modulus

    std::uint64_t Rdiv4 = (static_cast<std::uint64_t>(1) << 62);

    for (auto psp : psps) {
        if (psp >= Rdiv4)
            continue;
        hurchalla::MontgomeryQuarter<decltype(psp)> mf(psp);

// this bunch all succeeded.
        if (is_prime_mr<62, 6, 6>(mf)) {
            std::cout << "failure 28 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 5, 6>(mf)) {
            std::cout << "failure 29 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 4, 6>(mf)) {
            std::cout << "failure 30 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 3, 6>(mf)) {
            std::cout << "failure 31 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 2, 6>(mf)) {
            std::cout << "failure 32 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 1, 6>(mf)) {
            std::cout << "failure 33 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<62, 5, 5>(mf)) {
            std::cout << "failure 34 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 4, 5>(mf)) {
            std::cout << "failure 35 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 3, 5>(mf)) {
            std::cout << "failure 36 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 2, 5>(mf)) {
            std::cout << "failure 37 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 1, 5>(mf)) {
            std::cout << "failure 38 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<62, 4, 4>(mf)) {
            std::cout << "failure 39 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 3, 4>(mf)) {
            std::cout << "failure 40 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 2, 4>(mf)) {
            std::cout << "failure 41 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 1, 4>(mf)) {
            std::cout << "failure 42 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<62, 3, 3>(mf)) {
            std::cout << "failure 43 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 2, 3>(mf)) {
            std::cout << "failure 44 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 1, 3>(mf)) {
            std::cout << "failure 45 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<62, 2, 2>(mf)) {
            std::cout << "failure 46 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<62, 1, 2>(mf)) {
            std::cout << "failure 47 on " << psp << "\n";
            return 1;
        }
    }
#endif


    std::cout << "success\n";
    return 0;
}
