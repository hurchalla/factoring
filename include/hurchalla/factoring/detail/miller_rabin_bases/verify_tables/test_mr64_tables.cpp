
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

    for (auto psp : psps) {
        hurchalla::MontgomeryFull<decltype(psp)> mf(psp);

// this bunch all succeeded.
        if (is_prime_mr<64, 7, 7>(mf)) {
            std::cout << "failure 1 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 6, 7>(mf)) {
            std::cout << "failure 81 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 5, 7>(mf)) {
            std::cout << "failure 82 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 4, 7>(mf)) {
            std::cout << "failure 83 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 3, 7>(mf)) {
            std::cout << "failure 84 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 7>(mf)) {
            std::cout << "failure 85 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 7>(mf)) {
            std::cout << "failure 86 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded (for both of the ifdef'd tables).
        if (is_prime_mr<64, 6, 6>(mf)) {
            std::cout << "failure 2 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 5, 6>(mf)) {
            std::cout << "failure 3 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 4, 6>(mf)) {
            std::cout << "failure 4 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 3, 6>(mf)) {
            std::cout << "failure 5 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 6>(mf)) {
            std::cout << "failure 6 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 6>(mf)) {
            std::cout << "failure 7 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<64, 5, 5>(mf)) {
            std::cout << "failure 8 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 4, 5>(mf)) {
            std::cout << "failure 9 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 3, 5>(mf)) {
            std::cout << "failure 10 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 5>(mf)) {
            std::cout << "failure 11 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 5>(mf)) {
            std::cout << "failure 12 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<64, 4, 4>(mf)) {
            std::cout << "failure 13 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 3, 4>(mf)) {
            std::cout << "failure 14 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 4>(mf)) {
            std::cout << "failure 15 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 4>(mf)) {
            std::cout << "failure 16 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.  Both for the normal 25KB 2 hash version,
//   and for the ../alternative_tables 32KB 1 hash version
        if (is_prime_mr<64, 3, 3>(mf)) {
            std::cout << "failure 17 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 3>(mf)) {
            std::cout << "failure 18 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 3>(mf)) {
            std::cout << "failure 19 on " << psp << "\n";
            return 1;
        }

// this bunch all succeeded.
        if (is_prime_mr<64, 2, 2>(mf)) {
            std::cout << "failure 20 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 2>(mf)) {
            std::cout << "failure 21 on " << psp << "\n";
            return 1;
        }
    }


    std::cout << "success\n";
    return 0;
}
