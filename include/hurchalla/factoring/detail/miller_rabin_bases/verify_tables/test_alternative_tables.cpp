
// Unfortunately, the only way to test MillerRabinBases64_3_alt.h is to define
// the macro below, which causes is_prime_miller_rabin.h to #include the
// desired alternative header, and not to #include of MillerRabinBases64_3.h.
// The reason why this is necessary is that both those headers have an identical
// class name of MillerRabinBases<64, 3, DUMMY>, which would cause a compile
// error or undefined behavior (due to ODR) if both were included.

#define HURCHALLA_CHOOSE_ALTERNATE_MILLER_RABIN_BASES64_3 1


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
        if (is_prime_mr<64, 3, 3>(mf)) {
            std::cout << "failure 5 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 2, 3>(mf)) {
            std::cout << "failure 6 on " << psp << "\n";
            return 1;
        }
        if (is_prime_mr<64, 1, 3>(mf)) {
            std::cout << "failure 7 on " << psp << "\n";
            return 1;
        }
    }

    std::cout << "success\n";
    return 0;
}
