The classes in this folder are used by "is_prime_miller_rabin.h" to provide fast deterministic miller-rabin primality tests.  The header files use a naming scheme of MillerRabinBasesX_Y.  X indicates the largest integer type that can be tested with the bases in a class (e.g. 32 for uint32_t, 64 for uint64_t; the exception is 62 which indicates uint64_t restricted to values < (1<<62)).  Y indicates the total number of miller-rabin bases that the class uses.  Roughly speaking, the fewer the bases, the faster the class performs at identifying that a prime number is prime.  The total number of bases generally makes no difference in the performance of a class at identifying the primality of a composite number, with the rare exception of when a composite is a pseudoprime (pseudoprimes are much less common than primes).  The downside of using a class with fewer bases is that reducing the number of bases requires the class to use larger hash tables.  Large hash tables increase executable size and memory usage, and every access to a table can cause a CPU cache eviction, which can detract from program performance.  When you need to perform intensive primality testing, the classes with the fewest bases (and hence large tables) are usually worth using.  For other cases it may be hard to justify the large tables, and classes with a larger number of bases can be a better idea.  The classes with the largest number of bases use no hash table at all.

Note: you can find an alternative uint64_t 3 base (32KB hash table) class in the alternative_tables folder.  It uses a larger hash table than the version in the mainline miller_rabin_bases folder, and so it's generally not preferred.  However, it reads only one cache line per table access (the normal version reads two lines per access), and it uses a smaller and faster hash function - these differences are unlikely to be as important as the fact that it uses a larger hash table, but perhaps in some circumstance it could make it preferable.  To use it, predefine a C++ macro called HURCHALLA_CHOOSE_ALTERNATE_MILLER_RABIN_BASES64_3 when compiling the sources (typically you predefine a macro via a compiler command line option -D<macroname>).


For good background information, and also to see some different miller-rabin hash tables, see  https://miller-rabin.appspot.com/  particularly the Hashing section.
The website has great links to work on miller-rabin bases/hashes by (in no particular order) Jan Feitsma and William Galway, Bradley Berg, Michal Forisek and Jakub Jancina, Dana Jacobsen, Wojciech Izykowski and Marcin Panasiuk, Steve Worley, Gerhard Jaeschke, Jim Sinclair, and more.


I don't believe there would be much benefit to trying to generate better tables for primality testing uint32_t/uint64_t values than the ones here, since I expect it would be extremely difficult to find tables that are even 5-10% smaller.  The possible exceptions are the files MillerRabinBases62_3.h, MillerRabinBases64_4.h, and MillerRabinBases64_5.h, due to their use of equal size tables.  I expect it perhaps might be possible to reduce the size of tables used in MillerRabinBases64_4 or MillerRabinBases62_3 by as much as 20-25% (assuming that level of difference merits attention) by generating unequal size tables as described below.  I would expect less improvement for MillerRabinBases64_5.  Keep in mind that as also described below, I chose equal size tables in order to get the minor advantages equal size tables provide, usually when the size difference wouldn't be large.
If you wish to generate your own miller-rabin bases/tables for primality testing uint64_t values, you might find the following files interesting:
https://github.com/hurchalla/Extra-Files/blob/9db1c93fe912f39ad0c23c445750ec0826fbf4ac/bases_2.txt
https://github.com/hurchalla/Extra-Files/blob/9db1c93fe912f39ad0c23c445750ec0826fbf4ac/bases_2_15.txt
https://github.com/hurchalla/Extra-Files/blob/9db1c93fe912f39ad0c23c445750ec0826fbf4ac/base_combos_2_15.txt

If you'd still like to generate tables, here's some basic descriptions of what I did, and/or what seemed notable.
-----
For the tables designed for uint32_t and smaller deterministic primality testing:
First, for every uint16_t value>=2, I used each value as a miller-rabin base and calculated the uint32_t strong pseudoprimes for that base.  Upon completion of all bases I stored all the pseudoprimes (with corresponding base number) to file.  Next in a separate program, I looped through variations of simple hash functions involving left/right shifts and an xor, and for each attempted hash function I looped through the hash buckets, trying via brute-force to find a base for each bucket that correctly identified the primality of all uint32_t numbers that hash to the bucket.  The stored strong pseudoprimes made this fast, since I only needed to test that no pseudoprimes (for the base being attempted) mapped to the bucket.  If I found a good base for every bucket, those bases made up a finished hash table for use with the hash function.  If the program failed to find a good base for even one bucket, it started over using a different hash function.  The computationally intensive aspect was trying to successfully generate the smallest possible hash table (using hash functions with the fewest buckets possible), since the chance of getting all good bases decreases as we reduce the number of total buckets, and getting a good base takes much longer per bucket due to the decreased chance of any given base working for a bucket.  I noticed a fairly sharply edge between sizes of table attempts, where table generation would either easily succeed or always fail; for example a table with 144 buckets might be an easy success for nearly all hash functions, whereas trying 128 buckets might fail for all hash functions I tried.  It seemed there was about a 10-15% range above what we could imagine to be the smallest achievable table, in which generating a successful table went from very hard (but possible) to becoming easy, and anything larger was both fast and easy. 
I found there was little or no benefit to using "good" hash functions.  Very simple hash functions were just as good, using primarily left/right shifts.  The pseudoprime values themselves have enough of a random characteristic that they hash well across all buckets when employing left/right shifts, although I found for the hash tables designed for uint32_t that I needed to use a single xor in order to increase the randomness enough to get an excellent mapping across hash buckets.

For the tables designed for uint64_t deterministic primality testing:
In practice for us to generate miller-rabin bases for uint64_t primality testing, we must rely upon the Feitsma base-2 pseudoprimes database (see ./verify_tables/README.txt), in order to make the task computationally feasible.  Accordingly, we must always use the value 2 for our first base.  Thus we're required to use at least two total bases, since the first base is fixed at 2.
For the hash table for two total bases, I needed to find a hash table for the second base.  Similarly to generating hash tables for uint32_t, I looped through variations of simple hash functions, though it was possible to use even simpler hash functions (no xor needed) due to the inherent randomness of the base-2 pseudoprimes (which survive miller-rabin testing of the first base == 2) combined with the randomness of the pseudoprimes to any second base.  The surviving pseudoprimes to both bases tend to map out pretty evenly across hash buckets, and we look for the lucky occurence where a base being attemped does not have any psuedoprimes (to both bases) that map to the hash bucket being attempted.  The rest of the process was similar to the uint32_t description above.
For targeting three or more bases, it was (as always) required to use the number 2 for base[0] in the array of total bases.  I found it typically worked best to generate two hash tables rather than one. [I give credit to Bradley Berg's miller-rabin hashing code at techneon.com for the idea.]  I found that using a second hash table somewhere around 16x the size of the first table tended to succeed at creating the smallest possible combined hash tables.  When I generated one single hash table instead (with a fixed second base instead of a second hash table), the best table size I could achieve was almost always larger than the former approach.  Likewise, when I generated two equal sized hash tables (rather than dissimilar sizes), the combined best table size I could achieve tended to be larger than the first approach, yet still notably smaller than for a single hash table.  Perhaps a size difference around 8-24x might be a good range to explore for getting smallest combined hash table results; another range to explore might be second table sizes that are around the square of the first table's size.  Sometimes I preferred to use two equal size tables (primarily when there was little no difference in achieved table sizes), since two equal size tables need only one single hash function and allow me to store corresponding elements of both tables adjacent to each other in memory (via 2D array); this results in a lookup to the two tables needing to touch only one cache line, as opposed to a lookup to two unequal size arrays which almost always must touch two cache lines.  Regardless, the second hash table had to be a multiple of the size of the first table, and any value that mapped to a particular bucket of the second table had to always map to the same bucket of the first table.  Another way to look at this is that all values that hash to a particular bucket of the first table must hash to exactly one distinct subset of the buckets of the second table (with the total number of subsets equal to the total number of buckets in table1, and with every subset having size totalBucketsTable2/totalBucketsTable1).  As an example, a pair of extremely simple hash functions that satisfy these demands are:  bucket1 = (uint32_t)num >> 29  and  bucket2 = (uint32_t)num >> 25.  This example has 8 buckets for the first table and 128 buckets for the second, and any numbers that map to a particular second bucket always map to the same first bucket.  Also, to be even more specific, using the example hash functions, any value that maps to a bucket1 == 0 (out of the 8) can only map to possible bucket2 values 0 through 15 (out of the 128).
Here's the process I used to generate tables to achieve exactly 3 total bases, with two hash tables and the first base fixed at 2 (as always):
The basic idea we can use is to try a particular base for a particular bucket1 (for candidate bases I brute-force looped through all uint16_t values), and run miller-rabin trials with that base on all the base-two strong pseudoprimes that hash to that particular bucket1 (using the first hash function).  For any and all base-two psuedoprimes that map to bucket1 and survive the bucket1 base miller-rabin trial we just performed, we put the survivor into a new base1-bucket1-pseudoprimes list; and that list always tends to be quite small compared to the original base-two pseudoprimes list.  Then we do basically the same thing using this new pseudoprimes list, for a particular base2 for a particular bucket2, so long as all values that map to that bucket2 also map to the particular bucket1 that we are currently working on.  What we are trying to find is a base1 value for a bucket1, such that for every bucket2 that the bucket1 corresponds to, there exists at least one base2 value that will result in an empty pseudoprimes list (after filtering the original base-two pseudoprimes by hashing to bucket1, miller rabin testing with base1, and doing the same for bucket2/base2).  If we find a good base1 for a bucket1, we automatically have gotten the good base2s for the corresponding bucket2s in the process.

It may help to view code I wrote for the process:

**Disclaimer- I'd intended to essentially discard this generator code after I'd used it to generate the miller-rabin hash tables that are in this folder; it's quick/ugly for getting a job done.  But maybe it will be useful as a way to illustrate the process.  To make the code simpler to understand, I removed optimizations that sped it up 4x.

See the comments at the beginning of main() for info on the command line arguments.


// This code is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository

// This code uses my github repos: factoring, modular_arithmetic, and util
#include "hurchalla/factoring/detail/is_prime_miller_rabin.h"
#include "hurchalla/montgomery_arithmetic/montgomery_form_aliases.h"
#include "hurchalla/util/compiler_macros.h"

#include <vector>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <limits>

// these correspond to uint16_t value limits
#define BASE1_LIMIT 65536
#define BASE2_LIMIT 65536


template <typename T, typename B>
HURCHALLA_FORCE_INLINE
bool miller_rabin_trial(T number, B base)
{
    hurchalla::MontgomeryFull<T> mf(number);
    const std::array<T, 1> bases = { base };
    return hurchalla::detail::miller_rabin_trials<1>(mf, bases);
}

uint32_t get_hash_bucket2(uint64_t psp, uint64_t mask2, uint64_t multiplier, int shift)
{
    return static_cast<uint32_t>(((psp & mask2) * multiplier) >> shift);
}
uint32_t get_hash_bucket1(uint32_t hash_bucket2, uint32_t mask1, uint32_t shift1, bool use_DIT_form)
{
    if (use_DIT_form)
        return (hash_bucket2 & mask1);
    else
        return (hash_bucket2 >> shift1);
}

int generate_dual_hash_tables(
                 std::ofstream& out_file,
                 uint32_t total_buckets1,
                 uint32_t total_buckets2,
                 const std::vector<uint64_t>& pseudoprimes,
                 int shift_start = 6,
                 int shift_end = 40,
                 int shift_increment = 1)
{
    for (auto psp : pseudoprimes) {
        if (psp % 2 == 0) {
            std::cerr << "ERROR: unexpected even psp\n";
            return 1;
        }
    }
    if (total_buckets1 > total_buckets2) {
        std::cerr << "ERROR: total_buckets1 must be <= total_buckets2\n";
        return 2;
    }
    if (total_buckets1 == 0 || total_buckets2 == 0) {
        std::cerr << "ERROR: invalid zero value for total_buckets1 or total_buckets2\n";
        return 3;
    }
    if (total_buckets2 % total_buckets1 != 0) {
        std::cerr << "ERROR: total_buckets1 doesn't divide total_buckets2\n";
        return 4;
    }

    uint32_t bucket_ratio = (total_buckets2/total_buckets1);

    uint32_t multiplier = total_buckets2;
    uint32_t pow2exponent = 0;
    while (multiplier % 2 == 0) {
        ++pow2exponent;
        multiplier = static_cast<uint32_t>(multiplier/2);
    }

    uint32_t mask1 = static_cast<uint32_t>(total_buckets1 - 1); // ignored when we don't use DIT form
    uint32_t shift1 = 0;                                        // ignored when we do use DIT form
    bool use_DIT_form = false;                                     
    for (uint32_t n=1; n<bucket_ratio; n*=2) {
        ++shift1;
        if (n*2 > bucket_ratio)  // if bucket_ratio is not a power of 2, we must use DIT form
            use_DIT_form = true;
    }
    if (use_DIT_form) {
        // if we must use DIT form, then total_buckets1 must be a power of 2
        for (uint32_t n=1; n<total_buckets1; n*=2) {
            if (n*2 > total_buckets1) {  // if total_buckets1 is not a power of 2
                std::cerr << "ERROR: either bucket_ratio or total_buckets1 must be a power of 2\n";
                return 5;
            }
        }
    }

    std::vector<uint32_t> successful_bucketbases1(total_buckets1, 0);
    std::vector<uint32_t> successful_bucketbases2(total_buckets2, 0);
    std::vector<std::vector<uint64_t>> bucket1_psps(total_buckets1);
    std::vector<std::vector<uint64_t>> bucket2_psps(total_buckets2);


for (int shift = shift_start; shift <= shift_end; shift += shift_increment) {
    std::cout << "trying shift " << shift << "\n";
    uint64_t mask2 = (static_cast<uint64_t>(1) << (pow2exponent + shift)) - 1u;

    for (uint32_t bucketnum1 = 0; bucketnum1 < total_buckets1; ++bucketnum1)
        bucket1_psps[bucketnum1].clear();
    for (auto psp : pseudoprimes) {
        uint32_t hash_bucket2 = get_hash_bucket2(psp, mask2, multiplier, shift);
        uint32_t hash_bucket1 = get_hash_bucket1(hash_bucket2, mask1, shift1, use_DIT_form);
        bucket1_psps[hash_bucket1].push_back(psp);
    }

    bool complete_success = true;
    for (uint32_t bucketnum1 = 0; bucketnum1 < total_buckets1; ++bucketnum1) {
        bool got_good_bucketbase1 = false;
        for (uint32_t base1 = 3; base1 < BASE1_LIMIT; ++base1) {
            for (uint32_t bucketnum2 = 0; bucketnum2 < total_buckets2; ++bucketnum2)
                bucket2_psps[bucketnum2].clear();
            for (auto psp : bucket1_psps[bucketnum1]) {
                bool isProbablePrime = miller_rabin_trial(psp, base1);
                if (isProbablePrime) {
                    uint32_t hash_bucket2 = get_hash_bucket2(psp, mask2, multiplier, shift);
                    bucket2_psps[hash_bucket2].push_back(psp);
                }
            }
            bool subset_buckets2_good = true;

            uint32_t bucket2start = bucketnum1;
            uint32_t bucket2limit = total_buckets2;
            uint32_t bucket2increment = total_buckets1;
            if (!use_DIT_form) {
                bucket2start = bucketnum1 * bucket_ratio;
                bucket2limit = (bucketnum1 + 1) * bucket_ratio;
                bucket2increment = 1;
            }
            for (uint32_t bucketnum2 = bucket2start; bucketnum2 < bucket2limit; bucketnum2 += bucket2increment) {
                bool got_good_bucketbase2 = false;
                for (uint32_t base2 = 3; base2<BASE2_LIMIT; ++base2) {
                    bool success = true;
                    for (auto psp : bucket2_psps[bucketnum2]) {
                        bool isProbablePrime = miller_rabin_trial(psp, base2);
                        if (isProbablePrime) {
                            success = false;
                            break;
                        }
                    }
                    if (success) {
                        got_good_bucketbase2 = true;
                        successful_bucketbases2[bucketnum2] = base2;
                        break;
                    }
                }
                if (!got_good_bucketbase2) {
                    subset_buckets2_good = false;
                    break;
                }
            }
            if (subset_buckets2_good) {
                got_good_bucketbase1 = true;
                successful_bucketbases1[bucketnum1] = base1;
                break;
            }
        }
        if (!got_good_bucketbase1) {
            complete_success = false;
            break;
        }
    }
    if (complete_success) {
        std::cout << "success\n";
        out_file << "\n" << " shift " << shift << "\n";
        for (uint32_t bucketnum1 = 0; bucketnum1 < total_buckets1; ++bucketnum1)
            out_file << bucketnum1 << " " << successful_bucketbases1[bucketnum1] << "\n";
        for (uint32_t bucketnum2 = 0; bucketnum2 < total_buckets2; ++bucketnum2)
            out_file << bucketnum2 << " " << successful_bucketbases2[bucketnum2] << "\n";
    }
}
    return 0;
}



int main(int argc, char* argv[])
{
    // ********
    // To get the strong base-2 pseudoprimes file (SRCFILE_STRONG_BASE2_PSPS below),
    // download a two-part 7zip archive that contains the desired file,
    // "strong_psps_to_2_64.txt", from
    //    https://raw.githubusercontent.com/hurchalla/Extra-Files/main/psps.7z.001
    //    https://raw.githubusercontent.com/hurchalla/Extra-Files/main/psps.7z.002
    // For more details on this file, please read ./verify_tables/README.TXT
    // ********

    // Generally some good choices are SHIFT_START=10, SHIFT_END=35, SHIFT_INCREMENT=1.
    // Each different shift value results in a different hash function being used for table
    // generation.  If you want to evaluate only a single hash function to save time,
    // setting both SHIFT_START and SHIFT_END to a value around 22 seems to have some of the
    // best probabilities for table generation success.
    // SHIFT_INCREMENT exists to make it easier for you to run multiple instances of this
    // program in parallel on multiple CPU cores, without any overlap of work done.
    // For example, you could run the program in two two different shells (and hence
    // two different processes), and choose SHIFT_INCREMENT=2, and use SHIFT_START=10
    // for the first shell process, and SHIFT_START=11 for the second.

    // NUM_BUCKETS1 is the size (number of buckets) you wish to use for the first hash table.
    // NUM_BUCKETS2 is the size (number of buckets) you wish to use for the second hash table.
    //    NUM_BUCKETS2 must be a multiple of NUM_BUCKETS1.  And either NUM_BUCKETS1 must be
    //    a power of 2, or NUM_BUCKETS2/NUM_BUCKETS1 must be a power of 2.

    // The FIXED_BASE arguments are optional but important-
    // (For 3 bases) If you do not specify any FIXED_BASE arguments, then the two tables you
    // generate will be for deterministic miller-rabin testing of 64bit values using three
    // bases: base[0]=2, base[1]=table1_base, base[2]=table2_base.
    // (For 4 bases) If you specify only FIXED_BASE1, the two tables you generate will be for
    // four bases: base[0]=2, base[1]=FIXED_BASE1, base[2]=table1_base, base[3]=table2_base.
    // (For 5 bases) Specifying FIXED_BASE1 and FIXED_BASE2 will result in the two tables
    // being designed for five bases (in similar fashion).
    // (For 6 bases) Specifying FIXED_BASE1 and FIXED_BASE2 and FIXED_BASE3 will result in
    // two tables for six bases.

    if (argc < 8 || argc > 11) {
        std::cerr << "Usage: " << argv[0] << "  OUT_FILE  SRCFILE_STRONG_BASE2_PSPS"
                                          << "  NUM_BUCKETS1  NUM_BUCKETS2"
                                          << "  SHIFT_START  SHIFT_END  SHIFT_INCREMENT"
                                          << "  [FIXED_BASE1  FIXED_BASE2  FIXED_BASE3]\n";
        return 1;
    }

    static_assert(std::numeric_limits<unsigned long long>::digits ==
                  std::numeric_limits<uint64_t>::digits, "");
    uint64_t num_buckets1 = std::stoull(argv[3]);
    uint64_t num_buckets2 = std::stoull(argv[4]);
    uint64_t shift_start = std::stoull(argv[5]);
    uint64_t shift_end = std::stoull(argv[6]);
    uint64_t shift_increment = std::stoull(argv[7]);

    std::vector<uint64_t> fixed_bases;
    if (argc > 8)
        fixed_bases.push_back(std::stoull(argv[8]));
    if (argc > 9)
        fixed_bases.push_back(std::stoull(argv[9]));
    if (argc > 10)
        fixed_bases.push_back(std::stoull(argv[10]));

    std::ofstream out_file(argv[1]);
    if(!out_file) {  
        std::cerr << "Error: file open failed for " << argv[1] << "\n";
        return 2;
    }

    std::vector<uint64_t> surviving_pseudoprimes;
    {
        std::vector<uint64_t> base2_strong_psps;
        {
            std::ifstream in_file(argv[2]);
            if(!in_file) {
                std::cerr << "Error: file open failed for " << argv[2] << "\n";
                return 2;
            }
            uint64_t a;
            while (in_file >> a)
                base2_strong_psps.push_back(a);
            if (in_file.bad()) {
                std::cerr << "Error: I/O error while reading input file\n";
                return 4;
            }
            else if (in_file.eof()) {
            }
            else if (in_file.fail()) {
                std::cerr << "Error: non-integer data found in input file\n";
                return 5;
            }
            else {
                std::cerr << "Error: unknown file read error\n";
                return 6;
            }
        }
        for (auto psp : base2_strong_psps) {
            bool isProbablePrime = true;
            for (auto base : fixed_bases) {
                if (!miller_rabin_trial(psp, base))
                    isProbablePrime = false;
            }
            if (isProbablePrime)
                surviving_pseudoprimes.push_back(psp);
        }
    }

    generate_dual_hash_tables(
                 out_file,
                 static_cast<uint32_t>(num_buckets1),
                 static_cast<uint32_t>(num_buckets2),
                 surviving_pseudoprimes,
                 static_cast<int>(shift_start),
                 static_cast<int>(shift_end),
                 static_cast<int>(shift_increment));
    return 0;
}
