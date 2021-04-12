I generated all the miller-rabin hash tables myself (via brute-force searches, the details of which are beyond the scope of discussion here).
I verified their correctness using the verification scripts/cpp files in this folder.  If you wish to verify the hash tables and bases yourself, you can view and run the scripts/cpp files in this folder.  You will need a file that lists the strong base2 psuedoprimes for all numbers < (1<<64):

You can get a two-part 7zip archive that contains the desired file, "strong_psps_to_2_64.txt", from
https://raw.githubusercontent.com/hurchalla/Extra-Files/main/psps.7z.001
https://raw.githubusercontent.com/hurchalla/Extra-Files/main/psps.7z.002

I programmatically created the "strong_psps_to_2_64.txt" file using Jan Feitsma's Fermat base-2 pseudoprimes for all numbers < (1<<64), available at http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html  and see also http://www.janfeitsma.nl/math/psp2/index


For miller-rabin bases designed to primality test numbers < (1<<32), the verifications work by first creating a bit vector via the sieve of Eratosthenes, which serves as a lookup table to check primality of any uint32_t number.  Next, the verification loops over every possible uint32_t or uint16_t (for a 16 bit test) value, and for each value calls MillerRabinMontgomery::is_prime() (from is_prime_miller_rabin.h) to exercise the targetted bases from the folder miller_rabin_bases.  On each loop iteration, the verification compares the miller-rabin primality result with the bit vector primality lookup result.

For miller-rabin bases designed to primality test numbers < (1<<64),

There is no verification test for the miller-rabin bases designed for probabilistic primality testing of numbers < (1<<128), since there's no known way to truly verify this range.  Some testing of various primes and composites in this range does exist in the unit tests in the tests folder at the base of this repository, but that's not really verification.  We should keep in mind that even if a failure were were found as a result of using the bases for 128bit numbers, it would not exactly show the bases are wrong, since they were designed for probabilistic testing rather than deterministic testing.  Even so, assuming the expected probability of any failure of a primality test using these bases is even roughly close to realistic, failure due to the bases (rather than a programming bug) would be essentially an impossibility in any context other than theory.


I verified this using the Feitsma database psps-below-2-to-64.txt
from http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html
I first performed a miller rabin base 2 test on the database
entries to trim the database down to only strong pseudoprimes
base 2, which resulted in about 32 million entries.  I then ran
miller rabin tests on each entry using the remaining 6 bases,
and together those 6 bases correctly identified every entry as
composite.
So assuming that the Feitsma database is correct, I found that
the 7 bases correctly determine primality for all numbers < 2^64.
For additional info, see http://www.janfeitsma.nl/math/psp2/index


Note that the validity of all the 64bit and 62bit verification tests in this folder depends upon the correctness of the Feitsma base-2 Fermat pseudoprimes database.  I haven't attempted to verify the database, but at least on the surface it appears to be high quality work.  It also has been available for 8 years and I haven't seen any information to suggest it contains any mistakes.  This of course is not the same thing as an independent confirmation of correctness.  Additionally, as far as I know the work has not undergone peer review for an article publication.  Nevertheless I have assumed it is correct, and I found no disagreements with it when using trial-division primality testing over a very limited range of numbers (running for a few days, yet covering less than a millionth of the entire range of uint64_t).
