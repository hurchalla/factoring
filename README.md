# Factoring

A factoring and primality testing library for C++17, optimized for 64 and 32 bit integers and supporting sizes up to 128 bits.  This is a header-only library, designed for correctness and ideal performance when factoring or primality testing 64 bit and 32 bit integers.  128 bit integer factoring and primality testing is also supported, but you should note that other algorithms would be more suitable for values much above 2^64.  The algorithms used in this library are described further below.

This library requires a compiler that supports C++17 (if you are not using CMake, you may need to specify the option -std="c++17" when compiling).
For good performance you *must* ensure that the standard macro NDEBUG (see <cassert>) is defined when compiling.

## How to use the library

### With CMake

If you have your own CMake project and you want to add this factoring library, then clone this git repository to a folder on your system named 'factoring' (or another name if you prefer).  In your project's CMakeLists.txt file, add the line:  
add_subdirectory(factoring  PATH_TO_THIS_CLONED_FACTORING_REPOSITORY/factoring)

For good performance you *must* ensure that the standard macro NDEBUG (see <cassert>) is defined when compiling.

### Without CMake (mostly)

If you're not using CMake for your project, you'll need to install/copy these factoring headers and dependencies to some location so that you can use them.  To do so, clone this git repository to a folder on your system named 'factoring' (or another name if you prefer).  If you don't already have CMake on your system you'll need to install CMake, and then from your shell run the following commands:  

>cd *path_to_the_root_of_the_cloned_factoring_repo*  
>mkdir tmp  
>cd tmp  
>cmake -S.. -B.  
>cmake --install . --prefix *the_folder_you_want_to_install_to*  

This will copy all the header files needed for the factoring library to the installation folder of your choosing.
When compiling your project, you'll of course need to ensure that you have the installed folder as part of your include path.  

For good performance you *must* ensure that the standard macro NDEBUG (see <cassert>) is defined when compiling.

## The API

The API consists of five header files in total (all the headers that are not under the detail folder).  These files are the three general purpose files factorize.h, is_prime.h, greatest_common_divisor.h, and in the resource_intensive_api folder, the two special purpose files factorize_intensive_uint32.h and IsPrimeIntensive.h.  Please view these files for their documentation.  A quick summary of the functions is provided below; in all cases T is a template parameter of unsigned integral type.  

factorize(T x, int& num_factors).  Returns a std::array containing the factors of x.  
factorize_to_vector(T x).  Returns a std::vector containing the factors of x.  
greatest_common_divisor(T a, T b).  Returns the greatest common divisor of a and b.  
is_prime(T x).  Returns true if x is prime.  Otherwise returns false.  

(from the resource_intensive_api folder)  
IsPrimeIntensive(T x).  This is a functor that returns true if x is prime.  Otherwise returns false.  Depending on type T, this functor can use a very large amount of memory and can take many seconds to construct.  See IsPrimeIntensive.h for details.  
factorize_intensive_uint32(uint32_t x, int& num_factors, const IsPrimeIntensive<uint32_t,true>& ipi).  Returns a std::array containing the factors of x.  Note that the IsPrimeIntensive argument will usually take many seconds to construct and will use a large amount of memory.  See IsPrimeIntensive.h for details.  

## Example code

## Algorithms used

For factoring: Pollard-Rho Brent (https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants and "An Improved Monte Carlo Factorization Algorithm" by Richard Brent),
preceded by a stage of trial division using special algorithms for divisibility (Section 10-17 from Hacker's Delight 2nd edition by Henry Warren, and "ALGORITHM A: IS_DIV_A" from "Efficient long division via Montgomery" by Ernst W. Mayer).  

For primality testing: Deterministic Miller-Rabin (https://miller-rabin.appspot.com/ and https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants).  If using the resource_intensive_api, also Sieve of Eratosthenes for 32 bit and smaller types (https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes).

Special attention was paid to instruction level parallelism, to take advantage of typical pipelined/superscalar CPUs.  This resulted in modifications to some algorithms, particularly Miller-Rabin (and its associated modular exponentiation) to perform more than one trial/exponentiation at a time.

Interesting and near-optimal hash tables were generated for fast deterministic Miller-Rabin primality testing.  You can view the [README.md](https://github.com/hurchalla/factoring/blob/master/include/hurchalla/factoring/detail/miller_rabin_bases/README.TXT), and/or see the header files in the folder detail/miller_rabin_bases.  The general purpose is_prime and factoring functions use some of the smallest of these hash tables (8 to 100 byte), to minimize memory footprint and cache impact while still receiving a performance boost.  The resource intensive api functions use the largest of these hash tables.

The Pollard-Rho-Brent (and Pollard-Rho) algorithm benefited from a simple addition that I haven't seen elsewhere.  The addition is a "pre-loop" that advances as quickly as possible through a portion of the initial pseduo-random sequence before beginning the otherwise normal Pollard-Rho Brent algorithm.  Every Pollard-Rho pseudo-random sequence begins with a non-periodic segment, and trying to extract factors from that segment is mostly wasted work since the algorithm logic relies on a periodic sequence.  Using a "pre-loop" that does nothing but iterate for a set number of times through the sequence thus improves performance on average, since it quickly moves past some of the unwanted non-periodic segment.

## Status
Released, version 1.0.

## Miscellaneous
If you're interested in experimenting, predefining certain macros when compiling can improve performance - see [macros_for_performance.md](macros_for_performance.md).
