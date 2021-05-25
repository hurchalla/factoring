# Factoring

A factoring and primality testing library for C++, optimized for 64 and 32 bit integers and supporting sizes up to 128 bits.

## Design goals

A correct and flexible library with best possible performance for factoring and primality testing of 64 and 32 bit integer types.  128 bit integer factoring and primality testing is also supported, but note that other algorithms will usually be more suitable for values much above 2^64.

## How to use the library

### With CMake

If you have your own CMake project and you want to add this factoring library, it's easy.  Clone this repository to a folder on your system named 'factoring' (or some other name if you prefer).  In your project's CMakeLists.txt file, add the line:  
add_subdirectory(factoring  PATH_TO_THIS_CLONED_FACTORING_REPOSITORY/factoring)

### (Mostly) Without CMake

If you're not using CMake for your project, you'll need to install the factoring sources and dependencies to some location so that you can use them.  To do so, clone this repository to a folder on your system named 'factoring' (or another name if you prefer).  If you don't already have CMake on your system you'll need to install it, and then from your shell run the following commands:  

>cd *path_to_the_root_of_the_cloned_factoring_repo*  
>mkdir tmp  
>cd tmp  
>cmake -S.. -B.  
>cmake --install . [--prefix *the_folder_you_want_to_install_to*]  


### Include paths

### Example

## Algorithms used

For factoring: Pollard-Rho Brent (https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants and "An Improved Monte Carlo Factorization Algorithm" by Richard Brent),
preceded by a stage of trial division using special algorithms for divisibility (Section 10-17 from Hacker's Delight 2nd edition by Henry Warren, and "ALGORITHM A: IS_DIV_A" from "Efficient long division via Montgomery" by Ernst W. Mayer).  

For primality testing: Deterministic Miller-Rabin (https://miller-rabin.appspot.com/ and https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants).  If using the resource_intensive_api, also Sieve of Eratosthenes for 32 bit and smaller types.

## Status
Released.  All planned functions and unit tests are complete and work correctly.

## Miscellaneous
If you're interested in experimenting, predefining certain macros when compiling can improve performance - see [macros_for_performance.md](macros_for_performance.md).
