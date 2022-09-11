# The EPR Factoring Library

![Alt text](images/cycle6.jpg?raw=true "Greek rho cycle")

EPR is a high performance, easy to use factoring and primality checking C++ library (header-only) for any integer up to 128 bits in size.  At the time of this writing, EPR can provide you with one of the fastest factoring functions known for 64 bit integers (i.e. types int64_t and uint64_t).  Note that for good performance you *must* ensure that the standard macro NDEBUG is defined when compiling - see How to use the library.

The name EPR is an acronym for Ecm and Pollard-Rho, since those are the two main algorithms this library uses for factoring.  It's also a play on [Einstein-Podolsky-Rosen](https://en.wikipedia.org/wiki/EPR_paradox) for fun ([for now](https://en.wikipedia.org/wiki/Shor%27s_algorithm)).

Thanks to [Ben Buhrow](https://github.com/bbuhrow/yafu) for his great microecm code that implements ECM, which this library uses in an optimized form (exploiting [Clockwork](https://github.com/hurchalla/modular_arithmetic)) and extended to 128 bits.

## Design goals

The goal for EPR was to create a correct and easy to use library with extremely fast factoring (and primality checking) for native C++ integer types (i.e. 16/32/64 bit signed and unsigned ints).  Though it was optimized for native C++ integer types, it also efficiently supports 128 bit types - including the compiler extensions __int128_t and __uint128_t.  The EPR library uses variants of ECM and Pollard-Rho and deterministic Miller-Rabin algorithms, described [below](#algorithms), and has proofs for the new variants (with the exception of ECM) and extensive unit tests.  A portion of the optimizations (such as Montgomery arithmetic) that EPR uses are provided by the [Clockwork](https://github.com/hurchalla/modular_arithmetic) modular arithmetic library.

## Requirements

The EPR library requires compiler support for C++17 (if you are not using CMake, you may need to specify the option *-std="c++17"* when compiling).  Compilers that are confirmed to build the library without warnings or errors on x86 include clang6, clang10, gcc7, gcc10, intel compiler 19, and Microsoft Visual C++ 2017 and 2019.  The library is intended for use on all architectures (e.g. x86/64, ARM, RISC-V, Power), but has been tested only on x86/x64.  

For good performance you absolutely *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.

## Status

Released. All planned functionality and unit tests are finished and working correctly.

## Authors

* **Jeffrey Hurchalla**
* (microecm.h by **Ben Buhrow** and **Jeff Hurchalla**)

## License

This project is licensed under the MPL 2.0 License - see the [LICENSE.TXT](LICENSE.TXT) file for details

<br/>

## How to use the library

### With CMake

If you're using CMake for your project and you wish to add this library to it, then clone this git repository onto your system.  In your project's CMakeLists.txt file, add the following two lines with appropriate changes to their italic portions to match your project and paths ( an easy replacement for *your_binary_dir* is ${CMAKE_CURRENT_BINARY_DIR} ):  
add_subdirectory(*path_of_the_cloned_factoring_repository* &nbsp; *your_binary_dir*/factoring)  
target_link_libraries(*your_project_target_name* &nbsp; hurchalla_factoring)  

For good performance you absolutely *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.  You can do this by calling CMake with -DCMAKE_BUILD_TYPE=Release.  

It may help to see a simple [example project with CMake](examples/example_with_cmake).

### Without CMake

If you're not using CMake for your project, you'll need to install/copy the EPR library's headers and dependencies to some directory in order to use them.  To do this, first clone this git repository onto your system.  You'll need CMake on your system (at least temporarily), so install CMake if you don't have it.  Then from your shell run the following commands:  

>cd *path_of_the_cloned_factoring_repository*  
>mkdir tmp  
>cd tmp  
>cmake -S.. -B.  
>cmake --install . --prefix *the_folder_you_want_to_install_to*  
If you prefer, for the last command you could instead use CMake's default install location (on linux this is /usr/local) by omitting the --prefix and subsequent folder.  

This will copy all the header files needed for the factoring library to an "include" subfolder in the installation folder of your choosing.
When compiling your project, you'll of course need to ensure that you have that include subfolder as part of your include path.  

For good performance you absolutely *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.  You can generally do this by adding the option flag -DNDEBUG to your compile command.  

It may help to see a simple [example](examples/example_without_cmake).

## The API

The API consists of five header files in total (all the headers that are not under the *detail* folder).  These files are the three general purpose files *factorize.h*, *is_prime.h*, *greatest_common_divisor.h*, and in the resource_intensive_api folder, the two special purpose files *factorize_intensive_uint32.h* and *IsPrimeIntensive.h*.  Please view these files for their documentation.  A quick summary of the functions is provided below; in all cases T is a template parameter of integral type.  

*hurchalla::factorize(T x, int& num_factors)*.  Returns a std::array containing the factors of x.  
*hurchalla::factorize(T x, std::vector& factors)*.  Clears a std::vector and fills it with the factors of x.  
*hurchalla::greatest_common_divisor(T a, T b)*.  Returns the greatest common divisor of a and b.  
*hurchalla::is_prime(T x)*.  Returns true if x is prime.  Otherwise returns false.  

(from the resource_intensive_api folder)  
*hurchalla::IsPrimeIntensive(T x)*.  This is a functor that returns true if x is prime, and otherwise returns false.  Depending on the type T, this functor can use a very large amount of memory and can take many seconds to construct.  See IsPrimeIntensive.h for details.  
*hurchalla::factorize_intensive_uint32(uint32_t x, int& num_factors, const IsPrimeIntensive&lt;uint32_t,true&gt;& ipi)*.  Returns a std::array containing the factors of x.  Note that the IsPrimeIntensive argument will usually take many seconds to construct and will use a large amount of memory.  See IsPrimeIntensive.h for details.  

## Algorithms

For factoring: ECM (["Factoring Integers with Elliptic Curves"](https://www.jstor.org/stable/1971363)) by H.W. Lenstra Jr. (see also https://en.wikipedia.org/wiki/Lenstra_elliptic-curve_factorization).  For smaller integers: Pollard-Rho Brent (["An Improved Monte Carlo Factorization Algorithm"](https://maths-people.anu.edu.au/~brent/pub/pub051.html) by Richard Brent (see also https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants)).
Both ECM and Pollard-Rho preceded by a stage of trial division using special algorithms for divisibility (Section 10-17 from Hacker's Delight 2nd edition by Henry Warren, and "ALGORITHM A: IS_DIV_A" from ["Efficient long division via Montgomery multiply"](https://arxiv.org/abs/1303.0328) by Ernst W. Mayer).  

For primality testing: Deterministic Miller-Rabin (https://miller-rabin.appspot.com/ and https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants).  If using the resource_intensive_api, also Sieve of Eratosthenes for 32 bit and smaller types (https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes).

Special attention was paid to instruction level parallelism, to take advantage of typical pipelined/superscalar CPUs.  This resulted in modifications to some algorithms.  For example, Miller-Rabin is changed (along with its associated modular exponentiation) to perform more than one trial/exponentiation at a time.  And there is a changed Pollard-Rho Brent algorithm that simultaneously advances two independent sequences (from which it extracts factors), rather than just one.

Near-optimal hash tables are used for fast deterministic Miller-Rabin primality testing.  For information on the tables and how they were generated, you can view the [README.md](https://github.com/hurchalla/factoring/blob/master/include/hurchalla/factoring/detail/miller_rabin_bases/README.TXT), and see the [header files with the tables](include/hurchalla/factoring/detail/miller_rabin_bases).  These hash tables are the densest known at this time.  The general purpose functions *hurchalla::is_prime* and *hurchalla::factorize* use some of the smallest of these hash tables (8 to 320 byte), to minimize memory footprint and cache impact while still receiving a performance boost.  The resource_intensive_api functions use the largest of the hash tables for best possible performance.

The particular ECM variant used by this library originated from Ben Buhrow's microecm in [yafu](https://github.com/bbuhrow/yafu).

The Pollard-Rho-Brent algorithm uses an easy extra step that seems to be unmentioned in the literature.  The step is a "pre-loop" that advances as quickly as possible through a portion of the initial pseduo-random sequence before beginning the otherwise normal Pollard-Rho Brent algorithm.  The rationale for this is that every Pollard-Rho pseudo-random sequence begins with a non-periodic segment, and trying to extract factors from that segment is mostly wasted work since the algorithm logic relies on a periodic sequence.  Using a "pre-loop" that does nothing except iterate for a set number of times through the sequence thus improves performance on average, since it quickly gets past some of that unwanted non-periodic segment.  In particular, the "pre-loop" avoids calling the greatest common divisor, which would rarely find a factor during the non-periodic segment.  This optimization would likely help any form/variant of the basic Pollard-Rho algorithm.

## Performance Notes
If you're interested in experimenting, predefining certain macros when compiling can improve performance - see [macros_for_performance.md](macros_for_performance.md).
