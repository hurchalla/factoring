# Factoring

This is a factoring and primality testing library for C++17, optimized for 32 or 64 bit integer types and supporting integer types up to 128 bits.  It is a header-only library, designed for correctness and ideal performance.  Though 128 bit integers are supported, you should note that other algorithms would be more suitable for values much above (1<<64).  The algorithms used in this library are described further below.

This library requires a compiler that supports C++17 (if you are not using CMake, you may need to specify the option *-std="c++17"* when compiling).  Compilers that are confirmed to build the library without warnings or errors on x86 include clang6, clang10, gcc7, gcc10, intel compiler 19, and Microsoft Visual C++ 2017 and 2019.  The library is intended for use on all architectures (e.g. x86/64, ARM, RISC-V, Power), but has been tested only on x86/x64.  

For good performance you *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.

## Author

* **Jeffrey Hurchalla**

## How to use the library

### With CMake

If you're using CMake for your project and you wish to add this factoring library to it, then clone this git repository onto your system.  In your project's CMakeLists.txt file, add the following two lines with appropriate changes to their italic portions to match your project and paths ( an easy replacement for *your_binary_dir* is ${CMAKE_CURRENT_BINARY_DIR} ):  
add_subdirectory(*path_to_the_cloned_factoring_repository* &nbsp; *your_binary_dir*/factoring)  
target_link_libraries(*your_project_target_name* &nbsp; hurchalla_factoring)  

For best performance you *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.  You can do this by calling CMake with -DCMAKE_BUILD_TYPE=Release.  

It may help to see a simple [example project with CMake](example_with_cmake).

### Without CMake

If you're not using CMake for your project, you'll need to install/copy these factoring headers and dependencies to some directory in order to use them.  To do this, first clone this git repository onto your system.  You'll need CMake on your system (at least temporarily), so install CMake if you don't have it.  Then from your shell run the following commands:  

>cd *path_to_the_cloned_factoring_repository*  
>mkdir tmp  
>cd tmp  
>cmake -S.. -B.  
>cmake --install . --prefix *the_folder_you_want_to_install_to*  
If you prefer, for the last command you could instead use CMake's default install location (on linux this is /usr/local) by omitting the --prefix and subsequent folder.  

This will copy all the header files needed for the factoring library to an "include" subfolder in the installation folder of your choosing.
When compiling your project, you'll of course need to ensure that you have that include subfolder as part of your include path.  

For good performance you *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.  You can generally do this by adding the option flag -DNDEBUG to your compile command.  

It may help to see a simple [example](example_without_cmake).

## The API

The API consists of five header files in total (all the headers that are not under the *detail* folder).  These files are the three general purpose files *factorize.h*, *is_prime.h*, *greatest_common_divisor.h*, and in the resource_intensive_api folder, the two special purpose files *factorize_intensive_uint32.h* and *IsPrimeIntensive.h*.  Please view these files for their documentation.  A quick summary of the functions is provided below; in all cases T is a template parameter of unsigned integral type.  

*hurchalla::factorize(T x, int& num_factors)*.  Returns a std::array containing the factors of x.  
*hurchalla::factorize_to_vector(T x)*.  Returns a std::vector containing the factors of x.  
*hurchalla::greatest_common_divisor(T a, T b)*.  Returns the greatest common divisor of a and b.  
*hurchalla::is_prime(T x)*.  Returns true if x is prime.  Otherwise returns false.  

(from the resource_intensive_api folder)  
*hurchalla::IsPrimeIntensive(T x)*.  This is a functor that returns true if x is prime, and otherwise returns false.  Depending on the type T, this functor can use a very large amount of memory and can take many seconds to construct.  See IsPrimeIntensive.h for details.  
*hurchalla::factorize_intensive_uint32(uint32_t x, int& num_factors, const IsPrimeIntensive&lt;uint32_t,true&gt;& ipi)*.  Returns a std::array containing the factors of x.  Note that the IsPrimeIntensive argument will usually take many seconds to construct and will use a large amount of memory.  See IsPrimeIntensive.h for details.  

## Algorithms

For factoring: Pollard-Rho Brent (["An Improved Monte Carlo Factorization Algorithm"](https://maths-people.anu.edu.au/~brent/pub/pub051.html) by Richard Brent, and https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants),
preceded by a stage of trial division using special algorithms for divisibility (Section 10-17 from Hacker's Delight 2nd edition by Henry Warren, and "ALGORITHM A: IS_DIV_A" from ["Efficient long division via Montgomery multiply"](https://arxiv.org/abs/1303.0328) by Ernst W. Mayer).  

For primality testing: Deterministic Miller-Rabin (https://miller-rabin.appspot.com/ and https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants).  If using the resource_intensive_api, also Sieve of Eratosthenes for 32 bit and smaller types (https://en.wikipedia.org/wiki/Sieve_of_Eratosthenes).

Special attention was paid to instruction level parallelism, to take advantage of typical pipelined/superscalar CPUs.  This resulted in modifications to some algorithms, particularly Miller-Rabin (and its associated modular exponentiation) to perform more than one trial/exponentiation at a time.

Near-optimal hash tables are used for fast deterministic Miller-Rabin primality testing.  For information on the tables and how they were generated, you can view the [README.md](https://github.com/hurchalla/factoring/blob/master/include/hurchalla/factoring/detail/miller_rabin_bases/README.TXT), and see the [header files with the tables](include/hurchalla/factoring/detail/miller_rabin_bases).  The general purpose functions *hurchalla::is_prime* and *hurchalla::factoring* use some of the smallest of these hash tables (8 to 320 byte), to minimize memory footprint and cache impact while still receiving a performance boost.  The resource_intensive_api functions use the largest of the hash tables for best possible performance.

The Pollard-Rho-Brent algorithm uses an easy extra step that seems to be unmentioned in the literature (though the idea is obvious enough that it's unlikely to be new).  The step is a "pre-loop" that advances as quickly as possible through a portion of the initial pseduo-random sequence before beginning the otherwise normal Pollard-Rho Brent algorithm.  The rationale for this is that every Pollard-Rho pseudo-random sequence begins with a non-periodic segment, and trying to extract factors from that segment is mostly wasted work since the algorithm logic relies on a periodic sequence.  Using a "pre-loop" that does nothing except iterate for a set number of times through the sequence thus improves performance on average, since it quickly gets past some of that unwanted non-periodic segment.  This optimization would likely help any form/variant of the basic Pollard-Rho algorithm.

## Status
Released, version 1.0.

## Miscellaneous
If you're interested in experimenting, predefining certain macros when compiling can improve performance - see [macros_for_performance.md](macros_for_performance.md).
