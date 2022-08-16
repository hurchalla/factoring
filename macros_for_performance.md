
Optional macros to predefine to tune performance
------------------------------------------------
There are a number of macros you can optionally predefine to tune the
performance on your system for the factoring and primality testing functions.
You would predefine one or more of these macros when compiling the sources.  For
example, with CMake you would add the command "target_compile_definitions" to
the CMakeLists.txt, similarly to the following:
target_compile_features(hurchalla_factoring INTERFACE HURCHALLA_FACTORING_ECM_THRESHOLD_BITS=40)

If you are not using CMake, then with clang or gcc you would compile with the -D
compilation flag.  For example: 
clang++ -DHURCHALLA_FACTORING_ECM_THRESHOLD_BITS=40 
\
\
Macros for factorize() and factorize_intensive32():

HURCHALLA_FACTORING_EXPECT_LARGE_FACTORS - when this macro is defined, the
factoring functions expect only very large prime factors.  Therefore they
skip the normal trial division phase, as well as any other phase
designed to find small factors.  Defining this macro also slightly
lowers the size threshold for the factoring functions to choose the ECM
algorithm over Pollard-Rho.

HURCHALLA_FACTORING_ECM_THRESHOLD_BITS - this macro provides the cutoff between
using Pollard-Rho or ECM for factoring.  Any number to be factored that is
(after trial division) greater than (1 << HURCHALLA_FACTORING_ECM_THRESHOLD_BITS)
will be factored using ECM, and any number that is smaller will be factored using
Pollard-Rho.  The default value is 34 when HURCHALLA_FACTORING_EXPECT_LARGE_FACTORS
is defined, and 40 when HURCHALLA_FACTORING_EXPECT_LARGE_FACTORS is not defined.

HURCHALLA_TRIAL_DIVISION_SIZE_SMALL and HURCHALLA_TRIAL_DIVISION_SIZE_LARGE -
you can predefine these macros to integer values that together provide a range.
The defaults are 109 and 139, respectively, and the 'small' macro value must be
less than the 'large'.  This range influences the number of
small primes (starting at 2,3,5, etc) that will be trialed as potential factors
during the initial trial division stage.  When factoring a large number, a higher
value from this range will be used, and when factoring a small number, a lower
value will be used.  If you predefine one of these macros, you should also
predefine the other.

HURCHALLA_TRIAL_DIVISION_TEMPLATE - this is the name of the template that
performs the trial division; if you predefine it, you must set it to either
PrimeTrialDivisionWarren or PrimeTrialDivisionMayer.  PrimeTrialDivisionWarren
is the default.  If you have a CPU with very fast division instructions
(some CPUs from 2019 or later), you might be able to improve performance by
predefining this macro to PrimeTrialDivisionMayer while also predefining the
macro HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.  If you need to minimize resource
usage, then predefine this macro to PrimeTrialDivisionMayer, since the Mayer
template uses about a fifth of the memory of PrimeTrialDivisionWarren.  If you
do predefine this macro, you will very likely also want to predefine
HURCHALLA_TRIAL_DIVISION_SIZE_SMALL and HURCHALLA_TRIAL_DIVISION_SIZE_LARGE to
new values that work best with this macro choice.

HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME - the name of the algorithm/functor
that carries out a Pollard Rho factoring trial.  You will likely find that
either PollardRhoBrentSwitchingTrial or PollardRhoBrentTrialParallel is fastest
(the default is PollardRhoBrentSwitchingTrial).  See the experimental trial
folder's [README_pollard_rho.md](include/hurchalla/factoring/detail/experimental/README_pollard_rho.md)
for details.
\
\
Macros for is_prime():

HURCHALLA_ISPRIME_TRIALDIV_SIZE - the number of small primes (starting
at 2,3,5,7, etc) that will be trialed as potential factors via trial division
before testing primality with Miller-Rabin.  21 is the default.  You can
predefine this macro to a different number that is better tuned for your system.
\
\
Macros for IsPrimeIntensive:

HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE - the number of small primes
(starting at 2,3,5,7, etc) that will be trialed as potential factors via trial
division before testing primality with Miller-Rabin.  75 is the default.
You can predefine this macro to a different number that is better tuned for your
system.  Note that for IsPrimeIntensive, when testing primality of uint32_t and
smaller types, this macro has no effect because IsPrimeIntensive tests primality
of those types solely using the Sieve of Eratosthenes.

HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_TYPE - the name of the template
that performs the trial division; if you predefine it, you must set it to either
PrimeTrialDivisionWarren or PrimeTrialDivisionMayer.  PrimeTrialDivisionWarren
is the default.  If you have a CPU with very fast division instructions (some
CPUs from 2019 or later), you might be able to improve performance by
predefining this macro to PrimeTrialDivisionMayer while also predefining the
macro HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.  If you need to minimize resource
usage, then predefine this macro to PrimeTrialDivisionMayer, since the Mayer
template uses about a fifth of the memory of PrimeTrialDivisionWarren.  If you
do predefine this macro, you will very likely also want to predefine
HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE to a new value that works well with
this macro choice.
\
\
Miscellaneous macros:

HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE - predefine this macro if your CPU has very
fast division instructions (usually this is present only in CPUs from around
2019 or later).

HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE - you should usually predefine this macro if
your microprocessor lacks a division instruction.

HURCHALLA_FACTORIZE_NEVER_USE_MONTGOMERY_MATH - if predefined, this macro
causes the factorization functions to use standard division for modular
arithmetic, instead of montgomery arithmetic.  By default this macro is not
defined.  If you have a CPU with very fast division instructions, you might be
able to improve performance by predefining this macro while also predefining the
macro HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.

HURCHALLA_PRB_GCD_THRESHOLD - the maximum number of iterations of
the Pollard-Rho-Brent algorithm that will be allowed before calling the
greatest common divisor to try to extract a factor.  The current default is
around 600.

HURCHALLA_PRB_STARTING_LENGTH - the number of iterations of the
Pollard-Rho-Brent pseudo-random sequence to evaluate before beginning the normal
factoring algorithm.  19 is the default.  Some background: the expected
length of the starting nonperiodic segment of the Pollard-Rho sequence is the
same as the expected length of the sequence's periodic segment, and in a sense
any work we do during the nonperiodic segment is pointless.  Ideally we would
just like to move past the nonperiodic segment as quickly as possible.  So the
starting length (as given by this macro) lets us advance through a portion of
the nonperiodic segment while doing minimal work.  This increases performance,
so long as the 'starting_length' is not too much greater than the actual
nonperiodic segment's length.  Note: if you expect to have only large factors,
you will generally be able to improve the performance of Pollard-Rho-Brent
factoring by predefining this macro to a larger value than the default, because
both the periodic and nonperiodic segments will tend to be large when you have
large factors.

HURCHALLA_PREFER_EUCLIDEAN_GCD - predefining this macro allows you to use the
Euclidean algorithm for greatest common divisor (GCD), rather tham the default
Binary/Stein algorithm.  However, you will only get the Euclidean GCD if you
also predefine HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.  Note that even with both
macros predefined, the Binary GCD will still be used for any integer type T that
is larger than the CPU's native bit width.

HURCHALLA_ALLOW_INLINE_ASM_ALL - predefining this macro will enable all
available inline asm functions.  In some cases it may improve performance up to
25%, and in other cases it may make essentially no difference or harm
performance.  In all cases, inline asm is very difficult to thoroughly test,
since the code surrounding the inline asm under test will determine part of the
machine code generated from the inline asm.  Generally speaking, it is
[difficult to recommend inline asm](https://gcc.gnu.org/wiki/DontUseInlineAsm)
unless there is a large performance benefit or performance is critical.
\
\
Miscellaneous macros from dependencies:

Since this library depends upon the hurchalla modular_arithmetic library, it is
also affected by that library's optional [performance macros](https://github.com/hurchalla/modular_arithmetic/blob/master/macros_for_performance.md).

