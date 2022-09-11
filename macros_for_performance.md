
Optional macros to predefine to tune performance
------------------------------------------------
There are a number of macros you can optionally predefine to tune the
performance on your system for the factoring and primality testing functions.
You would predefine one or more of these macros when compiling *your* sources,
given that this is a header-only library.

For example, if you are compiling using clang or gcc from the command line, you would
specify the -D compilation flag, similarly to the following: 
clang++ -DHURCHALLA_TRIAL_DIVISION_SIZE=139 ...more arguments...
As another example, if you are using CMake you would add the command "target_compile_definitions"
to your CMakeLists.txt, similarly to the following: 
target_compile_definitions(&lt;your_target_name&gt;  PRIVATE  HURCHALLA_TRIAL_DIVISION_SIZE=139) 
\
\
Macros for factorize() and factorize_intensive32():

HURCHALLA_FACTORING_ECM_THRESHOLD_BITS - this macro specifies the cutoff between
using Pollard-Rho or ECM for factoring.  Any number to be factored that is
(after trial division) greater than (1 << HURCHALLA_FACTORING_ECM_THRESHOLD_BITS)
will be factored using ECM, and any number that is smaller will be factored using
Pollard-Rho.  The default value is 40.

HURCHALLA_TRIAL_DIVISION_SIZE - this macro specifies the number of small primes
(starting at 2,3,5, etc) to trial as potential factors during the initial trial
division stage of factoring.  The default is 139.

HURCHALLA_TRIAL_DIVISION_TEMPLATE - this is the name of the template that
performs the trial division during factoring; if you predefine it, you must set it to either
PrimeTrialDivisionWarren or PrimeTrialDivisionMayer.  PrimeTrialDivisionWarren
is the default.  If you have a CPU with very fast division instructions
(some CPUs from 2019 or later), you might be able to improve performance by
predefining this macro to PrimeTrialDivisionMayer while also predefining the
macro HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.  If you need to minimize resource
usage, then predefine this macro to PrimeTrialDivisionMayer, since the Mayer
template uses about a fifth of the memory of PrimeTrialDivisionWarren.  If you
do predefine this macro, you will very likely also want to predefine
HURCHALLA_TRIAL_DIVISION_SIZE to a new value that works best with this macro choice.

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

HURCHALLA_FACTORING_DISALLOW_INLINE_ASM - predefining this macro will prevent
any inline asm from being compiled.  If this macro is not predefined, this
library by default will use inline asm for improved performance.  However,
inline asm is difficult to thoroughly test, because the code surrounding the
inline asm under test in part determines what binary instructions are generated
from the inline asm.  Thus inline asm is at much greater risk of bugs than
ordinary code.  It can often be wise to 
[avoid inline asm](https://gcc.gnu.org/wiki/DontUseInlineAsm)
if performance is not critical.
\
\
Miscellaneous macros from dependencies:

Since this library depends upon the hurchalla modular_arithmetic library, it is
also affected by that library's optional [performance macros](https://github.com/hurchalla/modular_arithmetic/blob/master/macros_for_performance.md).

