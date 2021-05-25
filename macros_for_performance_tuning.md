
Optional macros to predefine to tune performance
------------------------------------------------
There are a number of macros you can optionally predefine to tune the
performance on your system for the factoring and primality testing functions.
You would predefine one or more of these macros when compiling the sources.  For
example, with clang you would compile with the -D compilation flag like this:
clang++ -DHURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD=300 
\
\
For tuning of factorize() and factorize_intensive32():

HURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD - this is the number of iterations of
the Pollard-Rho Brent algorithm that will be performed between calling the
greatest common divisor to try to extract a factor.  400 is the default
currently.  You can predefine this macro to a number that works best on your
system through experimentation.

HURCHALLA_POLLARD_RHO_BRENT_STARTING_LENGTH - upon starting Pollard-Rho Brent
factoring, this is the number of iterations that the 'hare' (of tortoise and
hare: https://en.wikipedia.org/wiki/Cycle_detection#Floyd.27s_Tortoise_and_Hare)
initially moves through the Pollard-Rho Brent pseudo-random sequence while the
'tortoise' stands still.  Once the hare reaches the end of this initial
distance, the distance doubles and the tortoise is set to the position of the
hare, and the process repeats again using the doubled distance, then repeats
again, over and over.  The default for this starting length macro is currently
20.  You can predefine this macro to a different number to set a different
initial length.  Every time the distance is doubled the greatest common divisor
(GCD) gets called, so this initial distance affects how often the GCD is called
early on in the process, thus impacting performance.  Note that the GCD is also
called after every HURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD iterations of the
sequence.

HURCHALLA_POLLARD_RHO_BRENT_PRE_LENGTH - this is the number of iterations of the
pseudo-random sequence to evaluate before beginning the normal Pollard-Rho Brent
factoring algorithm.  40 is the default currently.  You can predefine this macro
to a number that works best on your system.  Some background: the expected
length of the starting nonperiodic segment of the Pollard-Rho sequence is the
same as the expected length of the sequence's period, and in a sense any work we
do during the nonperiodic segment is wasted - ideally we just want to get
through the nonperiodic segment as quickly as possible.  Since each iteration
during the pre_length (as given by this macro) does nothing but advance the
sequence, it completes quickly compared to the Pollard-Rho Brent algorithm
(which does other work).  Thus we generally save time by specifying a pre_length
that we use to quickly pass over a portion of the nonperiodic segment, prior to
starting Pollard-Rho Brent factoring.

HURCHALLA_PR_TRIAL_DIVISION_SIZE - this is the number of small primes (starting
at 2,3,5,7, etc) that will be trialed as potential factors via trial division
before trying to find factors with Pollard Rho.  135 is the default currently.
You can predefine this macro to another number that might provide better
performance on your system.

HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE - this is the name of the template that
performs the trial division; if you predefine it, you must set it to either
PrimeTrialDivisionWarren or PrimeTrialDivisionMayer.  PrimeTrialDivisionWarren
is the default.  If you have a CPU with very fast division instructions
(some CPUs from 2019 or later), you might be able to improve performance by
predefining this macro to PrimeTrialDivisionMayer while also predefining the
macro HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.  If you need to minimize resource
usage, then predefine this macro to PrimeTrialDivisionMayer, since the Mayer
template uses about a fifth of the memory of PrimeTrialDivisionWarren.  If you
do predefine this macro, you will very likely also want to predefine
HURCHALLA_PR_TRIAL_DIVISION_SIZE to a new value that works well with this macro
choice.

HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH - if predefined, this macro
causes the factorization functions to use standard division for modular
arithmetic, instead of montgomery arithmetic.  By default this macro is not
defined.  If you have a CPU with very fast division instructions, you might be
able to improve performance by predefining this macro while also predefining the
macro HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.
\
\
For tuning of is_prime():

HURCHALLA_ISPRIME_TRIALDIV_SIZE - this is the number of small primes (starting
at 2,3,5,7, etc) that will be trialed as potential factors via trial division
before testing primality with Miller-Rabin.  15 is the current default.  You can
predefine this macro to a different number that is better tuned for your system.
\
\
For tuning of IsPrimeIntensive:

HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE - this is the number of small primes
(starting at 2,3,5,7, etc) that will be trialed as potential factors via trial
division before testing primality with Miller-Rabin.  75 is the current default.
You can predefine this macro to a different number that is better tuned for your
system.  Note that for IsPrimeIntensive, when testing primality of uint32_t and
smaller types, this macro has no effect because IsPrimeIntensive tests primality
of those types solely using the Sieve of Eratosthenes.

HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_TYPE - this is the name of the template
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

HURCHALLA_PREFER_EUCLIDEAN_GCD - predefining this macro allows you to use the
Euclidean algorithm for greatest common divisor (GCD), rather tham the default
Binary/Stein algorithm.  However, you will only get the Euclidean GCD if you
also predefine HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE.  Note that even with both
macros predefined, the Binary GCD will still be used for any integer type T that
is larger than the CPU's native bit width.

HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME - the name of the algorithm/functor
that carries out the Pollard Rho factoring trial.  Normally predefining this
will not help performance, but if you wish you can predefine it to either
PollardRhoTrial or PollardRhoBrentTrial.  The default is PollardRhoBrentTrial.
PollardRhoTrial will very likely be slower, but you can experiment.  If you
choose PollardRhoTrial, you should also try to predefine
HURCHALLA_POLLARD_RHO_GCD_THRESHOLD to an optimal number (see
HURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD for related details) for your system.
