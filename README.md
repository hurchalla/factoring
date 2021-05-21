# factoring

TODO
----
Define the results for the factorization functions when x < 2.

Assuming I can remove any need factorize_wheel210 might have for it,
get rid of the "next_prime" param in factorize_trialdivision().
Likewise get rid of it in is_prime_trialdivision() if possible.

Get rid of pIterationsPerformed in PollardRhoTrial functors

Use a faster gcd
Test whether the pollard-rho trials' use of a gcd functor is as fast as a direct call of gcd

Can I completely remove factorize_wheel210.h, or place it inside test or experimental?

Write tests for the headers that don't yet have tests


Document:
   -- trial_divide_mayer.h --
   HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE
   HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE

   -- impl_is_prime.h --
   HURCHALLA_ISPRIME_TRIALDIV_SIZE

   -- ImplIsPrimeIntensive.h --
   HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE
   HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_TYPE  (can be PrimeTrialDivisionMayer or PrimeTrialDivisionWarren)

   -- impl_factorize.h --
   HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR
   HURCHALLA_PR_TRIAL_DIVISION_TEMPLATE
   HURCHALLA_USE_PR_TRIAL_DIVISION
   HURCHALLA_PR_TRIAL_DIVISION_SIZE
   HURCHALLA_PR_TRIAL_DIVISION_INDEX_LIMIT

   -- pollard_rho_factorize.h --
   HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME   PollardRhoTrial or PollardRhoBrentTrial (or one of the trials in experimental)
   HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH

   -- PollardRhoTrial.h --
   HURCHALLA_POLLARD_RHO_GCD_THRESHOLD

   -- PollardRhoBrentTrial.h --
   HURCHALLA_POLLARD_RHO_BRENT_GCD_THRESHOLD
   HURCHALLA_POLLARD_RHO_BRENT_PRE_CYCLE_SIZE
   HURCHALLA_POLLARD_RHO_BRENT_INITIAL_CYCLE_SIZE

   -- factorize_wheel210.h --
   HURCHALLA_WHEELFACTOR_TRIAL_DIVISION_TEMPLATE

   -- is_prime_miller_rabin.h --
   HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS
   HURCHALLA_DEFAULT_TO_UNHASHED_MILLER_RABIN
   HURCHALLA_CHOOSE_ALTERNATE_MILLER_RABIN_BASES64_3


Find best value for HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE.

Find best value for HURCHALLA_ISPRIME_TRIALDIV_SIZE:
The best value for HURCHALLA_ISPRIME_TRIALDIV_SIZE
depends on whether the number (x) to test is very large or small.
If it's very large, for example close to UINT64_MAX, then we want
to use a somewhat large value for the macro.  If x is small, we
want to use a somewhat small value for the macro.
This is because is_prime_miller_rabin has worse performance for large
integer types when primality testing numbers that happen to be prime,
since it needs to test more bases than it would for a smaller type.
Testing uint128 numbers is extremely slow (relative to smaller types)
for miller rabin at the moment due to the insane number of bases it
uses for 128 bit primality testing.  However I'm not trying to make
program-wide changes just to help optimization for uint128_t, but I'm
willing to make simple isolated changes that help it.
Brief performance tests I did on Haswell led me to use a default of
HURCHALLA_ISPRIME_TRIALDIV_SIZE set to 16.

Add a (potential) factorization step to miller rabin.
Add check for 1 in inner loop of miller-rabin (I suspect this may be
the same TODO as the last item).


MontgomeryForm - force_inline on every function?

NDEBUG defined for all release builds?

Compare performance of asm to non-asm functions via the macros in
modular arithmetic README

publish my generalized form of the Dumas algorithm

maybe blog about gcc uint128_t bug- see end of test_REDC.cpp.

A possible perf improvement for impl_is_prime might be to scale the size of the
trial division at run time depending on how large x is (using the optional third
argument of is_prime_trialdivision).  The ideal scaling might be y = a * sqrt(x)
for some constant a, but sqrt would take too long.  A linear scaling might be
good enough, perhaps fit to the largest 80-90% of values of x.  It feels like
it's more trouble than it's worth.
The same thing is possible for ImplIsPrimeIntensive.h
-- Probably a better idea is to remove the size_limit parameter from
is_prime_trialdivision() since at least right now it's never used - this would
let me simplify the function, which would have the side benefit of making it
slightly faster too.

Maybe simplify is_prime_miller_rabin.h by moving stuff that's not really used
into an experimental folder
