# factoring

TODO
----
Define the results for the factorization functions when x < 2.

Assuming I can remove any need factorize_wheel210 might have for it,
get rid of the "next_prime" param in factorize_trialdivision().
Likewise get rid of it in is_prime_trialdivision() if possible.


Document:
   HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE

   HURCHALLA_ISPRIME_TRIALDIV_SIZE
   HURCHALLA_ISPRIME_INTENSIVE_TRIALDIV_SIZE

   -- factorize.h --
   HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME   PollardRhoTrial
   -- impl_factorize.h --
   HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR
   HURCHALLA_TRIAL_DIVISION_TEMPLATE
   HURCHALLA_USE_PR_TRIAL_DIVISION
   HURCHALLA_PR_TRIAL_DIVISION_SIZE
   HURCHALLA_PR_TRIAL_DIVISION_INDEX_LIMIT
   -- pollard_rho_factorize.h --
   HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH
   -- PollardRhoTrial.h --
   HURCHALLA_POLLARD_RHO_TRIAL_GCD_THRESHOLD

   -- factorize_wheel210.h --
   HURCHALLA_WHEELFACTOR_TRIAL_DIVISION_TEMPLATE

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
