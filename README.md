# factoring

TODO
----
Define the results for the factorization functions when x < 2.

Document:
   HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE
   HURCHALLA_ENABLE_EXPERIMENTAL_TRIAL_DIVIDE_UINT32_VIA_DOUBLE

   HURCHALLA_POLLARD_RHO_NEVER_USE_MONTGOMERY_MATH
   HURCHALLA_POLLARD_RHO_TRIAL_GCD_THRESHOLD
   HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR
   HURCHALLA_USE_WHEEL_FACTORIZATION
       HURCHALLA_POLLARD_RHO_MAX_TRIAL_FACTOR
   HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME   PollardRhoTrial

   HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048
   HURCHALLA_TEST_SMALL_TRIAL_DIVISION2048_INDEX_LIMIT

   HURCHALLA_MILLER_RABIN_ALLOW_EVEN_NUMBERS

   HURCHALLA_CHOOSE_ALTERNATE_MILLER_RABIN_BASES64_3


Find best value for HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR:
The best value for HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR
depends on whether the number (x) to test is very large or small.
If it's very large, for example close to UINT64_MAX, then we want
to use a somewhat large value for the macro.  If x is small, we
want to use a small value for the macro.
This is because is_prime_miller_rabin becomes slower for larger
x since it needs more witnesses.  Numbers close to uint128_max are
incredibly bad for miller rabin at the moment due to their insane
witness set - however I'm not trying to make program-wide changes
just to help optimization for uint128_t, but I'm willing to make
simple isolated changes that help it.

At the moment it appears to be fastest to set
HURCHALLA_ISPRIME_MAX_TRIAL_FACTOR to 16, essentially skipping
the is_prime_wheel.
Perhaps a replacement of the wheel with an
is_prime_small_divisions256 would make sense?
Possibly easiest way to do this is to change small_trial_division256()
to add a bool template param StopAtFirstFactor, and put the impl of 
small_trial_division256 inside a class (and use static member functions),
so as to allow partial specialization.
Or simpler, I could just use a function similar to the generic version in
smalL_trial_division256.


make miller_rabin more efficient overall, using tricks like hashing and ILP,
and using MF's pow with ILP.

Add a (potential) factorization step to miller rabin.
Add check for 1 in inner loop of miller-rabin (I suspect this may be the same TODO as the last item).


MontgomeryForm - force_inline on every function?

NDEBUG defined for all release builds?

Compare performance of asm to non-asm functions via the macros in modular arithmetic README

publish my generalized form of the Dumas algorithm

maybe blog about gcc uint128_t bug- see end of test_REDC.cpp.
