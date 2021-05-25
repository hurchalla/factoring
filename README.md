# factoring

You can tune the performance of these functions on your system by predefining macros: see [macros_for_performance_tuning.md](macros_for_performance_tuning.md).

TODO
----

do pathological numbers exist that fail most Pollard-Rho factoring trials (i.e. failing regardless of 'c' in the sequence x[i+1] = (x[i] * x[i]) + c)?

evaluate: a step to factor powers might help.
evaluate: a (potential) factorization step in miller rabin might help, but probably not.

MontgomeryForm - force_inline on every function?
MontgomeryFormCommon - show proof that fmadd and fmsub are correct.
NDEBUG defined for all release builds?
Compare performance of asm to non-asm functions via the macros in modular arithmetic README

publish the generalized form of the Dumas integer inverse algorithm
