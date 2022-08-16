
The PollardRho\*.h header files in this folder are experimental functors for Pollard-Rho trials.  
To use an experimental functor, predefine the macro HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME and give it the name of the experimental functor.  For example, if you want PollardRhoTrial and you are compiling with clang, you could invoke clang as follows:  
clang++ -DHURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME=PollardRhoTrial  ...more options and files...  

valid names to give HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME are:  
PollardRhoTrial  
PollardRhoBrentTrial  
PollardRhoBrentTrialParallel  
PollardRhoBrentSwitchingTrial  (currently this is the default, so it's not necessary to predefine the macro to this functor)  

The PollardRhoBrentSwitchingTrial functor or the PollardRhoBrentTrialParallel functor will likely perform best on your system, but you can try others.

The macro is used only in the file ../FactorizeStage2.h.

Functor Descriptions
--------------------

PollardRhoTrial.h:  
The original Pollard-Rho algorithm is by JM Pollard, as described at https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm and in his paper "A Monte Carlo method for factorization".  The implementation of it here uses a well-known improvement to the paper's algorithm, increasing the speed by using a large number of iterations before calling the greatest common divisor (GCD).  This functor would very likely be slower on a given system than Pollard-Rho-Brent, described next.

../PollardRhoBrentTrial.h:  
Pollard-Rho-Brent is the most well known variant of the Pollard-Rho algorithm and it is present in the mainline /detail folder (not in this experimental folder).  You can read about Pollard-Rho-Brent at https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants and in Richard Brent's paper ["An Improved Monte Carlo Factorization Algorithm"](https://maths-people.anu.edu.au/~brent/pub/pub051.html).

../PollardRhoBrentTrialParallel.h
A variant on Pollard-Rho-Brent that executes two different and independent Pollard-Rho-Brent trials in its main loop (note this functor is single-threaded, intentionally).  The two trials should execute at the same time due to modern processors' capabilities for instruction level parallelism via pipelining and/or superscalar execution units.  Running the two trials in parallel could in theory improve the average performance by up to 1.41x compared to a normal single trial of Pollard-Rho-Brent, but in practice due to extra overhead and imperfect instruction level parallelism, the improvement (if any) is much less.  For example, on my test system it has benchmarked to have essentially the same performance as the plain PollardRhoBrentTrial.  Nevertheless, it has the potential to be the fastest of all the functors when run on a recent processor that possesses more than one multiplication execution unit, since this functor provides the most opportunities for a CPU to execute multiple instructions in parallel.

../PollardRhoBrentSwitchingTrial.h:  
The default functor.  It is a novel variant on Pollard-Rho-Brent that interleaves two different trials - it always advances both trial sequences, but it switches back and forth between which of the two trials gets accumulated for a GCD.  The two parallel trials make it similar to PollardRhoBrentParallel, but it executes in its main loop exactly three Montgomery multiplications instead of four (PollardRhoBrentParallel).  These multiplications can in principle be executed at the same time (they do not depend upon one another), which maximizes instruction level parallelism on modern x86 processors - x86 processors typically can pipeline up to three independent multiplications.  Since it benchmarks noticeably fastest on my test system, it is the default.  However, we can reasonably expect that the trial PollardRhoBrentTrialParallel might benchmark faster on newer systems.

A Note on Montgomery's Variant
------------------------------
A further (unavailable) variant exists on the Pollard-Rho-Brent algorithm, by Peter Montgomery and described in section 3 (Brent's Improvement to Monte Carlo) of his paper ["Speeding the Pollard and Elliptic Curve Methods of Factorization"](https://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866113-7/S0025-5718-1987-0866113-7.pdf).  This variant uses slightly fewer modular multiplications, at the cost of extra complexity.  It is an interesting variation, in theory, to use when running two separate trials at the same time (like PollardRhoBrentTrialParallel.h).  The bottleneck when running two simultaneous trials of Pollard-Rho-Brent is typically the CPU's multiply execution unit, and so this variant's lower number of multiplications seems as though it might increase speed.  In practice, when I implemented it and a parallel version of it, it benchmarked noteably slower than the corresponding PollardRhoBrent variants described above, probably due to the overhead of the extra algorithm operations it requires.
