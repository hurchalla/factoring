
Each header file in this folder is an experimental functor for Pollard-Rho trials.  
To use an experimental functor, predefine the macro HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME and give it the name of the experimental functor.  For example, if you want PollardRhoTrial and you are compiling with clang, you could invoke clang as follows:  
clang++ -DHURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME=PollardRhoTrial  ...more options and files...  

valid names to give HURCHALLA_POLLARD_RHO_TRIAL_FUNCTOR_NAME are:  
PollardRhoTrial  
PollardRhoBrentTrial  (but since this is the default, it's not necessary to predefine the macro to this functor)  
PollardRhoBrentMontgomeryTrial  
PollardRhoBrentMontgomeryTrialParallel  

The default Pollard-Rho Brent functor will very likely perform best on your system, but you can try others.

The macro is used only by the file ../factorize_pollard_rho.h.

Functor Descriptions
--------------------

PollardRhoTrial.h:  
The original Pollard-Rho algorithm is by JM Pollard, as described at https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm and in his paper "A Monte Carlo method for factorization".  The trial functor within this experimental folder called PollardRhoTrial.h represents a slightly improved version of this algorithm, using a large number of iterations before calling the greatest common divisor (GCD) for improved efficiency.  This functor would very likely be slower on a given system than the default Pollard-Rho-Brent, described next.

(this next default functor is mainline code and not experimental, but it deserves description)  
../PollardRhoBrentTrial.h:  
Pollard-Rho-Brent is the most well known variant of the Pollard-Rho algorithm and it is the default and present and used in the mainline /detail folder (not in this experimental folder).  You can read about Pollard-Rho-Brent at https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants and in Brent's paper  "An Improved Monte Carlo Factorization Algorithm" at https://maths-people.anu.edu.au/~brent/pub/pub051.html.

PollardRhoBrentMontgomeryTrial.h:  
A variant on the Pollard-Rho algorithm is by Peter Montgomery, as described in section 3 (Brent's Improvement to Monte Carlo) of his paper "Speeding the Pollard and Elliptic Curve Methods of Factorization" at https://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866113-7/S0025-5718-1987-0866113-7.pdf  This variation uses slightly fewer modular multiplications, at the cost of extra complexity.  The trial functor representing this algorithm is PollardRhoBrentMontgomeryTrial.h.  In practice it benchmarked slower than Pollard-Rho Brent on my test system.

PollardRhoBrentMontgomeryTrialParallel.h:  
This is a second functor that also uses Montgomery's variant.  This functor performs two Pollard-Rho trials at the same time using Montgomery's variant, taking advantage of instruction level parallelism (note this functor is single-threaded, intentionally).  Running two Pollard-Rho trials in parallel in theory speeds up the answer by 1.41x (sqrt2) if we disregard costs of parallelism.  Unfortunately even with instruction level parallelism the almost doubling of instructions comes at a significant cost.  The end result on the tested CPU Haswell was close to break-even in performance to the simpler non-parallel PollardRhoBrentMontgomeryTrial.h, and as before, Pollard-Rho Brent performed better.  Montgomery's variant seems the most amenable of all the Pollard-Rho algorithms to running two trials while taking advantage of instruction level parallelism, since it has slightly fewer modular multiplications than the other algorithms.  Similar modifications to Pollard-Rho or Pollard-Rho-Brent thus do not seem promising to do better than break-even, and indeed a parallel version of Pollard-Rho-Brent that I quickly coded (but did not include here) provided no performance benefit on Haswell CPU.
