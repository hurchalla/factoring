// See LICENSE at bottom of this file.  FYI- this license is the same as
// https://github.com/danaj/Math-Prime-Util/blob/master/primality.c
// from which this file was derived.
//
// The contents of this experimental file were adapted from Dana Jacobsen's
// original file primality.c and his supporting header files, at
// https://github.com/danaj/Math-Prime-Util/blob/master/primality.c
// His primality.c file and all supporting header files are in his git
// repository Math-Prime-Util.

#ifndef EXPERIMENTAL_BPSW_DANA_LUCAS_H_INCLUDED
#define EXPERIMENTAL_BPSW_DANA_LUCAS_H_INCLUDED


#include <cstdint>
#include <math.h>
#include <assert.h>



// In this experimental project, the macros  LUCAS_FORCE_INLINE  NUM_LIMITS
// TRAIT_MAKE_SIGNED  TRAIT_MAKE_UNSIGNED  GCD_FUNCTION_NAME   are defined in
// benchmark_bpsw.cpp, which includes this header file.
// However -
// These macros could instead be defined by enabling the section below.
#if 0
#  include <limits>
#  include <type_traits>
#  include <numeric>
#  ifdef _MSC_VER
#    define LUCAS_FORCE_INLINE  __forceinline
#  else
#    define LUCAS_FORCE_INLINE  inline __attribute__((always_inline))
#  endif
#  define NUM_LIMITS  std::numeric_limits
#  define TRAIT_MAKE_SIGNED  std::make_signed
#  define TRAIT_MAKE_UNSIGNED  std::make_unsigned
#  define GCD_FUNCTION_NAME  std::gcd
#endif




template <typename UV>
LUCAS_FORCE_INLINE static UV gcd_ui(UV x, UV y) {
#if 1
  return GCD_FUNCTION_NAME(x, y);
#else
// dana's code
  UV t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
#endif
}


template <typename UV>
LUCAS_FORCE_INLINE static UV isqrt(UV n) {
  static_assert(NUM_LIMITS<UV>::is_integer, "");
  static_assert(!NUM_LIMITS<UV>::is_signed, "");
/*
#if BITS_PER_WORD == 32
  if (n >= UVCONST(4294836225)) return UVCONST(65535);
#else
  if (n >= UVCONST(18446744065119617025)) return UVCONST(4294967295);
#endif
*/
  static constexpr UV UV_DIGITS = NUM_LIMITS<UV>::digits;
  static_assert(UV_DIGITS >= 32);
  if constexpr (UV_DIGITS == 32) {
      if (n >= UINT32_C(4294836225)) return UINT32_C(65535);
  } else {
      static_assert(UV_DIGITS == 64);
      if (n >= UINT64_C(18446744065119617025)) return UINT64_C(4294967295);
  }
  UV root = (UV) sqrt((double)n);
  if (root*root > n)  root--;
  if ((root+1)*(root+1) <= n)  root++;
  return root;
}


template <typename UV>
LUCAS_FORCE_INLINE static int is_perfect_square(UV n)
{
  static_assert(NUM_LIMITS<UV>::is_integer, "");
  static_assert(!NUM_LIMITS<UV>::is_signed, "");
  /* Step 1, reduce to 18% of inputs */
  uint32_t m = n & 127;
  if ((m*0x8bc40d7d) & (m*0xa1e2f5d1) & 0x14020a)  return 0;
  /* Step 2, reduce to 7% of inputs (mod 99 reduces to 4% but slower) */
  m = n %240; if ((m*0xfa445556) & (m*0x8021feb1) & 0x614aaa0f) return 0;
  /* m = n % 99; if ((m*0x5411171d) & (m*0xe41dd1c7) & 0x80028a80) return 0; */
  /* Step 3, do the square root instead of any more rejections */
  m = isqrt(n);
  return (UV)m*(UV)m == n;
}


template <typename UV>
static int jacobi_iu(typename TRAIT_MAKE_SIGNED<UV>::type in, UV m) {
  static_assert(NUM_LIMITS<UV>::is_integer, "");
  static_assert(!NUM_LIMITS<UV>::is_signed, "");

  int j = 1;
  UV n = (in < 0) ? -in : in;

  if (m <= 0 || (m%2) == 0) return 0;
  if (in < 0 && (m%4) == 3) j = -j;
  while (n != 0) {
    while ((n % 2) == 0) {
      n >>= 1;
      if ( (m % 8) == 3 || (m % 8) == 5 )  j = -j;
    }
    { UV t = n; n = m; m = t; }
    if ( (n % 4) == 3 && (m % 4) == 3 )  j = -j;
    n = n % m;
  }
  return (m == 1) ? j : 0;
}



template <typename UV>
LUCAS_FORCE_INLINE static UV select_extra_strong_parameters(UV n, UV increment) {
  static_assert(NUM_LIMITS<UV>::is_integer, "");
  static_assert(!NUM_LIMITS<UV>::is_signed, "");
  int j;
  UV D, P = 3;
  while (1) {
    D = P*P - 4;
    j = jacobi_iu(D, n);
    if (j == 0) { UV g = gcd_ui(D,n);  if (g != 1 && g != n) return 0; }
    if (j == -1) break;
    if (P == (3+20*increment) && is_perfect_square(n)) return 0;
    P += increment;
    assert(P <= 65535); //      croak("lucas_extrastrong_params: P exceeded 65535");
  }
  if (P >= n)  P %= n;   /* Never happens with increment < 4 */
  return P;
}







// A generalization of Pari's shortcut to the extra-strong Lucas test.
//
// This only calculates and tests V, which means less work, but it does result
// in a few more pseudoprimes than the full extra-strong test.
//
// I've added a gcd check at the top, which needs to be done and also results
// in fewer pseudoprimes.  Pari always does trial division to 100 first so
// is unlikely to come up there.
//
// IAESLP_INCREMENT:  set to 1 for Baillie OEIS, 2 for Pari.
//
// With increment = 1, these results will be a subset of the extra-strong
// Lucas pseudoprimes.  With increment = 2, we produce Pari's results.


// Note by JH 1/26/25:
// I suspect it would be more accurate to call this function
// is_almost_extra_strong_lucas_probableprime(), because it appears to test
// whether n is a lucas probable prime, rather than a lucas pseudoprime.
// For consistency I retain the original name Dana used for it. 


template <typename MontType, unsigned char IAESLP_INCREMENT = 1>
LUCAS_FORCE_INLINE bool is_almost_extra_strong_lucas_pseudoprime(
      const MontType& mf, const typename MontType::CanonicalValue& mont2)
{
  using T = typename MontType::IntegerType;
  using UV = typename TRAIT_MAKE_UNSIGNED<T>::type;

  static_assert(NUM_LIMITS<UV>::is_integer, "");
  static_assert(!NUM_LIMITS<UV>::is_signed, "");
  UV n = mf.getModulus();

  auto UV_MAX = NUM_LIMITS<UV>::max();

// this is already checked by the function that calls us:
//    if (n < 13) return (n == 2 || n == 3 || n == 5 || n == 7 || n == 11);
//    if ((n % 2) == 0 || n == UV_MAX) return false;
  assert(n >= 13);
  assert((n % 2) != 0 && n != UV_MAX);

  static constexpr UV increment = IAESLP_INCREMENT;

//  if (increment < 1 || increment > 256)
//    croak("Invalid lucas parameter increment: %"UVuf"\n", increment);
  static_assert(increment == 1 || increment == 2, "");

  // Ensure small primes work with large increments.
//  if ( (increment >= 16 && n <= 331) || (increment > 148 && n <= 631) )
//    return !!is_prob_prime(n);


// --------------------------------------------------------------------

// start of test as in BPSW (which sets increment = 1):


  UV P = select_extra_strong_parameters(n, increment);
  if (P == 0) return false;

  assert(n % 2 == 1);
  UV d = n+1;
  UV s = 0;
  assert((d & 1) == 0);
  assert(d > 0);
//  while ( (d & 1) == 0 ) {  s++;  d >>= 1; }
  do { s++;  d >>= 1; } while ( (d & 1) == 0 );

  using M = typename MontType::MontgomeryValue;
  using C = typename MontType::CanonicalValue;
  using FV = typename MontType::FusingValue;

//  C mont2 = mf.add(mf.getUnityValue(), mf.getUnityValue());

    M V;
    {
//      const uint64_t montP = mont_geta(P, n);
      const M montP = mf.convertIn(P);
      const FV montPfv = mf.getFusingValue(montP);
//      W = submod(  mont_mulmod( montP, montP, n),  mont2, n);
      M W = mf.fusedSquareSub(montP, mont2);
      V = montP;

      UV b;
      { UV v = d; b = 0; while (v >>= 1) b++; }

      while (b > 0) {
        b--;

//        UV mT = submod(  mont_mulmod(V, W, n),  montP, n);
        M mT = mf.fmsub(V, W, montPfv);
#if 1
        if ( (d >> b) & 1 ) {
          V = mT;
//          W = submod(  mont_mulmod(W, W, n),  mont2, n);
          W = mf.fusedSquareSub(W, mont2);
        } else {
          W = mT;
//          V = submod(  mont_mulmod(V, V, n),  mont2, n);
          V = mf.fusedSquareSub(V, mont2);
        }
#else
// in theory this #else should be faster, but it is not with MSVC.

        M tmpW = mf.fusedSquareSub(W, mont2);
        M tmpV = mf.fusedSquareSub(V, mont2);
  #if 0
        UV shift = (d >> b);
        W = (shift & 1) ? tmpW : mT;
        V = (shift & 1) ? mT : tmpV;
  #else
        //UV shift = (d >> b);
        W = mT;
        W.cmov(( (d >> b) & 1 ), tmpW);
        V = tmpV;
        V.cmov(( (d >> b) & 1 ), mT);
  #endif
#endif
      }
    }

//    if (V == mont2 || V == (n-mont2))
    C cV = mf.getCanonicalValue(V);
    if (cV == mont2 || cV == mf.negate(mont2))
      return true;

    C zero = mf.getZeroValue();
    while (s > 1) {
      s--;
      if (cV == zero)
        return true;
//      V = submod(  mont_mulmod(V, V, n),  mont2, n);
      V = mf.fusedSquareSub(V, mont2);
      cV = mf.getCanonicalValue(V);

//#define IAESLP_USE_BPSW_LINES
// In Dana's original files, BPSW() has these lines, and
// is_almost_extra_strong_lucas_pseudoprime() does not.
// Note: these lines should make no difference to the return value, since
// if cV == mont2, then cV will forever remain == mont2 within this loop
// (and never == zero).
#ifdef IAESLP_USE_BPSW_LINES
      if (cV == mont2)
        return false;
#endif

    }

    return false;
}






// Lucas tests.  The following are the valid values for LUCAS_STRENGTH -
//  0: Standard
//  1: Strong
//  2: Stronger (Strong + page 1401 extra tests)
//  3: Extra Strong (Mo/Jones/Grantham)
// None of them have any false positives for the BPSW test.  Also see the
// "almost extra strong" test.

// Note by JH 1/26/25:
// I suspect it would be more accurate to call this function
// is_lucas_probableprime(), because it appears to test
// whether n is a lucas probable prime, rather than a lucas pseudoprime.
// For consistency I retain the original name Dana used.
//
template <unsigned int LUCAS_STRENGTH, typename MontType>
LUCAS_FORCE_INLINE bool is_lucas_pseudoprime(const MontType& mf,
                 const typename MontType::CanonicalValue& mont2)
{
  using T = typename MontType::IntegerType;
  using UV = typename TRAIT_MAKE_UNSIGNED<T>::type;
  using IV = typename TRAIT_MAKE_SIGNED<T>::type;
  static_assert(NUM_LIMITS<UV>::is_integer, "");
  static_assert(!NUM_LIMITS<UV>::is_signed, "");
  static_assert(NUM_LIMITS<IV>::is_integer, "");
  static_assert(NUM_LIMITS<IV>::is_signed, "");

  auto UV_MAX = NUM_LIMITS<UV>::max();

  // n is the number to test
  UV n = mf.getModulus();

  IV P, Q, D;
  UV Pu, Qu, Qk1, d, s;

  if (n < 5) return (n == 2 || n == 3);
  if ((n % 2) == 0 || n == UV_MAX) return 0;

  if (LUCAS_STRENGTH < 3) {
    UV Du = 5;
    IV sign = 1;
    int j;
    while (1) {
      D = Du * sign;
      j = jacobi_iu(D, n);
      if (j != 1 && Du != n) break;
      if (Du == 21 && is_perfect_square(n)) return 0;
      Du += 2;
      sign = -sign;
    }
    if (j != -1) return 0;
    P = 1;
    Q = (1 - D) / 4;
    if (LUCAS_STRENGTH == 2 && Q == -1) P=Q=D=5;  /* Method A* */
    /* Check gcd(n,2QD). gcd(n,2D) already done. */
    Qk1 = (Q >= 0)  ?  Q % n  :  n-(((UV)(-Q)) % n);
    if (gcd_ui(Qk1,n) != 1) return 0;
  } else {
    P = select_extra_strong_parameters(n, (UV)1);
    if (P == 0) return 0;
    Q = 1;
    D = P*P - 4;
  }
  assert(D == (P*P - 4*Q));  // "is_lucas_pseudoprime: incorrect DPQ"

#if 0   // Condition 2, V_n+1 = 2Q mod n
{ UV us, vs; lucasuvmod(&us, &vs, P, Q, n+1, n); return (vs == addmod(Q,Q,n)); }
#endif
#if 0   // Condition 3, n is a epsp(Q)
return is_euler_pseudoprime(n,Qk1);
#endif

  d = n+1;
  s = 0;
  if (LUCAS_STRENGTH > 0)
    while ( (d & 1) == 0 ) {  s++;  d >>= 1; }

  using M = typename MontType::MontgomeryValue;
  using C = typename MontType::CanonicalValue;
  using FV = typename MontType::FusingValue;

  {
//    const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);
//    const uint64_t mont2 = mont_get2(n);
    const C mont1 = mf.getUnityValue();
    const C montneg2 = mf.negate(mont2);
/*
    const uint64_t montP = (P == 1) ? mont1
                         : (P >= 0) ? mont_geta(P, n)
                         : n - mont_geta(-P, n);
    const uint64_t montQ = (Q == 1) ? mont1
                         : (Q >= 0) ? mont_geta(Q, n)
                         : n - mont_geta(-Q, n);
    const uint64_t montD = (D >= 0) ? mont_geta(D, n)
                         : n - mont_geta(-D, n);
*/
    const M montP = (P >= 0) ? mf.convertIn(P) : mf.negate(mf.convertIn(-P));
    const M montQ = (Q >= 0) ? mf.convertIn(Q) : mf.negate(mf.convertIn(-Q));
    const M montD = (D >= 0) ? mf.convertIn(D) : mf.negate(mf.convertIn(-D));

    UV b;
    { UV v = d; b = 0; while (v >>= 1) b++; }

    /* U, V, Qk, and mont* are in Montgomery space */
    M U = mont1;
    M V = montP;
    M Qk;

    if (Q == 1 || Q == -1) {   /* Faster code for |Q|=1, also opt for P=1 */
      int sign = Q;
      while (b--) {
//        U = mont_mulmod(U, V, n);
        U = mf.multiply(U, V);

//        if (sign == 1) V = submod( mont_sqrmod(V,n), mont2, n);
//        else           V = addmod( mont_sqrmod(V,n), mont2, n);
        C tmpMont = montneg2;
        tmpMont.cmov((sign == 1), mont2);
        V = mf.fusedSquareSub(V, tmpMont);

        sign = 1;
        if ( (d >> b) & 1 ) {
          // UV t2 = mont_mulmod(U, montD, n);
          M t2 = mf.multiply(U, montD);
          if (P == 1) {
//            U = addmod(U, V, n);
//            V = addmod(V, t2, n);
            U = mf.add(U, V);
            V = mf.add(V, t2);
          } else {
//            U = addmod( mont_mulmod(U, montP, n), V, n);
//            V = addmod( mont_mulmod(V, montP, n), t2, n);
            U = mf.fmadd(U, montP, mf.getFusingValue(V));
            V = mf.fmadd(V, montP, mf.getFusingValue(t2));
          }

// JH note 1/26/25:
// The translation Dana made of the two lines below into Mont space was
// surprising.  On the surface it seems like an inaccurate translation to Mont
// space, but in my limited testing, I got correct end-results with it.
// (Within Dana's Math-Prime-Util, you can see the equivalent non-montgomery
// code inside lucasuvmod(), in lucas_seq.c.).
// The reason the translation is surprising is that the original two lines do
// not have the same meaning if used directly in Montgomery space.
// The original two lines were:
//          if (U & 1) { U = (n>>1) + (U>>1) + 1; } else { U >>= 1; }
//          if (V & 1) { V = (n>>1) + (V>>1) + 1; } else { V >>= 1; }
// For his Montgomery code, Dana uses these lines unchanged (from non-montgomery
// space).
// ... (In non-Mont space) it looks like maybe(?) the two lines are meant to
// redistribute odd numbers to [n/2,n] and even numbers to [0,n/2], but as of
// yet I haven't studied the Lucas alg, so I don't know the reason for the code.
//
// A simple (and slow) way to retain the meaning of the original lines when
// translating to Mongomery space, would be to convert U and V out of Mont space
// and run the two lines, and then convert U and V back into Mont space.
// You can uncomment the following #define if you want to do this, which retains
// the exact original meaning (though note that converting out/in is slow):
//
//#define EXACT_TRANSLATION_TO_MONT_SPACE
//
#ifdef EXACT_TRANSLATION_TO_MONT_SPACE
// this retains the same meaning as the non-montgomery code, but it's slow.
          auto Utmp = mf.convertOut(U);
          auto Vtmp = mf.convertOut(V);
          if (Utmp & 1) { Utmp = (n>>1) + (Utmp>>1) + 1; } else { Utmp >>= 1; }
          if (Vtmp & 1) { Vtmp = (n>>1) + (Vtmp>>1) + 1; } else { Vtmp >>= 1; }
          U = mf.convertIn(Utmp);
          V = mf.convertIn(Vtmp);
#else
// this is fast and seems to work in my limited testing, but it doesn't retain
// the exact same meaning as the original code coming from non-montgomery space.
// Is it 100% correct?  It may be - I don't know.
          struct OpenMV : public MontType::MontgomeryValue {
            using MV = typename MontType::MontgomeryValue;
            LUCAS_FORCE_INLINE explicit OpenMV(MV x) : MV(x) {}
            LUCAS_FORCE_INLINE explicit OpenMV(T a) : MV(a) {}
            LUCAS_FORCE_INLINE T get() const { return MV::get(); }
          };
          T Ut = OpenMV(U).get();
          T Vt = OpenMV(V).get();
          if (Ut & 1) { Ut = (n>>1) + (Ut>>1) + 1; } else { Ut >>= 1; }
          if (Vt & 1) { Vt = (n>>1) + (Vt>>1) + 1; } else { Vt >>= 1; }
          U = OpenMV(Ut);
          V = OpenMV(Vt);
#endif
          sign = Q;
        }
      }
      Qk = (sign == 1) ? mont1 : mf.negate(mont1);
    } else {
      Qk = montQ;
      while (b--) {
//        U = mont_mulmod(U, V, n);
        U = mf.multiply(U, V);
//        V = submod( mont_sqrmod(V,n), addmod(Qk,Qk,n), n);
        V = mf.fusedSquareSub(V, mf.getCanonicalValue(mf.add(Qk,Qk)));
//        Qk = mont_sqrmod(Qk,n);
        Qk = mf.square(Qk);
        if ( (d >> b) & 1 ) {
//          UV t2 = mont_mulmod(U, montD, n);
          M t2 = mf.multiply(U, montD);
//          U = addmod( mont_mulmod(U, montP, n), V, n);
          U = mf.fmadd(U, montP, mf.getFusingValue(V));
//          V = addmod( mont_mulmod(V, montP, n), t2, n);
          V = mf.add(mf.multiply(V, montP), t2);

// JH note 1/26/25:
// See comments above for EXACT_TRANSLATION_TO_MONT_SPACE.
#ifdef EXACT_TRANSLATION_TO_MONT_SPACE
          auto Utmp = mf.convertOut(U);
          auto Vtmp = mf.convertOut(V);
          if (Utmp & 1) { Utmp = (n>>1) + (Utmp>>1) + 1; } else { Utmp >>= 1; }
          if (Vtmp & 1) { Vtmp = (n>>1) + (Vtmp>>1) + 1; } else { Vtmp >>= 1; }
          U = mf.convertIn(Utmp);
          V = mf.convertIn(Vtmp);
#else
          struct OpenMV : public MontType::MontgomeryValue {
            using MV = typename MontType::MontgomeryValue;
            LUCAS_FORCE_INLINE explicit OpenMV(MV x) : MV(x) {}
            LUCAS_FORCE_INLINE explicit OpenMV(T a) : MV(a) {}
            LUCAS_FORCE_INLINE T get() const { return MV::get(); }
          };
          T Ut = OpenMV(U).get();
          T Vt = OpenMV(V).get();
          if (Ut & 1) { Ut = (n>>1) + (Ut>>1) + 1; } else { Ut >>= 1; }
          if (Vt & 1) { Vt = (n>>1) + (Vt>>1) + 1; } else { Vt >>= 1; }
          U = OpenMV(Ut);
          V = OpenMV(Vt);
#endif

//          Qk = mont_mulmod(Qk, montQ, n);
          Qk = mf.multiply(Qk, montQ);
        }
      }
    }

    C zero = mf.getZeroValue();

    if (LUCAS_STRENGTH == 0) {
      if (mf.getCanonicalValue(U) == zero)
        return 1;
    } else if (LUCAS_STRENGTH == 1) {
      if (mf.getCanonicalValue(U) == zero)
        return 1;
      while (s--) {
        if (mf.getCanonicalValue(V) == zero)
          return 1;
        if (s) {
//          V = submod( mont_sqrmod(V,n), addmod(Qk,Qk,n), n);
          C Qc = mf.getCanonicalValue(Qk);
          V = mf.fusedSquareSub(V, mf.add(Qc,Qc));
//          Qk = mont_sqrmod(Qk,n);
          Qk = mf.square(Qk);
        }
      }
    } else if (LUCAS_STRENGTH == 2) {
      M Ql = zero;
      C Qj = zero;
      int qjacobi, is_slpsp = 0;
      if (mf.getCanonicalValue(U) == zero)
        is_slpsp = 1;
      while (s--) {
        if (mf.getCanonicalValue(V) == zero)
          is_slpsp = 1;
        Ql = Qk;
//        V = submod( mont_sqrmod(V,n), addmod(Qk,Qk,n), n);
        C Qc = mf.getCanonicalValue(Qk);
        V = mf.fusedSquareSub(V, mf.add(Qc,Qc));
//        Qk = mont_sqrmod(Qk,n);
        Qk = mf.square(Qk);
      }
      if (!is_slpsp)                  return 0; /* slpsp */
      C montQc = mf.getCanonicalValue(montQ);
      if (mf.getCanonicalValue(V) != mf.add(montQc,montQc)) return 0; /* V_{n+1} != 2Q mod n */
      qjacobi = jacobi_iu(Q,n);
      Qj = (qjacobi == 0) ? zero : (qjacobi == 1) ? montQc : mf.negate(montQc);
      if (mf.getCanonicalValue(Ql) != Qj)                   return 0; /* n is epsp base Q */
      return 1;
    } else {
      C Vc = mf.getCanonicalValue(V);
      if ( mf.getCanonicalValue(U) == zero && (Vc == mont2 || Vc == montneg2) )
        return 1;
      s--;
      while (s--) {
        if (Vc == zero)
          return 1;
        if (s) {
//          V = submod( mont_sqrmod(V,n), mont2, n);
          V = mf.fusedSquareSub(V, mont2);
          Vc = mf.getCanonicalValue(V);
        }
      }
    }
    return 0;
  }
}




// Lucas tests.  The following are the valid values for LUCAS_STRENGTH -
//  0: Standard
//  1: Strong
//  2: Stronger (Strong + page 1401 extra tests)
//  3: Extra Strong (Mo/Jones/Grantham)
//  4: Almost-Extra-Strong
template <unsigned int LUCAS_STRENGTH, typename MontType>
LUCAS_FORCE_INLINE bool is_lucas_probable_prime(
        const MontType& mf, const typename MontType::CanonicalValue& mont2)
{
  static_assert(LUCAS_STRENGTH >= 0 && LUCAS_STRENGTH <= 4);
  if constexpr (LUCAS_STRENGTH == 4) {
    constexpr unsigned char IAESLP_INCREMENT = 1;
    return is_almost_extra_strong_lucas_pseudoprime<MontType, IAESLP_INCREMENT>(mf, mont2);
  } else {
    return is_lucas_pseudoprime<LUCAS_STRENGTH>(mf, mont2);
  }
}




// -------------------------------------------------------------
// LICENSE
// -------------------------------------------------------------

/*
This software is Copyright (c) 2011-2021 by Dana Jacobsen,
and (c) 2025 by Jeff Hurchalla.

This is free software; you can redistribute it and/or modify it under
the same terms as the Perl 5 programming language system itself.

Terms of the Perl programming language system itself

a) the GNU General Public License as published by the Free
   Software Foundation; either version 1, or (at your option) any
   later version, or
b) the "Artistic License"

--- The GNU General Public License, Version 1, February 1989 ---

This software is Copyright (c) 2011-2021 by Dana Jacobsen.

This is free software, licensed under:

  The GNU General Public License, Version 1, February 1989

                    GNU GENERAL PUBLIC LICENSE
                     Version 1, February 1989

 Copyright (C) 1989 Free Software Foundation, Inc.
 51 Franklin St, Suite 500, Boston, MA  02110-1335  USA

 Everyone is permitted to copy and distribute verbatim copies
 of this license document, but changing it is not allowed.

                            Preamble

  The license agreements of most software companies try to keep users
at the mercy of those companies.  By contrast, our General Public
License is intended to guarantee your freedom to share and change free
software--to make sure the software is free for all its users.  The
General Public License applies to the Free Software Foundation's
software and to any other program whose authors commit to using it.
You can use it for your programs, too.

  When we speak of free software, we are referring to freedom, not
price.  Specifically, the General Public License is designed to make
sure that you have the freedom to give away or sell copies of free
software, that you receive source code or can get it if you want it,
that you can change the software or use pieces of it in new free
programs; and that you know you can do these things.

  To protect your rights, we need to make restrictions that forbid
anyone to deny you these rights or to ask you to surrender the rights.
These restrictions translate to certain responsibilities for you if you
distribute copies of the software, or if you modify it.

  For example, if you distribute copies of a such a program, whether
gratis or for a fee, you must give the recipients all the rights that
you have.  You must make sure that they, too, receive or can get the
source code.  And you must tell them their rights.

  We protect your rights with two steps: (1) copyright the software, and
(2) offer you this license which gives you legal permission to copy,
distribute and/or modify the software.

  Also, for each author's protection and ours, we want to make certain
that everyone understands that there is no warranty for this free
software.  If the software is modified by someone else and passed on, we
want its recipients to know that what they have is not the original, so
that any problems introduced by others will not reflect on the original
authors' reputations.

  The precise terms and conditions for copying, distribution and
modification follow.

                    GNU GENERAL PUBLIC LICENSE
   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION

  0. This License Agreement applies to any program or other work which
contains a notice placed by the copyright holder saying it may be
distributed under the terms of this General Public License.  The
"Program", below, refers to any such program or work, and a "work based
on the Program" means either the Program or any work containing the
Program or a portion of it, either verbatim or with modifications.  Each
licensee is addressed as "you".

  1. You may copy and distribute verbatim copies of the Program's source
code as you receive it, in any medium, provided that you conspicuously and
appropriately publish on each copy an appropriate copyright notice and
disclaimer of warranty; keep intact all the notices that refer to this
General Public License and to the absence of any warranty; and give any
other recipients of the Program a copy of this General Public License
along with the Program.  You may charge a fee for the physical act of
transferring a copy.

  2. You may modify your copy or copies of the Program or any portion of
it, and copy and distribute such modifications under the terms of Paragraph
1 above, provided that you also do the following:

    a) cause the modified files to carry prominent notices stating that
    you changed the files and the date of any change; and

    b) cause the whole of any work that you distribute or publish, that
    in whole or in part contains the Program or any part thereof, either
    with or without modifications, to be licensed at no charge to all
    third parties under the terms of this General Public License (except
    that you may choose to grant warranty protection to some or all
    third parties, at your option).

    c) If the modified program normally reads commands interactively when
    run, you must cause it, when started running for such interactive use
    in the simplest and most usual way, to print or display an
    announcement including an appropriate copyright notice and a notice
    that there is no warranty (or else, saying that you provide a
    warranty) and that users may redistribute the program under these
    conditions, and telling the user how to view a copy of this General
    Public License.

    d) You may charge a fee for the physical act of transferring a
    copy, and you may at your option offer warranty protection in
    exchange for a fee.

Mere aggregation of another independent work with the Program (or its
derivative) on a volume of a storage or distribution medium does not bring
the other work under the scope of these terms.

  3. You may copy and distribute the Program (or a portion or derivative of
it, under Paragraph 2) in object code or executable form under the terms of
Paragraphs 1 and 2 above provided that you also do one of the following:

    a) accompany it with the complete corresponding machine-readable
    source code, which must be distributed under the terms of
    Paragraphs 1 and 2 above; or,

    b) accompany it with a written offer, valid for at least three
    years, to give any third party free (except for a nominal charge
    for the cost of distribution) a complete machine-readable copy of the
    corresponding source code, to be distributed under the terms of
    Paragraphs 1 and 2 above; or,

    c) accompany it with the information you received as to where the
    corresponding source code may be obtained.  (This alternative is
    allowed only for noncommercial distribution and only if you
    received the program in object code or executable form alone.)

Source code for a work means the preferred form of the work for making
modifications to it.  For an executable file, complete source code means
all the source code for all modules it contains; but, as a special
exception, it need not include source code for modules which are standard
libraries that accompany the operating system on which the executable
file runs, or for standard header files or definitions files that
accompany that operating system.

  4. You may not copy, modify, sublicense, distribute or transfer the
Program except as expressly provided under this General Public License.
Any attempt otherwise to copy, modify, sublicense, distribute or transfer
the Program is void, and will automatically terminate your rights to use
the Program under this License.  However, parties who have received
copies, or rights to use copies, from you under this General Public
License will not have their licenses terminated so long as such parties
remain in full compliance.

  5. By copying, distributing or modifying the Program (or any work based
on the Program) you indicate your acceptance of this license to do so,
and all its terms and conditions.

  6. Each time you redistribute the Program (or any work based on the
Program), the recipient automatically receives a license from the original
licensor to copy, distribute or modify the Program subject to these
terms and conditions.  You may not impose any further restrictions on the
recipients' exercise of the rights granted herein.

  7. The Free Software Foundation may publish revised and/or new versions
of the General Public License from time to time.  Such new versions will
be similar in spirit to the present version, but may differ in detail to
address new problems or concerns.

Each version is given a distinguishing version number.  If the Program
specifies a version number of the license which applies to it and "any
later version", you have the option of following the terms and conditions
either of that version or of any later version published by the Free
Software Foundation.  If the Program does not specify a version number of
the license, you may choose any version ever published by the Free Software
Foundation.

  8. If you wish to incorporate parts of the Program into other free
programs whose distribution conditions are different, write to the author
to ask for permission.  For software which is copyrighted by the Free
Software Foundation, write to the Free Software Foundation; we sometimes
make exceptions for this.  Our decision will be guided by the two goals
of preserving the free status of all derivatives of our free software and
of promoting the sharing and reuse of software generally.

                            NO WARRANTY

  9. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
REPAIR OR CORRECTION.

  10. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES,
INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING
OUT OF THE USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED
TO LOSS OF DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY
YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER
PROGRAMS), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.

                     END OF TERMS AND CONDITIONS

        Appendix: How to Apply These Terms to Your New Programs

  If you develop a new program, and you want it to be of the greatest
possible use to humanity, the best way to achieve this is to make it
free software which everyone can redistribute and change under these
terms.

  To do so, attach the following notices to the program.  It is safest to
attach them to the start of each source file to most effectively convey
the exclusion of warranty; and each file should have at least the
"copyright" line and a pointer to where the full notice is found.

    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 19yy  <name of author>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 1, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston MA  02110-1301 USA


Also add information on how to contact you by electronic and paper mail.

If the program is interactive, make it output a short notice like this
when it starts in an interactive mode:

    Gnomovision version 69, Copyright (C) 19xx name of author
    Gnomovision comes with ABSOLUTELY NO WARRANTY; for details type `show w'.
    This is free software, and you are welcome to redistribute it
    under certain conditions; type `show c' for details.

The hypothetical commands `show w' and `show c' should show the
appropriate parts of the General Public License.  Of course, the
commands you use may be called something other than `show w' and `show
c'; they could even be mouse-clicks or menu items--whatever suits your
program.

You should also get your employer (if you work as a programmer) or your
school, if any, to sign a "copyright disclaimer" for the program, if
necessary.  Here a sample; alter the names:

  Yoyodyne, Inc., hereby disclaims all copyright interest in the
  program `Gnomovision' (a program to direct compilers to make passes
  at assemblers) written by James Hacker.

  <signature of Ty Coon>, 1 April 1989
  Ty Coon, President of Vice

That's all there is to it!


--- The Artistic License 1.0 ---

This software is Copyright (c) 2011-2021 by Dana Jacobsen.

This is free software, licensed under:

  The Artistic License 1.0

The Artistic License

Preamble

The intent of this document is to state the conditions under which a Package
may be copied, such that the Copyright Holder maintains some semblance of
artistic control over the development of the package, while giving the users of
the package the right to use and distribute the Package in a more-or-less
customary fashion, plus the right to make reasonable modifications.

Definitions:

  - "Package" refers to the collection of files distributed by the Copyright
    Holder, and derivatives of that collection of files created through
    textual modification. 
  - "Standard Version" refers to such a Package if it has not been modified,
    or has been modified in accordance with the wishes of the Copyright
    Holder. 
  - "Copyright Holder" is whoever is named in the copyright or copyrights for
    the package. 
  - "You" is you, if you're thinking about copying or distributing this Package.
  - "Reasonable copying fee" is whatever you can justify on the basis of media
    cost, duplication charges, time of people involved, and so on. (You will
    not be required to justify it to the Copyright Holder, but only to the
    computing community at large as a market that must bear the fee.) 
  - "Freely Available" means that no fee is charged for the item itself, though
    there may be fees involved in handling the item. It also means that
    recipients of the item may redistribute it under the same conditions they
    received it. 

1. You may make and give away verbatim copies of the source form of the
Standard Version of this Package without restriction, provided that you
duplicate all of the original copyright notices and associated disclaimers.

2. You may apply bug fixes, portability fixes and other modifications derived
from the Public Domain or from the Copyright Holder. A Package modified in such
a way shall still be considered the Standard Version.

3. You may otherwise modify your copy of this Package in any way, provided that
you insert a prominent notice in each changed file stating how and when you
changed that file, and provided that you do at least ONE of the following:

  a) place your modifications in the Public Domain or otherwise make them
     Freely Available, such as by posting said modifications to Usenet or an
     equivalent medium, or placing the modifications on a major archive site
     such as ftp.uu.net, or by allowing the Copyright Holder to include your
     modifications in the Standard Version of the Package.

  b) use the modified Package only within your corporation or organization.

  c) rename any non-standard executables so the names do not conflict with
     standard executables, which must also be provided, and provide a separate
     manual page for each non-standard executable that clearly documents how it
     differs from the Standard Version.

  d) make other distribution arrangements with the Copyright Holder.

4. You may distribute the programs of this Package in object code or executable
form, provided that you do at least ONE of the following:

  a) distribute a Standard Version of the executables and library files,
     together with instructions (in the manual page or equivalent) on where to
     get the Standard Version.

  b) accompany the distribution with the machine-readable source of the Package
     with your modifications.

  c) accompany any non-standard executables with their corresponding Standard
     Version executables, giving the non-standard executables non-standard
     names, and clearly documenting the differences in manual pages (or
     equivalent), together with instructions on where to get the Standard
     Version.

  d) make other distribution arrangements with the Copyright Holder.

5. You may charge a reasonable copying fee for any distribution of this
Package.  You may charge any fee you choose for support of this Package. You
may not charge a fee for this Package itself. However, you may distribute this
Package in aggregate with other (possibly commercial) programs as part of a
larger (possibly commercial) software distribution provided that you do not
advertise this Package as a product of your own.

6. The scripts and library files supplied as input to or produced as output
from the programs of this Package do not automatically fall under the copyright
of this Package, but belong to whomever generated them, and may be sold
commercially, and may be aggregated with this Package.

7. C or perl subroutines supplied by you and linked into this Package shall not
be considered part of this Package.

8. The name of the Copyright Holder may not be used to endorse or promote
products derived from this software without specific prior written permission.

9. THIS PACKAGE IS PROVIDED "AS IS" AND WITHOUT ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED WARRANTIES OF
MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.

The End
*/


#endif
