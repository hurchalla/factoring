/*
Copyright (c) 2022, Jeff Hurchalla

Original source file prior to modifications was:
https://github.com/bbuhrow/yafu/blob/25b65990d6501b0a71e69963fb59c1fc4ab28df1/factor/gmp-ecm/microecm.c

IMPORTANT
Currently (and temporarily) this file and all of its modifications from the
original file are unavailable under any license.  For now, you are explicitly
*NOT* allowed to use, copy, distribute, modify, or share this software for any
purpose, without permission from the author.

You can expect this file will soon have a normal permissive license.
*/

/*
The following is reproduced to comply with the original source file:

Copyright (c) 2014, Ben Buhrow
All rights reserved.

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.
*/


#ifndef HURCHALLA_FACTORING_MICRO_ECM_CPP_H_INCLUDED
#define HURCHALLA_FACTORING_MICRO_ECM_CPP_H_INCLUDED


#include <stddef.h>
#include <stdint.h>
#include <type_traits>

#ifdef HURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM
#include "DualMontgomeryForm.h"
#endif
#include "hurchalla/modular_arithmetic/modular_multiplicative_inverse.h"
#include "hurchalla/factoring/greatest_common_divisor.h"
#include "hurchalla/util/count_trailing_zeros.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"


namespace hurchalla { namespace detail {



template <class MF>
struct uecm_mfpt
{
    using MV = typename MF::MontgomeryValue;
    MV X;
    MV Z;
};


// Implementation note on DUMMY param:
// The DUMMY template parameter ensures that if this struct is never used, then
// its static array will not exist and thus will use no memory.  We enforce
// DUMMY to always be void via static_assert, which ensures that only one
// instantiation is possible, and thus no more than one copy of the static array
// can ever exist.
template <typename DUMMY=void>
struct uecm_steps80 {
    static_assert(std::is_same<DUMMY, void>::value, "");
    static constexpr uint8_t steps[] = {
        9,                // 5
        11,9,             // 7
        10,9,             // 11
        8,8,9,            // 13
        9,                // 5 (2nd time)
        8,10,9,           // 17
        8,11,8,9,         // 19
        10,10,9,          // 23
        10,8,8,9,         // 29
        8,8,11,8,9,       // 31

        10,11,8,11,9,     // 37
//        11,11,8,10,9,     // 37 alternate- either one seems fine
        10,11,8,8,9,      // 41
        8,11,8,10,9,      // 43
        10,8,8,8,9,       // 47
        11,9,             // 7 (2nd time)
        10,11,11,8,8,9,   // 53
        10,11,8,10,9,     // 59
        10,10,8,8,9,      // 61
        9,                // 5 (3rd time)
        11,8,8,8,10,9,    // 67

        11,8,10,8,8,9,    // 71
        8,8,8,8,10,9,     // 73
        8,8,8,11,8,8,9,   // 79  (can potentially swap with the line below)
        10,9              // 11 (2nd time)
    };
    static constexpr int num_primes = 24;
};


struct uecm_prime_and_ratio {
    uint8_t prime;
    float ratio;
};

template <typename DUMMY=void>  //for DUMMY rationale, see uecm_steps80 comments
struct uecm_upracparams {
    static_assert(std::is_same<DUMMY, void>::value, "");
    static constexpr uecm_prime_and_ratio params[] = {
        {  83, 0.548409048446403258 },
        {  89, 0.618033988749894903 },
        {  97, 0.723606797749978936 },
        { 101, 0.556250337855490828 },
        { 103, 0.632839806088706269 },
        { 107, 0.580178728295464130 },
        { 109, 0.548409048446403258 },
        { 113, 0.618033988749894903 },
        { 127, 0.548409048446403258 },
        { 131, 0.618033988749894903 },

        {  13, 0.618033988749894903 },   // 2nd time using 13
        { 137, 0.548409048446403258 },
        { 139, 0.520000000000000000 },
        { 149, 0.580178728295464130 },
        { 151, 0.650000000000000000 },
        { 157, 0.640157392785047019 },
        { 163, 0.551390822543526449 },
        { 167, 0.580178728295464130 },
        { 173, 0.612429949509495031 },
        {   5, 0.618033988749894903 },   // 4th time using 5

        { 179, 0.618033988749894903 },
        { 181, 0.551390822543526449 },
        { 191, 0.618033988749894903 },
        { 193, 0.618033988749894903 },
        { 197, 0.540000000000000000 },
        { 199, 0.551390822543526449 },
        {   7, 0.618033988749894903 },   // 3rd time using 7
        { 211, 0.504000000000000000 },
// spot for bits 64
    };

//    static constexpr int num_params = 28;
    static constexpr int num_params = sizeof(params)/sizeof(params[0]);

    // Note that there were some composites that were used in the original code,
    // that we don't use.  They didn't help performance in my tests and they
    // require more complex code.  These are the composites, for reference:
    //    uprac(mf, P,  7747, 0.552188778811121, s);  //  61 x 127
    //    uprac(mf, P, 14111, 0.632839806088706, s);  // 103 x 137
    //    uprac(mf, P, 20989, 0.620181980807415, s);  // 139 x 151
    //    uprac(mf, P, 29353, 0.580178728295464, s);  // 149 x 197
};



// extra stage1 primes array for 128 bit types (unused for <= 64 bit types)
template <typename DUMMY=void>  //for DUMMY rationale, see uecm_steps80 comments
struct uecm_upracparams2 {
    static_assert(std::is_same<DUMMY, void>::value, "");
    static constexpr uint16_t primes[] = {
        223, 227, 229, 233, 239, 241, 17,
//7.  spot for bits 66
                                           251, 257, 263,
        269, 271, 277, 281, 283, 293, 307,  19, 311, 313,
        317, 331, 337, 347, 349, 353, 359, 367, 373,  23,
        379, 383, 389, 397, 401,  11, 409,
//37.  spot for bits 74
                                           419, 421, 431,
        433, 439, 443, 449, 457, 461, 463, 467, 479,   5,
        487, 491, 499, 503, 509, 521, 523, 541, 547, 557,
        29,  563, 569, 571, 577, 587, 593, 599, 601, 607,
        613, 617, 619, 631, 641, 643, 647, 653,  31,
//79.  spot for bits 82
                                                     659,
        661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
        739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
        811, 821,  37, 823, 827, 829, 839, 853, 857, 859,
        863, 877, 881, 883, 887, 907, 911, 919, 929, 937,
        941,  41, 947, 953, 967, 971, 977,   7, 983, 991,
        997,  13, 1009, 1013, 1019, 1021, 1031, 1033, 1039,
//139.  spot for bits 90
                                                              1049,
        1051, 1061, 1063, 1069,   43, 1087, 1091, 1093, 1097, 1103,
        1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187,
        1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259,
        1277, 1279,   47, 1283, 1289, 1291, 1297, 1301, 1303, 1307,
        1319, 1321, 1327, 1361,   5,  1367, 1373, 1381, 1399, 1409,
        1423, 1427, 1429, 1433, 1439, 1447, 1451, 1453, 1459, 1471,
        1481, 1483, 1487, 1489, 1493, 1499, 1511,   53, 1523, 1531,
        1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 1597, 1601,
        1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 1663, 1667,
        1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 1741, 1747,
          59, 1753, 1759,
//242.  spot for bits 100
                          1777, 1783, 1787, 1789, 1801, 1811, 1823,
        1831, 1847, 1861, 1867, 1871,   61, 1873, 1877, 1879,   17,
          11, 1889, 1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973,
        1979, 1987, 1993, 1997, 1999, 2003, 2011, 2017, 2027,   67,
        2029, 2039, 2053, 2063, 2069, 2081, 2083, 2087, 2089, 2099,
        2111, 2113, 2129, 2131, 2137, 2141, 2143, 2153, 2161,   71,
        2179, 2203, 2207, 2213, 2221, 2237, 2239, 2243, 2251, 2267,
        2269, 2273, 2281, 2287, 2293,   73, 2297, 2309, 2311, 2333,
        2339, 2341, 2347, 2351, 2357, 2371, 2377, 2381, 2383, 2389,
        2393, 2399, 2411, 2417, 2423, 2437, 2441, 2447, 2459,   79,
        2467,   19, 2473,    7, 2477, 2503,
//346.  spot for bits 108
                                            2521, 2531, 2539, 2543,
        2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 2621,   83,
        2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 2689,
        2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, 2749,
        2753, 2767, 2777, 2789, 2791,   89, 2797, 2801, 2803, 2819,
        2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
        2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
        3001, 3011, 3019, 3023, 3037, 3041,   97, 3049, 3061, 3067,
        3079, 3083, 3089,    5, 3109, 3119, 3121, 3137, 3163, 3167,
        3169, 3181, 3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251,
         101, 3253, 3257, 3259, 3271, 3299, 3301, 3307, 3313, 3319,
        3323, 3329, 3331, 3343,   13, 3347, 3359, 3361,  103, 3371,
        3373, 3389, 3391, 3407,   23, 3413, 3433, 3449, 3457, 3461,
        3463, 3467, 3469, 3491, 3499, 3511, 3517, 3527, 3529, 3533,
        3539, 3541,  107, 3547, 3557, 3559, 3571, 3581, 3583, 3593,
        3607, 3613, 3617, 3623,
//494.  spot for bits 116
                                3631, 3637, 3643, 3659, 3671, 3673,
        3677, 3691, 3697,  109, 3701, 3709, 3719, 3727, 3733, 3739,
        3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, 3823, 3833,
        3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, 3911, 3917,
        3919,  113, 3923, 3929, 3931, 3943, 3947, 3967, 3989, 4001,
        4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, 4073,
        4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, 4153,
        4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, 4241,
        4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, 4327,
        4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, 4421,
        4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 4507,
        4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,  127,
        4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
        4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
        4759, 4783, 4787, 4789, 4793, 4799, 4801,  131, 4813, 4817,
        4831, 4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933,
        4937, 4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999,
        5003, 5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081,
        5087, 5099,  137, 5101, 5107, 5113, 5119, 5147, 5153, 5167,
        5171, 5179, 5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261,
        5273, 5279, 5281, 5297,  139, 5303, 5309, 5323, 5333, 5347,
        5351, 5381, 5387, 5393, 5399,   29, 5407, 5413, 5417, 5419,
//710.  spot for bits 124
        5431, 5437, 5441, 5443, 5449, 5471, 5477, 5479, 5483,   11,
        5501, 5503, 5507, 5519, 5521, 5527, 5531, 5557, 5563, 5569,
        5573, 5581, 5591, 5623, 5639, 5641, 5647, 5651, 5653, 5657,
        5659, 5669, 5683,  149, 5689, 5693, 5701, 5711, 5717, 5737,
        5741, 5743, 5749, 5779, 5783, 5791, 5801, 5807, 5813, 5821,
        5827, 5839, 5843, 5849, 5851, 5857, 5861, 5867, 5869, 5879,
        5881, 5897,  151, 5903, 5923, 5927, 5939, 5953, 5981, 5987,
        6007, 6011, 6029, 6037, 6043, 6047, 6053, 6067, 6073, 6079,
        6089, 6091, 6101, 6113, 6121, 6131, 6133, 6143, 6151, 6163,
        6173, 6197, 6199, 6203, 6211, 6217, 6221, 6229, 6247, 6257,
        6263, 6269, 6271, 6277, 6287, 6299, 6301, 6311, 6317, 6323,
        6329, 6337, 6343, 6353, 6359, 6361, 6367, 6373, 6379, 6389,
        6397, 6421,  157, 6427, 6449, 6451, 6469, 6473, 6481, 6491,
        6521, 6529,   31, 6547, 6551, 6553, 6563, 6569, 6571, 6577,
        6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 6679,   // 7(??),
//860.  spot for bits 128
    };

    static constexpr int num_params = sizeof(primes)/sizeof(primes[0]);
    static_assert(num_params == 860);
};






// we use a struct with static member functions, in order to disallow ADL
struct micro_ecm {


static
inline uint32_t lcg_rand_32B(uint32_t lower, uint32_t upper, uint64_t *ploc_lcg)
{
    constexpr double INV_2_POW_32 = 1.0 / (double)((uint64_t)(1) << 32);
    *ploc_lcg = 6364136223846793005ULL * (*ploc_lcg) + 1442695040888963407ULL;
    return lower + (uint32_t)(
        (double)(upper - lower) * (double)((*ploc_lcg) >> 32) * INV_2_POW_32);
}



template <class MF>
static HURCHALLA_FORCE_INLINE
void inlined_uadd(const MF& mf, const uecm_mfpt<MF> P1, const uecm_mfpt<MF> P2,
                  const uecm_mfpt<MF> Pin, uecm_mfpt<MF> *Pout)
{
    using MV = typename MF::MontgomeryValue;
    // compute:
    //x+ = z- * [(x1-z1)(x2+z2) + (x1+z1)(x2-z2)]^2
    //z+ = x- * [(x1-z1)(x2+z2) - (x1+z1)(x2-z2)]^2
    // where:
    //x- = original x
    //z- = original z
    // given the sums and differences of the original points
    MV diff1 = mf.subtract(P1.X, P1.Z);
    MV sum1 = mf.add(P1.X, P1.Z);
    MV diff2 = mf.subtract(P2.X, P2.Z);
    MV sum2 = mf.add(P2.X, P2.Z);

    MV tt1 = mf.multiply(diff1, sum2); //U
    MV tt2 = mf.multiply(sum1, diff2); //V

    MV tt3 = mf.add(tt1, tt2);
    MV tt4 = mf.subtract(tt1, tt2);
    tt1 = mf.square(tt3);   //(U + V)^2
    tt2 = mf.square(tt4);   //(U - V)^2

    Pout->X = mf.multiply(tt1, Pin.Z);     //Z * (U + V)^2
    Pout->Z = mf.multiply(tt2, Pin.X);     //x * (U - V)^2

    return;
}
template <class MF>
static void uadd(const MF& mf, const uecm_mfpt<MF> P1, const uecm_mfpt<MF> P2,
                 const uecm_mfpt<MF> Pin, uecm_mfpt<MF> *Pout)
{
    return inlined_uadd(mf, P1, P2, Pin, Pout);
}


template <class MF>
static HURCHALLA_FORCE_INLINE
void inlined_udup(const MF& mf, typename MF::MontgomeryValue s,
                  uecm_mfpt<MF> point, uecm_mfpt<MF> *P)
{
    using MV = typename MF::MontgomeryValue;
    using CV = typename MF::CanonicalValue;
    MV indiff = mf.subtract(point.X, point.Z);
    MV insum = mf.add(point.X, point.Z);
    MV tt1 = mf.square(indiff);          // U=(x1 - z1)^2
    MV tt2 = mf.square(insum);           // V=(x1 + z1)^2
    P->X = mf.multiply(tt1, tt2);        // x=U*V

    CV cv1 = mf.getCanonicalValue(tt1);

    MV tt3 = mf.subtract(tt2, tt1);      // w = V-U
    tt2 = mf.fmadd(tt3, s, cv1);         // w = (A+2)/4 * w + U
    P->Z = mf.multiply(tt2, tt3);        // Z = w*(V-U)

    return;
}
template <class MF>
static void udup(const MF& mf, typename MF::MontgomeryValue s,
                 uecm_mfpt<MF> point, uecm_mfpt<MF> *P)
{
    return inlined_udup(mf, s, point, P);
}



template <class MF>
static void uprac(const MF& mf, uecm_mfpt<MF> *P, uint64_t c, double v,
                  typename MF::MontgomeryValue s)
{
    HPBC_PRECONDITION2(c > 0);

    using MV = typename MF::MontgomeryValue;
    uint64_t d, e;
    int shift = hurchalla::count_trailing_zeros(c);
    c = c >> shift;
    d = c;
    {
        uint64_t r = (uint64_t)((double)d * v + 0.5);
        d = c - r;
        e = 2 * r - c;
    }
    uecm_mfpt<MF> pt1, pt2, pt3;

    // the first one is always a doubling
    // point1 is [1]P
    pt1.X = pt2.X = pt3.X = P->X;
    pt1.Z = pt2.Z = pt3.Z = P->Z;

    // point2 is [2]P
    udup(mf, s, pt1, &pt1);

    while (d != e)
    {
        if (d < e)
        {
            uint64_t tmp = d;
            d = e;
            e = tmp;
            MV swpX = pt1.X;
            pt1.X = pt2.X;
            pt2.X = swpX;
            MV swpZ = pt1.Z;
            pt1.Z = pt2.Z;
            pt2.Z = swpZ;
        }

        if (d - e <= e / 4 && ((d + e) % 3) == 0)
        {
            d = (2 * d - e) / 3;
            e = (e - d) / 2;

            uecm_mfpt<MF> pt4;
            inlined_uadd(mf, pt1, pt2, pt3, &pt4); // T = A + B (C)
            uecm_mfpt<MF> pt5;
            inlined_uadd(mf, pt4, pt1, pt2, &pt5); // T2 = T + A (B)
            inlined_uadd(mf, pt2, pt4, pt1, &pt2); // B = B + T (A)

            pt1.X = pt5.X;
            pt1.Z = pt5.Z;
        }
        else if (d - e <= e / 4 && (d - e) % 6 == 0)
        {
            d = (d - e) / 2;
            inlined_uadd(mf, pt1, pt2, pt3, &pt2);        // B = A + B (C)
            inlined_udup(mf, s, pt1, &pt1);        // A = 2A
        }
        else if ((d + 3) / 4 <= e)
        {
            d -= e;

            uecm_mfpt<MF> pt4;
            inlined_uadd(mf, pt2, pt1, pt3, &pt4);        // T = B + A (C)

            pt3.X = pt2.X;
            pt2.X = pt4.X;
            pt3.Z = pt2.Z;
            pt2.Z = pt4.Z;
        }
        else if ((d + e) % 2 == 0)
        {
            d = (d - e) / 2;
            inlined_uadd(mf, pt2, pt1, pt3, &pt2);        // B = B + A (C)
            inlined_udup(mf, s, pt1, &pt1);        // A = 2A
        }
        else if (d % 2 == 0)
        {
            d /= 2;
            inlined_uadd(mf, pt3, pt1, pt2, &pt3);        // C = C + A (B)
            inlined_udup(mf, s, pt1, &pt1);        // A = 2A
        }
        else if (d % 3 == 0)
        {
            d = d / 3 - e;
            uecm_mfpt<MF> pt4;
            inlined_udup(mf, s, pt1, &pt4);        // T = 2A
            uecm_mfpt<MF> pt5;
            inlined_uadd(mf, pt1, pt2, pt3, &pt5);        // T2 = A + B (C)
            inlined_uadd(mf, pt4, pt1, pt1, &pt1);        // A = T + A (A)
            inlined_uadd(mf, pt4, pt5, pt3, &pt4);        // T = T + T2 (C)

            pt3.X = pt2.X;
            pt2.X = pt4.X;
            pt3.Z = pt2.Z;
            pt2.Z = pt4.Z;
        }
        else if ((d + e) % 3 == 0)
        {
            d = (d - 2 * e) / 3;
            uecm_mfpt<MF> pt4;
            inlined_uadd(mf, pt1, pt2, pt3, &pt4);        // T = A + B (C)
            inlined_uadd(mf, pt4, pt1, pt2, &pt2);        // B = T + A (B)
            inlined_udup(mf, s, pt1, &pt4);        // T = 2A
            inlined_uadd(mf, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else if ((d - e) % 3 == 0)
        {
            d = (d - e) / 3;
            uecm_mfpt<MF> pt4;
            inlined_uadd(mf, pt1, pt2, pt3, &pt4);        // T = A + B (C)
            inlined_uadd(mf, pt3, pt1, pt2, &pt3);        // C = C + A (B)

            pt2.X = pt4.X;
            pt2.Z = pt4.Z;

            inlined_udup(mf, s, pt1, &pt4);        // T = 2A
            inlined_uadd(mf, pt1, pt4, pt1, &pt1);        // A = A + T (A) = 3A
        }
        else
        {
            e /= 2;
            inlined_uadd(mf, pt3, pt2, pt1, &pt3);        // C = C + B (A)
            inlined_udup(mf, s, pt2, &pt2);        // B = 2B
        }
    }
    uadd(mf, pt1, pt2, pt3, P);     // A = A + B (C)

    for (int i = 0; i < shift; i++)
    {
        udup(mf, s, *P, P);     // P = 2P
    }
    return;
}



template <class MF>
static void uprac_precalc_80(const MF& mf, uecm_mfpt<MF> *P,
                             typename MF::MontgomeryValue s, int target_bits)
{
// For reference, here are some unused composites:
//        10,8,8,8,8,8,8,8,8,11,8,8,9,      // 37 * 53
//        11,8,11,8,10,8,8,8,8,8,8,8,8,9,   // 37 * 83
//        8,8,8,8,8,8,8,8,8,11,8,11,8,9,    // 41 * 53
// With the best choices I could make, I got less than 0.5% perf boost by using
// the best of these composites.  Since they help only slightly but increase
// complexity, I don't use them.

    int num_array_primes;
    if (target_bits <= 32)
        num_array_primes = 4;
    else if (target_bits <= 50)
        num_array_primes = target_bits - 29 + (target_bits <= 36);
    else
        num_array_primes = uecm_steps80<>::num_primes;


    uecm_mfpt<MF> pt1, pt2, pt3;

    pt2.X = pt3.X = P->X;
    pt2.Z = pt3.Z = P->Z;
    udup(mf, s, *P, &pt1);

    int primes_completed = 0;
    for (int i = 0; true; i++)
    {
        uint8_t step = uecm_steps80<>::steps[i];
        if (step == 8)
        {
            uecm_mfpt<MF> pt4;
            inlined_uadd(mf, pt2, pt1, pt3, &pt4);
            pt3.X = pt2.X;
            pt2.X = pt1.X;
            pt1.X = pt4.X;
            pt3.Z = pt2.Z;
            pt2.Z = pt1.Z;
            pt1.Z = pt4.Z;
        }
        else if (step == 9)
        {
            uecm_mfpt<MF> pt4;
            inlined_uadd(mf, pt2, pt1, pt3, &pt4);
            inlined_uadd(mf, pt1, pt4, pt2, P);

            primes_completed++;
            if (primes_completed == num_array_primes)
                break;

            pt2.X = pt3.X = P->X;
            pt2.Z = pt3.Z = P->Z;
            inlined_udup(mf, s, *P, &pt1);
        }
        else if (step == 10)
        {
            inlined_uadd(mf, pt2, pt1, pt3, &pt2);
            inlined_udup(mf, s, pt1, &pt1);
        }
        else // if (step == 11)
        {
            uecm_mfpt<MF> pt4;
            inlined_uadd(mf, pt2, pt1, pt3, &pt4);
            pt3.X = pt2.X;
            pt2.X = pt4.X;
            pt3.Z = pt2.Z;
            pt2.Z = pt4.Z;
        }
    }

    return;
}





template <class MF>
static void uecm_stage1(const MF& mf, uecm_mfpt<MF> *P,
                        typename MF::MontgomeryValue s, int target_bits)
{
    using T = typename MF::IntegerType;
    HPBC_PRECONDITION2(0 < target_bits && target_bits <= 128);

    // loops for the prime numbers 2 and 3
    int prime2_iterations = 7 + (target_bits > 41) + (target_bits > 48)
                            + (target_bits > 55);
    int prime3_iterations = 4 + (target_bits > 41) + (target_bits > 55);
    if constexpr (hurchalla::ut_numeric_limits<T>::digits > 64) {
        if (target_bits > 64) {
            prime2_iterations += ((target_bits - 64 + 7) >> 3);
            prime3_iterations += ((5*(target_bits - 64 + 11)) >> 6);
        }
    }
    for (int i = 0; i < prime2_iterations; i++) {
        inlined_udup(mf, s, *P, P);
    }
    for (int i = 0; i < prime3_iterations; i++) {
        uecm_mfpt<MF> tmp;
        udup(mf, s, *P, &tmp);
        uadd(mf, tmp, *P, *P, P);
    }

    uprac_precalc_80(mf, P, s, target_bits);

    // uprac_precalc_80 was by itself sufficient for bit widths up to 51
    if (target_bits <= 51)
        return;
    int limit;
    // We want 'limit' to be an ok fit to the settings that performed best
    // in tests:
    // target_bits   limit
    // -------------------
    // 52            2
    // 54            5
    // 56            8
    // 58            11
    // 60            16  15?
    // 62            20
    // 64            27
    // This is fairly linear up to 58.
    // We can add in a second linear fit from 60 to 63.
    if (target_bits >= 64) {
        limit = uecm_upracparams<>::num_params;
        static_assert(uecm_upracparams<>::num_params == 28);
    }
    else {
        HPBC_ASSERT2(target_bits >= 52);
        limit = ((5 + 3*(target_bits - 52)) >> 1)
                + (target_bits >= 60)*(target_bits - 59);
        HPBC_ASSERT2(limit <= uecm_upracparams<>::num_params);
    }
    for (int i=0; i<limit; ++i)
    {
        uprac(mf, P, uecm_upracparams<>::params[i].prime,
              uecm_upracparams<>::params[i].ratio, s);
    }

    if constexpr (hurchalla::ut_numeric_limits<T>::digits > 64) {
        // the uprac loop above was sufficient for bit widths up to 64
        if (target_bits <= 64)
            return;
/*
        Set limit2 to values that performed best (roughly measured) in tests...
        Roughly best performing values for limit2 were
        66  bits: 7
        74  bits: 37
        82  bits: 79
        90  bits: 139
        100 bits: 242
        108 bits: 346
        116 bits: 494
        124 bits: 710
        The final point (129 bits, 880)  is a rough extrapolation.
        An exponential function is probably best fit to data, but linear
        interpolation is fine and easier.
*/
        struct BitPoint {
            int bits;
            int indexlimit;
        };
        constexpr std::array<BitPoint, 10> bitpoints = {{
                      {64, 0}, {66, 7}, {74, 37}, {82, 79}, {90, 139},
                      {100, 242}, {108, 346}, {116, 494}, {124, 710}, {129, 881}
                  }};
        float slope = 0.0f;
        float intercept = 0.0f;
        for (decltype(bitpoints.size()) i = 1; i < bitpoints.size(); ++i) {
            if (target_bits < bitpoints[i].bits) {
                slope = (static_cast<float>(bitpoints[i].indexlimit)
                           - static_cast<float>(bitpoints[i-1].indexlimit))
                        / (static_cast<float>(bitpoints[i].bits)
                           - static_cast<float>(bitpoints[i-1].bits));
                intercept = static_cast<float>(bitpoints[i-1].indexlimit)
                       - slope*(static_cast<float>(bitpoints[i-1].bits) - 0.4f);
                break;
            }
        }
        int limit2 = static_cast<int>(slope*target_bits + intercept + 0.5f);
        if (limit2 > uecm_upracparams2<>::num_params)
            limit2 = uecm_upracparams2<>::num_params;

        for (int i = 0; i < limit2; ++i)
            uprac(mf, P, uecm_upracparams2<>::primes[i], 0.54, s);
    }
    return;
}




template <int bitsT>
static int get_stage2_num_giant_steps(int target_bits)
{
    // These settings for num_giant_steps were empirically determined.  We use
    // simple linear interpolation between best performing values found in perf
    // testing.
    int num_giant_steps;
    if (target_bits < 38)
        num_giant_steps = 4;   // 3 steps benchmarked worse until bits <= 32
    else if (target_bits < 53)
        num_giant_steps = 4 + (target_bits - 38) + (target_bits==38);
    else if (target_bits < 64) {
        num_giant_steps = 20 + (7*(target_bits - 53) >> 2)
                          + 2*(target_bits>61)*(target_bits - 61);
    }
    else {
      if constexpr (bitsT <= 64)
        num_giant_steps = 44;
      else {
/*
        Roughly best performing values for num_giant_steps were
        64  bits: 43
        66  bits: 48
        74  bits: 78
        82  bits: 116
        90  bits: 186
        100 bits: 273
        108 bits: 388
        116 bits: 577
        124 bits: 897
        The final point (129 bits, 1164)  is a rough extrapolation.
        An exponential function is probably best fit to data, but linear
        interpolation is fine and easier.
*/
        struct BitPoint2 {
            int bits;
            int numsteps;
        };
        constexpr std::array<BitPoint2, 10> bitpoints = {{
                     {64, 43}, {66, 48}, {74, 78}, {82, 116}, {90, 186},
                     {100, 273}, {108, 388}, {116, 577}, {124, 897}, {129, 1164}
                  }};
        float slope = 0.0f;
        float intercept = 0.0f;
        for (decltype(bitpoints.size()) i = 1; i < bitpoints.size(); ++i) {
            if (target_bits < bitpoints[i].bits) {
                slope = (static_cast<float>(bitpoints[i].numsteps)
                           - static_cast<float>(bitpoints[i-1].numsteps))
                        / (static_cast<float>(bitpoints[i].bits)
                           - static_cast<float>(bitpoints[i-1].bits));
                intercept = static_cast<float>(bitpoints[i-1].numsteps)
                       - slope*(static_cast<float>(bitpoints[i-1].bits) - 0.4f);
                break;
            }
        }
        num_giant_steps = static_cast<int>(slope*target_bits + intercept + 0.5f);
        HPBC_ASSERT2(num_giant_steps > 0);
      }
    }
    return num_giant_steps;
}




#ifdef HURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM
template <class T>
struct trait_is_dual_monty {
    static constexpr bool value = false;
};
template <class T>
struct trait_is_dual_monty<hurchalla::DualMontgomeryForm<T>> {
    static constexpr bool value = true;
};
template <class MF>
inline constexpr bool is_dual_monty = trait_is_dual_monty<MF>::value;
#endif




template <class MF>
static
typename MF::MontgomeryValue uecm_stage2(const MF& mf, const uecm_mfpt<MF>& P,
                                int target_bits, typename MF::MontgomeryValue s)
{
    static constexpr uint32_t map[61] = {
        0, 1, 2, 0, 0, 0, 0, 3, 0, 0,
        0, 4, 0, 5, 0, 0, 0, 6, 0, 7,
        0, 0, 0, 8, 0, 0, 0, 0, 0, 9,
        0, 10, 0, 0, 0, 0, 0, 11, 0, 0,
        0, 12, 0, 13, 0, 0, 0, 14, 0, 15,
        0, 0, 0, 16, 0, 0, 0, 0, 0, 17,
        18 };

    using MV = typename MF::MontgomeryValue;
    using CV = typename MF::CanonicalValue;

    constexpr int ECM_PARAM_D = 60;

    //stage 2 init
    //Q = P = result of stage 1
    //compute [d]Q for 0 < d <= D
    uecm_mfpt<MF> Pb[20];
    uecm_mfpt<MF> *Pd = &Pb[map[ECM_PARAM_D]];

    CV Pbprod[20];

    // [1]Q
    Pb[1].Z = P.Z;
    Pb[1].X = P.X;
    Pbprod[1] = mf.getCanonicalValue(mf.multiply(Pb[1].X, Pb[1].Z));

    // [2]Q
    udup(mf, s, P, &Pb[2]);
    Pbprod[2] = mf.getCanonicalValue(mf.multiply(Pb[2].X, Pb[2].Z));

    /*
    D is small in tinyecm, so it is straightforward to just enumerate the needed
    points.  We can do it efficiently with two progressions mod 6.
    Pb[0] = scratch
    Pb[1] = [1]Q;
    Pb[2] = [2]Q;
    Pb[3] = [7]Q;   prog2
    Pb[4] = [11]Q;  prog1
    Pb[5] = [13]Q;  prog2
    Pb[6] = [17]Q;  prog1
    Pb[7] = [19]Q;  prog2
    Pb[8] = [23]Q;  prog1
    Pb[9] = [29]Q;  prog1
    Pb[10] = [30 == D]Q;
    Pb[11] = [31]Q; prog2
    Pb[12] = [37]Q; prog2
    Pb[13] = [41]Q; prog1
    Pb[14] = [43]Q; prog2
    Pb[15] = [47]Q; prog1
    Pb[16] = [49]Q; prog2
    Pb[17] = [53]Q; prog1
    Pb[18] = [59]Q; prog1

    two progressions with total of 17 adds to get 15 values of Pb.
    6 + 5(1) -> 11 + 6(5) -> 17 + 6(11) -> 23 + 6(17) -> 29 + 6(23) -> 35 + 6(29) -> 41 + 6(35) -> 47 + 6(41) -> 53 + 6(47) -> 59
    6 + 1(5) -> 7 + 6(1) -> 13 + 6(7) -> 19 + 6(13) -> 25 + 6(19) -> 31 + 6(25) -> 37 + 6(31) -> 43 + 6(37) -> 49

    we also need [2D]Q = [60]Q
    to get [60]Q we just need one more add:
    compute [60]Q from [31]Q + [29]Q([2]Q), all of which we
    have after the above progressions are computed.

    we also need [A]Q = [((B1 + D) / (2D) * 2D]Q
    which is equal to the following for various common B1:
    B1      [x]Q
    65      [120]Q
    85      [120]Q
    125     [180]Q
    165     [180]Q
    205     [240]Q

    and we need [x-D]Q as well, for the above [x]Q.
    So far we are getting [x]Q and [x-2D]Q each from prac(x,Q).
    There is a better way using progressions of [2D]Q
    [120]Q = 2*[60]Q
    [180]Q = [120]Q + [60]Q([60]Q)
    [240]Q = 2*[120]Q
    [300]Q = [240]Q + [60]Q([180]Q)
    ...
    etc.
    */

    uecm_mfpt<MF> pt1, pt3;

    // Calculate all Pb: the following is specialized for D=60
    // [2]Q + [1]Q([1]Q) = [3]Q
    uadd(mf, Pb[1], Pb[2], Pb[1], &Pb[3]);        // <-- temporary

    // 2*[3]Q = [6]Q
    udup(mf, s, Pb[3], &pt3);   // pt3 = [6]Q

    // [3]Q + [2]Q([1]Q) = [5]Q
    uadd(mf, Pb[3], Pb[2], Pb[1], &pt1);    // <-- pt1 = [5]Q
    Pb[3].X = pt1.X;
    Pb[3].Z = pt1.Z;

    // [6]Q + [5]Q([1]Q) = [11]Q
    uadd(mf, pt3, pt1, Pb[1], &Pb[4]);    // <-- [11]Q

    int h = 3;
    int k = 4;
    int j = 5;
    while ((j + 12) < ECM_PARAM_D)
    {
        // [j+6]Q + [6]Q([j]Q) = [j+12]Q
        uadd(mf, pt3, Pb[k], Pb[h], &Pb[map[j + 12]]);
        h = k;
        k = map[j + 12];
        j += 6;
    }

    // [6]Q + [1]Q([5]Q) = [7]Q
    uadd(mf, pt3, Pb[1], pt1, &Pb[3]);    // <-- [7]Q
    h = 1;
    k = 3;
    j = 1;
    while ((j + 12) < ECM_PARAM_D)
    {
        // [j+6]Q + [6]Q([j]Q) = [j+12]Q
        uadd(mf, pt3, Pb[k], Pb[h], &Pb[map[j + 12]]);
        h = k;
        k = map[j + 12];
        j += 6;
    }

    // Pd = [2w]Q
    // [31]Q + [29]Q([2]Q) = [60]Q
    uadd(mf, Pb[9], Pb[10], Pb[2], Pd);   // <-- [60]Q

    // make all of the Pbprod's
    for (int i = 3; i < 19; i++)
    {
        Pbprod[i] = mf.getCanonicalValue(mf.multiply(Pb[i].X, Pb[i].Z));
    }


    //initialize info needed for giant step
    // temporary - make [4]Q
    udup(mf, s, Pb[2], &pt3);   // pt3 = [4]Q

    uecm_mfpt<MF> Pad;

    // Pd = [w]Q
    // [17]Q + [13]Q([4]Q) = [30]Q
    uadd(mf, Pb[map[17]], Pb[map[13]], pt3, &Pad);    // <-- [30]Q


    uecm_mfpt<MF> Pa1;
    uecm_mfpt<MF> *Pa = &Pa1;

    // [60]Q + [30]Q([30]Q) = [90]Q
    uadd(mf, *Pd, Pad, Pad, Pa);
    pt1.X = Pa->X;
    pt1.Z = Pa->Z;

    // [90]Q + [30]Q([60]Q) = [120]Q
    uadd(mf, *Pa, Pad, *Pd, Pa);
    Pd->X = Pa->X;
    Pd->Z = Pa->Z;

    // [120]Q + [30]Q([90]Q) = [150]Q
    uadd(mf, *Pa, Pad, pt1, Pa);


    // adjustment of Pa and Pad for larger B1.
    // Currently we have Pa=150, Pd=120, Pad=30
    if (target_bits > 58)
    {
        if (target_bits <= 62)
        {
            // need Pa = 180, Pad = 60
            // [150]Q + [30]Q([120]Q) = [180]Q
            uadd(mf, *Pa, Pad, *Pd, Pa);
            udup(mf, s, Pad, &Pad);   // Pad = [60]Q
        }
        else  // if (target_bits > 62)
        {
            // need Pa = 210, Pad = 90.
            // have pt1 = 90
            udup(mf, s, Pad, &Pad);   // Pad = [60]Q
            // [150]Q + [60]Q([90]Q) = [210]Q
            uadd(mf, *Pa, Pad, pt1, Pa);
            Pad.X = pt1.X;
            Pad.Z = pt1.Z;
        }
//TODO: We haven't adjusted for target_bits > 64 yet.  Presumably we should...
//TODO: Is the earlier setup in this function appropriate for target_bits > 64?
    }

    //initialize accumulator and Paprod
    MV acc = mf.getUnityValue();
    MV Paprod = mf.multiply(Pa->X, Pa->Z);


    using T = typename MF::IntegerType;
    int num_giant_steps = get_stage2_num_giant_steps<
                          hurchalla::ut_numeric_limits<T>::digits>(target_bits);

    MV acc2 = mf.getUnityValue();

    uecm_mfpt<MF> Pb2[16];
    CV Pbprod2[16];
    Pb2[0] = Pb[1];
    Pbprod2[0] = Pbprod[1];
    for (int i = 1; i < 16; ++i) {
        Pb2[i] = Pb[i+2];
        Pbprod2[i] = Pbprod[i+2];
    }

    int g = 0;
    while (true)
    {
        for (int i = 0; i < 16; ++i)
        {
            //we accumulate XrZd - XdZr = (Xr - Xd) * (Zr + Zd) + XdZd - XrZr
            //in CP notation, Pa -> (Xr,Zr), Pb -> (Xd,Zd)
            // accumulate the cross product  (zimmerman syntax).
            // page 342 in C&P
            MV tt1 = mf.subtract(Pa->X, Pb2[i].X);
            MV tt2 = mf.add(Pa->Z, Pb2[i].Z);
            tt1 = mf.fmadd(tt1, tt2, Pbprod2[i]);
            tt2 = mf.subtract(tt1, Paprod);
            bool isZero;
            MV tmp = mf.multiply(acc, tt2, isZero);
            if (isZero)
                break;
            acc = tmp;
# ifdef HURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM
            // the dual_monty form already provides us with instruction level
            // parallelism in this loop.  We don't need to create more ILP via
            // a separate accumulator in the loop unroll below.
            if constexpr (!is_dual_monty<MF>)
# endif
            {
                ++i;
                MV tt3 = mf.subtract(Pa->X, Pb2[i].X);
                MV tt4 = mf.add(Pa->Z, Pb2[i].Z);
                tt3 = mf.fmadd(tt3, tt4, Pbprod2[i]);
                tt4 = mf.subtract(tt3, Paprod);
                bool isZero2;
                MV tmp2 = mf.multiply(acc2, tt4, isZero2);
                if (isZero2)
                    break;
                acc2 = tmp2;
            }
        }
        if (g == num_giant_steps)
            break;
        ++g;

        //giant step - use the addition formula for ECM
        pt1.X = Pa->X;
        pt1.Z = Pa->Z;
        //Pa + Pd
        inlined_uadd(mf, *Pa, *Pd, Pad, Pa);
        //Pad holds the previous Pa
        Pad.X = pt1.X;
        Pad.Z = pt1.Z;
        //and Paprod
        Paprod = mf.multiply(Pa->X, Pa->Z);
    }

# ifdef HURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM
    if constexpr (!is_dual_monty<MF>)
# endif
    {
        bool isZero;
        MV tmp = mf.multiply(acc, acc2, isZero);
        if (!isZero)
            acc = tmp;
    }
    return acc;
}





template <typename T>
static inline T modinverse(T a, T modulus, T& gcd)
{
    return hurchalla::modular_multiplicative_inverse(a, modulus, gcd);
}


#if !defined(HURCHALLA_TARGET_CPU_HAS_FAST_DIVIDE)

// Overload for uint64_t, for performance.
// jeff: "likely_gcd" is probably always the correct gcd, but I have not made
// any proof to show that it is indeed always the gcd.  See comments within this
// function for more info.  It's possible that likely_gcd might be wrong
// sometime despite passing all my tests so far.
static inline uint64_t modinverse(uint64_t a, uint64_t p, uint64_t& likely_gcd)
{
    /* thanks to the folks at www.mersenneforum.org */

    uint64_t ps1, ps2, parity, dividend, divisor, rem, q, t;

    q = 1;
    rem = a;
    dividend = p;
    divisor = a;
    ps1 = 1;
    ps2 = 0;
    parity = 0;

    while (divisor > 1) {
        rem = dividend - divisor;
        t = rem - divisor;
        if (rem >= divisor) {
            q += ps1; rem = t; t -= divisor;
            if (rem >= divisor) {
                q += ps1; rem = t; t -= divisor;
                if (rem >= divisor) {
                    q += ps1; rem = t; t -= divisor;
                    if (rem >= divisor) {
                        q += ps1; rem = t; t -= divisor;
                        if (rem >= divisor) {
                            q += ps1; rem = t; t -= divisor;
                            if (rem >= divisor) {
                                q += ps1; rem = t; t -= divisor;
                                if (rem >= divisor) {
                                    q += ps1; rem = t; t -= divisor;
                                    if (rem >= divisor) {
                                        q += ps1; rem = t;
                                        if (rem >= divisor) {
                                            q = dividend / divisor;
                                            rem = dividend % divisor;
                                            q *= ps1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        q += ps2;
        parity = ~parity;
        dividend = divisor;
        divisor = rem;
        ps2 = ps1;
        ps1 = q;
    }

    // jeff: added "likely_gcd".  this function seems to be a variant on the
    // extended euclidean algorithm, and thus the gcd likely equals the dividend
    // as calculated below.  However I'm doing this by analogy with the extended
    // euclidean and by educated guess, not by proof.  It appears to work in all
    // tests, so I suspect it is correct, but it's possible this could be wrong
    // in some cases.
    if (divisor == 1)
        dividend = divisor;
    likely_gcd = dividend;

    if (parity == 0)
        return ps1;
    else
        return p - ps1;
}

#endif




template <class MF>
static typename MF::IntegerType
ubuild(const MF& HURCHALLA_RESTRICT mf,
       uecm_mfpt<MF>& HURCHALLA_RESTRICT P,
       uint64_t& HURCHALLA_RESTRICT loc_lcg,
       typename MF::MontgomeryValue& HURCHALLA_RESTRICT s)
{
    using MV = typename MF::MontgomeryValue;
    using CV = typename MF::CanonicalValue;
    using T = typename MF::IntegerType;

    T n = mf.getModulus();
    uint32_t sigma = lcg_rand_32B(7, (uint32_t)-1, &loc_lcg);
    MV u;

    // failing the static assert might not be a problem, but it's unexpected
    // given there's likely no reason to use ecm with such small types, so I
    // haven't checked to see if it would be ok (or not).
    static_assert(hurchalla::ut_numeric_limits<T>::digits >= 32);
#ifdef HURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM
    if constexpr (is_dual_monty<MF>) {
        uint32_t sigma2 = lcg_rand_32B(7, (uint32_t)-1, &loc_lcg);
        u = mf.convertIn(static_cast<T>(sigma), static_cast<T>(sigma2));
    }
    else
#endif
        u = mf.convertIn(static_cast<T>(sigma));

    //printf("sigma = %" PRIu64 ", u = %" PRIu64 ", n = %" PRIu64 "\n", sigma, u, n);

    CV cu = mf.getCanonicalValue(u);
    CV cv = mf.add(cu, cu);
    cv = mf.add(cv, cv);                 // v = 4*sigma

    CV one = mf.getUnityValue();
    CV two = mf.add(one, one);
    CV four = mf.add(two, two);
    CV five = mf.add(four, one);

    u = mf.fusedSquareSub(cu, five);         // u = sigma^2 - 5
    cu = mf.getCanonicalValue(u);

    MV mvx = mf.multiply(mf.square(u), u);     // x = u^3

    CV cv2 = mf.add(cv, cv);             // 2*v
    CV cv4 = mf.add(cv2, cv2);           // 4*v
    CV cv8 = mf.add(cv4, cv4);           // 8*v
    CV cv16 = mf.add(cv8, cv8);          // 16*v
    MV t5 = mf.multiply(cv16, mvx);         // 16*v*u^3

    MV mvz = mf.multiply(mf.square(cv), cv);   // z = v^3

    //compute parameter A
    MV t1 = mf.subtract(cv, cu);
    MV t4 = mf.multiply(mf.square(t1), t1);  // (v - u)^3

    CV t7 = mf.add(cu, cu);
    CV t8 = mf.add(cu, cv);
    MV t3 = mf.add(t7, t8);           // 3u + v

    t1 = mf.multiply(t3, t4);         // t1 = (v-u)^3 * (3u + v)

    T likely_gcd = 1;
#ifdef HURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM
    if constexpr (is_dual_monty<MF>) {
#  if 0
  // We normally disable this section since the #else should be roughly twice
  // as fast.  However, this section more clearly shows what we want to do.
        T s4a, s4b;
        mf.convertOut(s4a, s4b, t5);
        T likely_gcdA;
        T s3a = modinverse(s4a, n, likely_gcdA);
        T s3b = modinverse(s4b, n, likely_gcd);
        if (likely_gcdA > 1)
            likely_gcd = likely_gcdA;
        t3 = mf.convertIn(s3a, s3b);
#  else
        // For the following description, let s4a and s4b be the values that
        // would be obtained by calling
        //   mf.convertOut(s4a, s4b, t5);
        // We will use the facts that
        //    inv(s4a) == inv(s4a * s4b) * s4b  (mod n).
        //    inv(s4b) == inv(s4a * s4b) * s4a  (mod n).
        // For details, see section 2.2 of "Application of Montgomery's Trick to Scalar Multiplication for Elliptic and Hyperelliptic Curves Using a Fixed Base Point"
        //     https://www.iacr.org/archive/pkc2004/29470042/29470042.pdf
        //   See also https://twitter.com/vitalikbuterin/status/1246213886338048000?lang=en
        //   And the original paper by Montgomery, "Speeding the Pollard and Elliptic Curve Methods of Factorization"
        //   https://wstein.org/edu/Fall2001/124/misc/montgomery.pdf

        // The next line sets s4ab = (s4a * s4b) mod n.
        T s4ab = mf.crossMultiplyAndConvertOut(t5);

        // get the inverse of (s4a * s4b) mod n.
        T inv_ab = modinverse(s4ab, n, likely_gcd);

        // set t6 = mf.convertIn(inv_ab, inv_ab)
        MV t6 = mf.convertInAndCopy(inv_ab);

        // set t5swap = mf.convertIn(s4b, s4a).
        MV t5swap = mf.swapChannels(t5);

        // Thus the following multiply gives us the equivalent of
        // t3 = mf.convertIn(modinverse(s4a, n), modinverse(s4b, n))
        t3 = mf.multiply(t6, t5swap);
#  endif
    }
    else
#endif
    {
        T s4 = mf.convertOut(t5);
        T s3 = modinverse(s4, n, likely_gcd);
        t3 = mf.convertIn(s3);
    }
    // accomplish the division by multiplying by the modular inverse
    MV mvs = mf.multiply(t3, t1);

    P.X = mvx;
    P.Z = mvz;
    s = mvs;

    return likely_gcd;
}





template <class MF>
static bool ucheck_factor(const MF& HURCHALLA_RESTRICT mf,
                          typename MF::MontgomeryValue z,
                          typename MF::IntegerType n,
                          typename MF::IntegerType& HURCHALLA_RESTRICT f)
{
//    f = hurchalla::greatest_common_divisor(mf.convertOut(z), n);
    auto gcdfunc = [](auto x, auto y)
                   { return hurchalla::greatest_common_divisor(x, y); };

    bool found_factor = false;
    f = mf.gcd_with_modulus(z, gcdfunc);
    if (f > 1) {
        if (f == n) {
            f = 0;
            found_factor = false;
        }
        else
            found_factor = true;
    }
    return found_factor;
}



template <class MF>
static bool microecm(const MF& HURCHALLA_RESTRICT mf,
                     typename MF::IntegerType& HURCHALLA_RESTRICT f,
                     int curves,
                     uint64_t& HURCHALLA_RESTRICT loc_lcg,
                     int target_bits)
{
    using MV = typename MF::MontgomeryValue;
    using T = typename MF::IntegerType;
    T n = mf.getModulus();
/*
    uint32_t B1;
    if (target_bits <= 42)
        B1 = 35;
    else if (target_bits <= 49) // was once 48 I believe
        B1 = 70;
    else if (target_bits <= 52)
        B1 = 85;
    else if (target_bits <= 58)
        B1 = 125;
    else if (target_bits <= 62)
        B1 = 165;
    else if (target_bits <= 64)
        B1 = 205;
    uint32_t B2 = 25*B1;
    uint32_t stg1_max = B1;
*/

    //attempt to factor n with the elliptic curve method
    //following brent and montgomery's papers, and CP's book

    bool found_factor = false;

    int curve_advance = 1;
#ifdef HURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM
    if constexpr (is_dual_monty<MF>)
        curve_advance = 2;
#endif

    f = 1;
    for (int i = 0; i < curves; i += curve_advance)
    {
        MV s;
        uecm_mfpt<MF> P;
        T likely_gcd = ubuild(mf, P, loc_lcg, s);
        if (likely_gcd > 1)
        {
            // If the gcd gave us a factor, we're done.  If not, since gcd != 1,
            // the inverse calculated in ubuild would be bogus and so this curve
            // is probably set up for failure (hence we continue).
            if (likely_gcd != n && n % likely_gcd == 0) {
                f = likely_gcd;
                return true;
            }
            continue;
        }

        uecm_stage1(mf, &P, s, target_bits);
        found_factor = ucheck_factor(mf, P.Z, n, f);
        if (found_factor)
            return true;

        MV stg2acc = uecm_stage2(mf, P, target_bits, s);
        found_factor = ucheck_factor(mf, stg2acc, n, f);
        if (found_factor)
            return true;
    }

    return found_factor;
}


template <typename T>
static int ecm_getbits(T n)
{
    static_assert(hurchalla::ut_numeric_limits<T>::is_integer, "");
    int i = 0;
    while (n != 0)
    {
        n >>= 1;
        i++;
    }
    return i;
}


// Prior to your first call of get_ecm_factor(), set loc_lcg = 0 (or some
// arbitrary value).  After that, don't change loc_lcg.
// FYI: loc_lcg is used within this file by a random number generator, and
// holds the current value of a pseudo random sequence.  Your first assigment
// to loc_lcg seeds the sequence, and after seeding it you don't want to
// change loc_lcg, since that would interrupt the sequence.
template <class MF>
static typename MF::IntegerType get_ecm_factor(const MF& HURCHALLA_RESTRICT mf,
                                           uint64_t& HURCHALLA_RESTRICT loc_lcg)
{
    using T = typename MF::IntegerType;
    int curves;
    T f64 = 1;
    bool found_factor = false;

    int targetBits = ecm_getbits(mf.getModulus());

#ifndef HURCHALLA_FACTORING_EXPECT_LARGE_FACTORS
    // since we don't expect large factors, try fast attempts to find
    // potential small factors.
    curves = 1;
    if constexpr (hurchalla::ut_numeric_limits<T>::digits <= 64) {
        if (targetBits > 39) {
            int tmp_bits = 39;
            found_factor = microecm(mf, f64, curves, loc_lcg, tmp_bits);
            if (found_factor)
                return f64;
            if (targetBits > 45) {
                tmp_bits = 45;
                found_factor = microecm(mf, f64, curves, loc_lcg, tmp_bits);
                if (found_factor)
                    return f64;

                if (targetBits > 51) {
                    tmp_bits = 51;
                    found_factor = microecm(mf, f64, curves, loc_lcg, tmp_bits);
                    if (found_factor)
                        return f64;

                    if (targetBits > 58) {
                        tmp_bits = 58;
                        found_factor = microecm(mf, f64, curves, loc_lcg, tmp_bits);
                        if (found_factor)
                            return f64;
                    }
                }
            }
        }
    }
    else {
        for (int tmp_bits = 34; targetBits > tmp_bits; tmp_bits += 6) {
            found_factor = microecm(mf, f64, curves, loc_lcg, tmp_bits);
            if (found_factor)
                return f64;
        }
    }
#endif

    curves = 16*targetBits;

#ifndef HURCHALLA_ECM_ALLOW_DUAL_MONTGOMERY_FORM
    found_factor = microecm(mf, f64, curves, loc_lcg, targetBits);
#else
    if (targetBits <= 62)
        found_factor = microecm(mf, f64, curves, loc_lcg, targetBits);
    else
    {
        if constexpr (hurchalla::ut_numeric_limits<T>::digits <=
                      HURCHALLA_TARGET_BIT_WIDTH)
        {
            hurchalla::DualMontgomeryForm<MF> dmf(mf);
            found_factor = microecm(dmf, f64, curves, loc_lcg, targetBits);
        }
        else
            found_factor = microecm(mf, f64, curves, loc_lcg, targetBits);
    }
#endif

    return f64;
}


}; // end struct micro_ecm



}} // end namespace


#endif // include guard
