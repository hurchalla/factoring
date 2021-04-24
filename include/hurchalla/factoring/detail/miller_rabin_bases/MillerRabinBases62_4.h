// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES62_4_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES62_4_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 1472 byte hash table lets you determine the primality of any unsigned
// int number less than (1<<62), via miller-rabin primality testing using 4
// bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<62, 4, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 64 bit number being tested for primality.
    // As a precondition, 'num' must be less than (1<<62).
    // This function returns 4 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 4> get(std::uint64_t num)
    {
        HPBC_PRECONDITION2(num < (static_cast<std::uint64_t>(1) << 62));
        // I generated/verified the hash table and bases.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 4> bases;
        bases[0] = 2;
        bases[1] = 15;
        uint32_t mask = (static_cast<uint32_t>(1) << 22) - 1u;
        uint32_t hash_bucket = ((static_cast<uint32_t>(num) & mask) * 23) >> 18;
        bases[2] = table[hash_bucket][0];
        bases[3] = table[hash_bucket][1];
        return bases;
    }
private:
    static constexpr std::size_t SIZE = 368;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE][2] = {
#else
    static constexpr std::uint16_t table[][2] = {
#endif
        { 32269, 37115 },
        { 64613, 62921 },
        { 43367, 62829 },
        { 31840, 46195 },
        { 32625, 19817 },
        { 21728, 36728 },
        { 36203, 11681 },
        { 62413, 43389 },
        { 57427, 57223 },
        { 54978, 23082 },
        { 51301, 52403 },
        {  1053, 54158 },
        { 36403, 16491 },
        { 63518, 18312 },
        { 21312, 27582 },
        { 35263, 40467 },
        { 43983, 25581 },
        { 44225, 64985 },
        { 61329, 13893 },
        {  2493, 23623 },
        { 50377, 36542 },
        { 54925, 58567 },
        { 64055, 27643 },
        { 43940, 31284 },
        { 60214, 48665 },
        { 52631, 61663 },
        { 65247, 56253 },
        { 24657, 16679 },
        {  1714, 21392 },
        { 43986,  4705 },
        { 47598, 56232 },
        { 40483, 43790 },
        { 28431, 47009 },
        { 34558,  9306 },
        { 19294, 36849 },
        { 58031, 55609 },
        { 40849, 32831 },
        {  9763, 25152 },
        { 12574, 57909 },
        { 52014, 23236 },
        { 17965, 24983 },
        { 26473, 43481 },
        {  1427, 22357 },
        { 62568, 18626 },
        { 60260, 23779 },
        { 32470, 62077 },
        { 61723, 58455 },
        { 60408, 45428 },
        { 51573, 38679 },
        { 38859, 63761 },
        { 32845, 32121 },
        { 14595, 39422 },
        { 18785, 37821 },
        {  8125, 11701 },
        { 24229, 65103 },
        { 16858, 45384 },
        { 49787, 46583 },
        { 65533, 57217 },
        { 45861, 43971 },
        { 51377, 58887 },
        { 32929, 58213 },
        {  6877, 50869 },
        { 36331, 27610 },
        {  6978, 43025 },
        { 31969, 30358 },
        {  1755, 63912 },
        { 17797, 48616 },
        { 16620, 59678 },
        { 39359,  7228 },
        { 34986,  6290 },
        { 54657,  5373 },
        { 32583, 25791 },
        { 37664, 23801 },
        { 40535, 10127 },
        { 30078, 28410 },
        {  3757, 25931 },
        { 26113, 58898 },
        { 57647,  8177 },
        {  6877, 10556 },
        { 21863, 30415 },
        {   915, 61531 },
        {  3625, 65073 },
        { 35809, 50790 },
        { 28962, 30741 },
        { 50951, 33681 },
        { 26454, 24750 },
        { 62234, 42860 },
        { 27195, 56889 },
        { 31853, 22461 },
        { 56709,  6733 },
        { 43313, 23458 },
        { 60298, 15037 },
        { 29952, 34354 },
        { 20109, 52060 },
        { 47473, 29365 },
        { 14905, 38041 },
        { 21853, 51485 },
        { 64566, 47435 },
        { 54374, 16830 },
        { 51475, 48607 },
        { 45449, 49631 },
        { 44298, 46053 },
        { 54730, 39905 },
        { 31478, 24714 },
        { 10007, 51398 },
        { 12835, 34476 },
        { 51305, 45170 },
        { 35964, 62121 },
        { 21737,  7660 },
        {  1804, 47819 },
        { 59143, 18509 },
        { 50519, 32134 },
        { 38125, 55445 },
        { 65533,  1569 },
        { 61396,  2857 },
        { 41289,  9649 },
        { 38556, 45306 },
        { 44955, 52977 },
        { 60510, 28189 },
        { 28676, 22279 },
        { 16582, 54998 },
        { 41734, 37614 },
        {  4284, 28007 },
        {  5265, 34925 },
        { 60214, 62181 },
        { 19041, 35787 },
        { 29064, 33818 },
        { 45096, 44278 },
        { 50522, 38228 },
        { 34420, 57081 },
        { 27133, 49487 },
        { 28962, 48577 },
        { 13431, 49367 },
        { 22337,  5236 },
        { 29835, 48482 },
        { 26426, 49355 },
        { 12137, 40733 },
        { 31213, 32238 },
        { 14376, 36347 },
        {  6881, 55438 },
        {  5950, 52869 },
        { 15925, 18682 },
        { 30686, 12248 },
        { 63289,  8507 },
        { 58231,  4855 },
        { 32612, 17822 },
        { 62170, 32498 },
        { 28120, 36561 },
        { 44054, 22742 },
        { 58357, 23817 },
        { 58603, 50346 },
        { 32012, 63266 },
        { 52771, 52019 },
        {  5200, 14804 },
        { 39546, 26659 },
        { 14940, 19825 },
        { 42959, 27348 },
        { 30169, 47709 },
        { 49937, 22263 },
        { 54710, 37090 },
        { 61893, 49869 },
        { 24705, 34785 },
        { 39444, 34939 },
        { 55500, 57099 },
        { 65533, 40370 },
        { 45889, 16918 },
        {  6830, 27517 },
        { 31943, 62957 },
        { 36026, 29560 },
        { 31727, 46905 },
        {  4693, 16529 },
        { 59851, 47635 },
        { 32079, 23581 },
        { 53946, 47873 },
        { 11640, 58395 },
        { 16419, 59105 },
        { 34313, 55133 },
        { 37518, 42883 },
        { 16470, 53588 },
        { 44668, 32249 },
        { 62426, 42185 },
        { 58154, 47042 },
        { 36260, 59642 },
        { 32625, 31943 },
        { 44080, 42914 },
        { 53805, 35555 },
        {  2925, 12552 },
        { 46284, 55738 },
        { 53391, 55172 },
        { 42237, 44866 },
        { 35546, 58562 },
        { 32079, 43965 },
        { 29952, 49341 },
        { 43486, 27963 },
        { 63942, 44523 },
        { 19573, 62347 },
        {  6207, 17662 },
        { 23146, 39387 },
        { 56355, 26318 },
        { 18785, 53333 },
        { 16919, 52015 },
        { 54925, 64341 },
        {  4693,  4341 },
        { 13702, 63323 },
        { 56709,  5050 },
        { 21970, 48063 },
        { 42237, 36321 },
        {  6292, 18733 },
        { 21612, 51915 },
        { 31265, 55315 },
        { 52670, 28421 },
        { 33022, 27179 },
        { 39239, 53394 },
        { 52345, 43505 },
        { 34385, 42553 },
        { 35016, 53558 },
        { 28385,  9502 },
        { 54925, 26216 },
        { 26210, 64470 },
        { 53833, 20455 },
        { 46904, 21469 },
        { 12029, 51282 },
        { 37583, 22897 },
        { 36371, 44698 },
        { 41242,  7565 },
        { 45325, 65509 },
        { 33901, 41521 },
        { 46133, 42428 },
        { 12636,   514 },
        { 45287, 45661 },
        { 57202, 39890 },
        { 18380, 49158 },
        { 28935, 39396 },
        { 41111, 49174 },
        { 59292,   881 },
        { 57434, 48535 },
        {  6702, 16518 },
        { 40875, 38123 },
        { 54282, 41667 },
        { 11837, 48726 },
        { 35964, 14357 },
        { 63518, 39385 },
        {  7857, 26525 },
        { 62465, 46462 },
        { 18028, 26815 },
        { 59415, 14855 },
        { 40293, 34797 },
        {  6292,  5951 },
        { 31961, 33808 },
        { 65425, 30667 },
        { 58621, 33373 },
        { 30798, 34558 },
        { 56595, 27003 },
        { 21522, 36268 },
        { 31574, 38417 },
        { 55532, 50988 },
        { 15726, 23479 },
        { 64982, 52806 },
        { 64386, 19906 },
        { 50915, 57484 },
        { 38322, 56885 },
        { 21853, 25436 },
        { 37894, 11756 },
        {  7865, 18734 },
        { 46023, 39441 },
        { 22621, 43830 },
        { 36517, 19001 },
        { 20782, 63993 },
        { 29250, 63913 },
        {  8530, 56189 },
        { 31507, 16655 },
        { 63002, 48210 },
        { 20605, 43609 },
        { 43367, 49261 },
        { 53367, 18068 },
        { 40975, 42072 },
        { 24653, 60178 },
        { 52335, 16671 },
        { 49599, 51110 },
        { 25537, 57521 },
        { 36030, 22357 },
        { 32634, 40781 },
        { 43255, 18663 },
        { 52962,  4135 },
        { 25203, 21770 },
        { 21866, 45501 },
        { 54665, 42731 },
        { 27452, 38361 },
        { 41137, 43933 },
        { 26533, 37842 },
        { 58981, 55678 },
        { 58627, 26173 },
        { 63625, 64590 },
        { 24065, 14082 },
        { 44021,  9266 },
        { 28206, 33451 },
        { 50533, 57181 },
        { 55629, 26391 },
        { 54949, 36509 },
        {   195, 37407 },
        { 41688, 20532 },
        { 42237, 23674 },
        {  8663, 25519 },
        { 29444, 42230 },
        { 27799, 45221 },
        { 34850, 42223 },
        { 65379, 23075 },
        { 25892,  5386 },
        { 62200, 41910 },
        { 56868, 54033 },
        { 54899, 46879 },
        { 34262, 41102 },
        { 62814, 59754 },
        { 45964, 61750 },
        { 20902,  1788 },
        { 28717, 25945 },
        { 21825, 35962 },
        { 54899, 35074 },
        { 53805, 43558 },
        { 59949, 40748 },
        { 59292, 40913 },
        {  8651, 11678 },
        { 50522, 65395 },
        { 52685, 55207 },
        { 45597, 54755 },
        { 62570, 31143 },
        { 36030, 19158 },
        { 58093, 42546 },
        {  5718, 63090 },
        { 32799,  6434 },
        { 19110, 10348 },
        { 56355, 36478 },
        { 28120, 29065 },
        { 19110, 50071 },
        { 60426, 63280 },
        { 38851, 37573 },
        { 26986, 48091 },
        { 35557, 26175 },
        { 35082, 46429 },
        { 37908, 14082 },
        {  1769, 16637 },
        { 29195, 52823 },
        { 31965, 21299 },
        { 36505, 60278 },
        { 60408, 39669 },
        { 44955,  4820 },
        { 25009, 43912 },
        { 27239, 24188 },
        { 58118, 36283 },
        { 58615, 32771 },
        { 41807, 18390 },
        { 31990, 43802 },
        {  4693, 56777 },
        { 48251, 60927 },
        { 65533, 58994 },
        { 33786, 35838 },
        { 17521, 39621 },
        {  8555, 33443 },
        { 38709, 17374 },
        { 29545, 41087 },
        { 26705, 19037 },
        { 59851, 39871 },
        { 39881, 41205 },
        { 29523,  9668 },
        {  6877, 59538 },
        { 43415, 34823 },
        { 62856, 31310 },
        { 21374, 33687 }
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
    static_assert(sizeof(table[0])/sizeof(table[0][0]) == 2u, "");
    static_assert(sizeof(table)/sizeof(table[0][0]) == SIZE*2, "");
};
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<62, 4, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<62, 4, DUMMY>::table[][2];


}}  // end namespace

#endif
