// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES64_4_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES64_4_H_INCLUDED


#include "hurchalla/factoring/detail/miller_rabin_bases/MillerRabinBases.h"
#include "hurchalla/util/compiler_macros.h"
#include <cstdint>
#include <cstddef>
#include <type_traits>
#include <array>

namespace hurchalla { namespace detail {


// This 2.75 KB hash table lets you determine the primality of any 64 bit
// unsigned int number, via miller-rabin primality testing using 4 bases.
// This test permits numbers that are even, as well as of course odd numbers.
// Note however that montgomery arithmetic (which is one way to implement the
// miller-rabin test) always requires odds.

// see MillerRabinBases.h for why this template uses a DUMMY parameter.
template <typename DUMMY>
struct MillerRabinBases<64, 4, DUMMY> {
    static_assert(std::is_same<DUMMY, void>::value, "");
public:
    // 'num' is the unsigned 64 bit number being tested for primality.
    // This function returns 4 bases that can be used by miller-rabin testing
    // to correctly determine (non-probabilistically) the primality of num.
    static HURCHALLA_FORCE_INLINE
    std::array<std::uint16_t, 4> get(std::uint64_t num)
    {
        // I generated/verified the hash table and bases.  See README.TXT
        using std::uint32_t;
        std::array<std::uint16_t, 4> bases;
        bases[0] = 2;
        bases[1] = 15;
        uint32_t mask = (static_cast<uint32_t>(1) << 24) - 1u;
        uint32_t hash_bucket = ((static_cast<uint32_t>(num) & mask) * 11) >> 18;
        bases[2] = table[hash_bucket][0];
        bases[3] = table[hash_bucket][1];
        return bases;
    }
private:
    static constexpr std::size_t SIZE = 704;
#if defined(__INTEL_COMPILER) || defined(_MSC_VER)
    // unless array size is explicit, icc and msvc fail on sizeof() later
    static constexpr std::uint16_t table[SIZE][2] = {
#else
    static constexpr std::uint16_t table[][2] = {
#endif
        {    83, 42707 },
        {  8293, 22045 },
        {  1194, 16331 },
        {   246, 10073 },
        {   161, 29013 },
        {    17, 24719 },
        {   963, 58923 },
        {    13, 12877 },
        {  1187, 46205 },
        {  7091, 49049 },
        {   239,  1989 },
        {  7998, 64035 },
        {   907, 10209 },
        {   946,  8607 },
        {  2435, 34035 },
        {   143,  9395 },
        {  2106, 65483 },
        {  3214, 41366 },
        {   106, 40214 },
        {  1115, 36107 },
        {   823,  5669 },
        {   637,  1294 },
        {  2034, 32793 },
        {   975, 14171 },
        {  1799,  6341 },
        {  1986, 11658 },
        {  2851,  7246 },
        {   389, 50907 },
        {   610, 53019 },
        {  1109, 24506 },
        {  1625, 44982 },
        {  9129, 54955 },
        {  3306, 27333 },
        {  2599, 38053 },
        {  1239,  2314 },
        {  2617, 52999 },
        {  3797, 55246 },
        {  2481,  8799 },
        {  7462, 23622 },
        {  1053, 26006 },
        {    71, 55807 },
        {    13, 25369 },
        {   630, 19427 },
        { 22359, 22514 },
        {    78, 14702 },
        {    22, 25537 },
        {  2823, 29123 },
        {  4214, 32971 },
        {  3529, 30233 },
        {   249, 47815 },
        {    78, 21473 },
        { 10485, 42603 },
        {  5977, 43819 },
        {   502, 26727 },
        { 12123, 48562 },
        {   861, 58666 },
        {  1366, 26403 },
        { 24315, 42309 },
        {   353, 48403 },
        {   214, 12207 },
        {   118, 51925 },
        {   435,  2763 },
        {  1284, 13289 },
        {  3134, 41427 },
        {  1294, 15711 },
        {   293, 16893 },
        {  2842, 25503 },
        {  3547, 25966 },
        {   458, 33786 },
        {   921, 36379 },
        {   542, 33475 },
        {  1939, 40149 },
        {  1711, 16926 },
        {   283, 59854 },
        {   403, 24061 },
        {  2971, 12114 },
        {   426,  2504 },
        {  2402, 26737 },
        {  1353,  2342 },
        {  2545, 29167 },
        {   445, 35566 },
        {  2882, 60539 },
        {   770, 65087 },
        {   909, 28419 },
        {   883, 59838 },
        {   130, 15239 },
        {  1071, 37266 },
        {   101, 57778 },
        {  4991, 34635 },
        {  1951, 19614 },
        {  3206, 53863 },
        {   835, 43146 },
        {  3253, 33751 },
        {   309, 45309 },
        {   267, 46683 },
        {   557, 26701 },
        {  1094,  7171 },
        {  6378, 48563 },
        {  4658, 53086 },
        {  1233, 16393 },
        {   223, 14350 },
        {  1418, 50759 },
        {  1011,  3322 },
        {  1113, 64159 },
        {  1966, 54862 },
        {    56, 46381 },
        {   534, 63338 },
        {   259, 50461 },
        {  4006, 35134 },
        {   988,  4246 },
        {    33, 24306 },
        {  2249, 31955 },
        {  1150, 37779 },
        {   521, 52374 },
        {  6307, 49386 },
        {   261,  4523 },
        {  2843, 48423 },
        {  3295, 42087 },
        {  2421, 64983 },
        {   247,  9618 },
        {  7867, 63350 },
        {   403,  4906 },
        {  4075, 20067 },
        {   327, 53621 },
        {  2482, 18511 },
        {   319, 26626 },
        {   235, 37250 },
        {   851, 57775 },
        {  5277,  6625 },
        {  3306, 36465 },
        {   582, 40369 },
        {  1246, 37939 },
        {  2247, 48717 },
        {  1967, 33710 },
        {  1685, 52125 },
        {  1987, 33711 },
        {  3830, 19105 },
        {   695,  7486 },
        {   899, 11077 },
        {  4875, 65327 },
        {   271, 51875 },
        {    11, 11038 },
        {  2713, 49266 },
        {  2145, 51973 },
        {  2478, 63035 },
        {  2106, 17654 },
        {  6099, 14195 },
        {   754, 28474 },
        {   333, 40353 },
        {  1237, 32610 },
        {   869, 14425 },
        { 15175, 39646 },
        {  1405, 28158 },
        {  5273, 14029 },
        {  7386, 28371 },
        {   582,  1069 },
        {  2021, 11589 },
        {   319, 15006 },
        {  2139,  9194 },
        {   913, 19233 },
        {  4589, 27587 },
        {   341, 11573 },
        {   627, 11553 },
        {   989, 16421 },
        {  1137, 63597 },
        {  5434, 16275 },
        {  1326,  2765 },
        {   827,  2129 },
        {   137, 34190 },
        {  1211, 37233 },
        {   313, 54830 },
        {   502, 25729 },
        {   827, 48585 },
        {  2371, 23114 },
        {   223, 29823 },
        {   325, 43043 },
        {  2026, 36229 },
        {   325, 25391 },
        {  1421, 11350 },
        {   774, 54073 },
        {    34, 30073 },
        {   701, 49531 },
        {   511, 42210 },
        {   742,  9965 },
        {  4174, 25115 },
        {  8045, 63353 },
        {  5421, 47395 },
        {   335, 53650 },
        {  7026, 42897 },
        {  1094, 29363 },
        {   574, 25879 },
        {   398, 51911 },
        {  1867,  2947 },
        {    39,   318 },
        {  3725, 31870 },
        {  8333, 51397 },
        {  1207, 31318 },
        {    37, 21429 },
        {  1730, 28010 },
        {   351, 27219 },
        {  1154,  5838 },
        {   371, 40131 },
        {   291,  9313 },
        {   719, 58310 },
        {  7975, 17918 },
        {   223, 33130 },
        {   959, 52521 },
        {   979, 31930 },
        {   534, 55270 },
        {   970, 35926 },
        {   754, 52153 },
        {   611, 32338 },
        {  5981,  7869 },
        {   886, 53885 },
        {  1221, 49046 },
        {  2017, 59938 },
        { 17153, 25326 },
        {  2142,  4709 },
        {    83, 35895 },
        {   218, 19587 },
        {   367, 50117 },
        {   658, 20151 },
        {  3370, 23650 },
        {  3201, 55598 },
        {   337, 42269 },
        {   463, 48297 },
        {   266, 41646 },
        {  2173, 59310 },
        {   498,  6971 },
        {  1631, 58604 },
        {  4629, 59979 },
        {   278, 56570 },
        {   630, 15522 },
        {   905,  3942 },
        {   181, 59739 },
        {  4379, 49151 },
        {  1613, 24474 },
        {  2562,  7367 },
        {   439, 45314 },
        {   238,  3531 },
        {  2455, 24965 },
        {  1450, 20017 },
        {  2145,  4981 },
        {   951, 43402 },
        {   215,  7254 },
        {  1190, 18199 },
        {  1111, 50673 },
        {  1322, 48911 },
        {   198, 12067 },
        {   446, 11677 },
        {  1218, 22629 },
        {   639, 14394 },
        {  3557, 13416 },
        {    62, 28403 },
        {   238, 62138 },
        {   557, 52053 },
        {  1505, 60851 },
        {   154, 38778 },
        {  3155, 46342 },
        {  3754, 37262 },
        {   475, 11614 },
        {  2959, 15174 },
        {  1181, 64406 },
        {    26, 25555 },
        {  4545, 29049 },
        {  1921, 60574 },
        {   459, 49958 },
        {   589, 24723 },
        {  1877, 20538 },
        {   763, 55179 },
        {    78, 30083 },
        {   583, 64013 },
        {   711, 37573 },
        {   769, 23461 },
        {   471, 44133 },
        {    77, 17021 },
        {  3230, 14717 },
        {  8951, 33909 },
        {  2134, 19101 },
        {  3605,  4295 },
        {    74, 38983 },
        {  2065, 29571 },
        {   743, 42823 },
        {  8006, 40026 },
        {   138, 19383 },
        {  1526, 46623 },
        {    41, 38883 },
        {  2324, 48797 },
        {   309, 33471 },
        {  5261, 53286 },
        {   183,  2301 },
        {  1062,  1857 },
        {   142, 20717 },
        {  1790, 11198 },
        {   785, 17974 },
        {   805, 28619 },
        {    91, 21475 },
        {  1975, 27746 },
        {   554, 34078 },
        {  1430, 54942 },
        {  3337, 11793 },
        {    87, 39209 },
        {  2233, 63814 },
        {  4294,  9555 },
        {   741, 64354 },
        {   211, 39418 },
        {  1145, 27241 },
        {  1207, 44246 },
        {   259, 36073 },
        {   237, 36971 },
        {   926, 42118 },
        {   946, 26586 },
        {  8359, 39115 },
        {  4479, 54426 },
        {  1327, 39325 },
        {   583, 56887 },
        { 10223, 13397 },
        {  1515, 36905 },
        {  3810, 22795 },
        {   934, 54130 },
        {  1135, 34921 },
        {  1081,  4087 },
        {   235,  3993 },
        {  2221,  5613 },
        {   685, 21243 },
        {  3222, 48618 },
        {    13, 33459 },
        {   109, 65381 },
        {  2639, 50631 },
        {   457, 59462 },
        {   179, 34902 },
        {  1354, 51562 },
        {  3759, 65282 },
        {   444, 50274 },
        {    26, 26066 },
        {   691, 48911 },
        {   554, 33817 },
        {  2850, 61134 },
        {  1357, 20067 },
        {  2909, 28690 },
        {   130, 40905 },
        {  1091, 49131 },
        {    57, 56943 },
        {   678, 41295 },
        {   658,  6178 },
        {   569, 22126 },
        {  2713, 36945 },
        {  2902,  6905 },
        {   485, 24627 },
        {    51,  4606 },
        {   326, 45149 },
        {  5683, 41055 },
        {  4873, 21629 },
        {  1429, 15587 },
        {   964, 24998 },
        {  1435, 29810 },
        {  1867, 29958 },
        {  5405, 26865 },
        {   631, 30485 },
        {  1053, 35015 },
        {    39, 19262 },
        {   163, 47745 },
        {   191, 48474 },
        {   733, 44439 },
        {  1046, 42113 },
        {  1603, 48833 },
        {    47, 47275 },
        {  1361, 37485 },
        {    34, 38069 },
        {  1913, 56981 },
        {   831, 53901 },
        {   442,  3730 },
        {  3838, 18929 },
        {  1382, 57046 },
        {  5943, 14578 },
        {  4273, 21585 },
        {  2278, 42533 },
        {  4619, 37723 },
        {  1195, 33238 },
        {  1305, 61447 },
        {  1762, 27583 },
        {  3113, 50767 },
        {  9967, 51951 },
        {  2175, 22406 },
        {   655, 17427 },
        {  2117, 29273 },
        {  5761, 31727 },
        {  2901,  6218 },
        {  1106,  7541 },
        {  1581, 16002 },
        {  1614, 28207 },
        {  2551, 11647 },
        {   251, 25541 },
        {  1303, 49529 },
        {  4402, 47794 },
        {  2959, 50959 },
        {   541, 12418 },
        {    38, 47123 },
        {   894,  7066 },
        { 12386, 65315 },
        {  1270, 22021 },
        {  2183, 25321 },
        {  2854,  3430 },
        {   919, 56613 },
        {  1474, 31483 },
        {   574, 10801 },
        {  1301, 33410 },
        {   195,  6775 },
        {  1507, 14828 },
        {  5975, 28471 },
        {  7329, 22510 },
        {    74, 30518 },
        {   351, 44923 },
        {  3283, 19027 },
        {  3082, 14295 },
        {  2659, 44733 },
        {  6461, 58511 },
        {   501, 11922 },
        {   603,  1238 },
        {   646,  4505 },
        {  2339, 10150 },
        {   318, 14447 },
        {  6098, 20675 },
        {   487, 23423 },
        {  3498, 45574 },
        {  4247, 16297 },
        {  2948, 61339 },
        {  1657, 14125 },
        {  6547, 37341 },
        {  4122, 20503 },
        {   963, 22746 },
        {  1211,  8982 },
        {   813, 59331 },
        {   123, 49465 },
        {  3705, 41155 },
        {   263, 16894 },
        {  6174, 30967 },
        {  1085, 11275 },
        {  6155, 47549 },
        {  2479, 51721 },
        {   198,  5675 },
        {  5861, 26267 },
        {   231, 30782 },
        {  1021,  8930 },
        {   230, 12166 },
        {     7, 58246 },
        {  3018, 45833 },
        {   365, 53877 },
        {  1971, 13357 },
        {   117,  8089 },
        {   239, 39374 },
        {  1157, 34741 },
        {   487, 56362 },
        {  2171, 26027 },
        {  1061,  7898 },
        {  1050, 27069 },
        {  1685, 38875 },
        {  1498, 27721 },
        { 15557, 16647 },
        {   197, 23182 },
        {  2617, 55889 },
        {  1053, 13187 },
        {  1343, 25499 },
        {   118, 47206 },
        {   885, 33046 },
        {   685, 59477 },
        {  3163, 25963 },
        {  4387, 35807 },
        {  1018, 35034 },
        {    65, 14241 },
        {  6150, 64874 },
        {  1582, 27558 },
        {  3162, 40911 },
        {  1466, 56973 },
        {   602,  9699 },
        {  9138, 40879 },
        {  2935, 23721 },
        { 11403, 14049 },
        {  4183, 34135 },
        {   130, 17453 },
        {   333, 53431 },
        {  1557, 10486 },
        {    77, 56649 },
        {   183, 35715 },
        {   763, 54622 },
        {   355, 59950 },
        {    82, 50695 },
        {  1333, 57975 },
        {  3295, 15081 },
        {  1422, 37841 },
        {   786, 41754 },
        {   197, 58659 },
        {  1185,  6123 },
        { 11743, 35837 },
        {  2314, 36157 },
        {   766,  5159 },
        {   374, 16270 },
        {  1921, 54565 },
        { 10215, 41357 },
        {  1571, 43921 },
        {   497, 59315 },
        {  2158, 18694 },
        {  1281, 56115 },
        {   986, 51622 },
        {   517,  8099 },
        {  1549, 48309 },
        {  1237, 13385 },
        {  1383,  2797 },
        {   830, 27019 },
        {   546,  6399 },
        {   847, 54051 },
        {  3318, 44995 },
        {  1418, 23709 },
        {  1018, 35377 },
        {  3419, 10941 },
        {  3153, 30197 },
        {   930, 55222 },
        {  4361, 49901 },
        {   157,  5231 },
        {  3983, 37165 },
        {   157, 61719 },
        {  2945,  3326 },
        {  4211, 36629 },
        {  2225, 30495 },
        {   694, 33797 },
        {   809, 56490 },
        {  3905, 52233 },
        {   145, 42846 },
        {   378,  6549 },
        {  1827,  8351 },
        {  4749, 62655 },
        {  1146,  3387 },
        {   386, 50338 },
        {  8406, 16404 },
        {    94, 33362 },
        {  7323, 19639 },
        { 10741, 36482 },
        {  7834, 24170 },
        {   495, 65431 },
        {  4121, 48753 },
        {  2185, 50090 },
        {  3809,  8207 },
        {  1147, 16783 },
        {   305, 39865 },
        {  1647, 44917 },
        {  1882, 51837 },
        {    82,  7097 },
        {   663, 55877 },
        {  5737, 23417 },
        {  3166, 52913 },
        {  1101, 57298 },
        {  1201, 56001 },
        {   265, 16353 },
        {  1221, 44566 },
        {  2009, 21970 },
        {   262,  3759 },
        { 14866, 37881 },
        {  1148,  4571 },
        {  2426, 25074 },
        {   482, 15999 },
        {  4575, 25315 },
        {   681, 33798 },
        {  1294, 27915 },
        {  1595, 18763 },
        {   204, 36611 },
        {   269, 57551 },
        {   219, 32074 },
        {   227,  6145 },
        {   409, 30546 },
        {   842, 58875 },
        {  8361, 14399 },
        {  1477, 18333 },
        {   846, 64093 },
        {  3505, 21589 },
        {  6877, 46886 },
        {   466, 13627 },
        {  1294, 36047 },
        {   901, 16649 },
        {   103, 30005 },
        {   751, 19313 },
        {  1171, 16317 },
        {  1311, 11521 },
        {  1427,  8871 },
        {    58, 29031 },
        { 20223, 63959 },
        {  2787, 12826 },
        {  6835, 62086 },
        {   483, 32238 },
        {   919, 43683 },
        {  1267, 44146 },
        {   583,  5690 },
        {  1550, 51701 },
        {  2267, 24227 },
        {  1795, 53187 },
        {  2191, 46006 },
        {  1201, 45163 },
        {  3821, 27009 },
        {  2019, 10217 },
        {  7069, 20587 },
        {  1635, 58961 },
        {  3682, 64967 },
        {  2478, 62902 },
        {   154, 46986 },
        {  4963, 56799 },
        {  2202, 43437 },
        {  3578, 57385 },
        {   746, 11487 },
        {  2277, 15565 },
        {   227, 60925 },
        {  3153, 59849 },
        {  1669, 51321 },
        {   174, 12666 },
        { 15894, 31194 },
        {  1913, 56647 },
        {   910, 61541 },
        {   871, 61891 },
        {  1576, 42309 },
        {  1495, 45134 },
        {  1509, 54575 },
        {   563, 59417 },
        {   305, 28574 },
        {  1141, 62510 },
        {   914, 54630 },
        {  5671, 45411 },
        { 12071, 17921 },
        {   227, 18607 },
        {   183, 16233 },
        {  2570, 15079 },
        {  2959, 64007 },
        {   307, 19925 },
        { 11457, 16254 },
        {   218, 25327 },
        {   325, 55321 },
        {  2035, 22754 },
        {   303, 10389 },
        {   408, 17733 },
        {   287, 20389 },
        {    39, 55106 },
        {  2006,  5375 },
        {  2981, 49017 },
        {  3250, 58277 },
        {   869, 35973 },
        {  1625, 34138 },
        {  5914, 64455 },
        {  1842,  4171 },
        {  6837, 44723 },
        {  1955, 18479 },
        {  3705, 65278 },
        {  1434, 29198 },
        {  1526, 22770 },
        {  3679, 65317 },
        {  1142, 28531 },
        {   899,  6098 },
        {  7566, 32977 },
        {   791, 13471 },
        {  1537, 37743 },
        {  1322, 55913 },
        {   146, 19870 },
        {   346, 58182 },
        {   962, 47078 },
        {   677,  5027 },
        {  1446, 13811 },
        {  1615, 49954 },
        {   333, 20022 },
        {  2718, 11563 },
        {  3329, 48641 },
        {  4499, 39250 },
        {  2347, 39169 },
        {  3715, 26758 },
        {  1022, 36218 },
        {   175, 15501 },
        {   345, 48562 },
        {  2291, 52951 },
        {  1863, 13383 },
        {  1231, 34390 },
        {   314, 21935 },
        {   978, 20469 },
        {  1159, 38503 },
        {  3843, 24895 },
        {    13,  2915 },
        {   117, 61939 },
        {   767, 18867 },
        {  2489, 51409 },
        {   463, 61261 },
        {  1626, 49230 },
        {    77, 25793 },
        {  1534, 42238 },
        {   685, 21923 },
        {  1491, 20335 },
        {   339, 15531 },
        {  4295, 18737 },
        {  1747, 56253 },
        {  1166, 56449 },
        {  5798, 54265 },
        {  2043,  3646 },
        {  6753, 53697 },
        {   837, 53725 },
        {  9769, 32183 },
        {  1409, 12317 },
        {  1854, 54687 },
        {   261, 36154 },
        { 47981, 55781 },
        {  3739, 21526 },
        {   561, 37373 }
    };
    static_assert(sizeof(table)/sizeof(table[0]) == SIZE, "");
    static_assert(sizeof(table[0])/sizeof(table[0][0]) == 2u, "");
    static_assert(sizeof(table)/sizeof(table[0][0]) == SIZE*2, "");
};
template <typename DUMMY>
constexpr std::size_t MillerRabinBases<64, 4, DUMMY>::SIZE;
template <typename DUMMY>
constexpr std::uint16_t MillerRabinBases<64, 4, DUMMY>::table[][2];


}}  // end namespace

#endif
