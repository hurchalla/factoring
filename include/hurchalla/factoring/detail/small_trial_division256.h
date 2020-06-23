// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_SMALL_TRIAL_DIVISION256_H_INCLUDED
#define HURCHALLA_FACTORING_SMALL_TRIAL_DIVISION256_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#include <type_traits>
#include <limits>

namespace hurchalla { namespace factoring {


// small_trial_division256 guarantees it will try at least all possible prime
// factors less than 256.  It might someday try factors greater than 256 also.
//
// small_trial_division256 return values-
// 0: The function couldn't find any factors and couldn't tell if x is composite
//     --fca is left empty.
// 1: Either the function determined x is not composite (and placed x as the
//     single element in fca), or it found all factors (and placed them all in
//     fca).  Either way, the return value is the quotient of x divided by all
//     elements in fca.
//   *Note this indicates all possible factoring is complete.
//     --fca.size() == 1 indicates x is non-composite.  The single element is x.
//     --fca.size() > 1 indicates all factors were found and placed in fca.
// >1: This return value is the quotient of x divided by all the factors in fca.
//     This return value may be factorable, or it may be prime (we don't know).
//     --fca has all the factors found so far.
//
// Postcondition of small_trial_division256(): if the argument x passed in is
// <= 65535, then the returned quotient == 1 and the final fca.size() >= 1.



// The more specific overloads of small_trial_division256() below (if available
// for your type) should provide better performance than this generic version.
template <class C, typename T>
typename std::enable_if<std::is_same<typename C::value_type, T>::value, T>::type
small_trial_division256(C& fca, const T x)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!ma::ma_numeric_limits<T>::is_signed, "");

    // We'll populate small_primes with all primes less than 256.
    const uint8_t small_primes[] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
        41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
        113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
        193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251 };
    constexpr size_t array_len = sizeof(small_primes)/sizeof(small_primes[0]);

    if (x < 4) {
        fca.push(x);
        return 1;  // x is not composite
    }

    T q = x;
    HPBC_INVARIANT2(q > 0);
    for (size_t i=0; i<array_len; ++i) {
        const unsigned int prime = small_primes[i];
        // If no primes <= sqrt(q) are factors, q is prime or 0 or 1
        if (prime * prime > q) {
            if (q > 1)
                fca.push(q);
            return 1;
        }
        while (q % prime == 0) {
            fca.push(static_cast<T>(prime));
            q = static_cast<T>(q / prime);
        }
        HPBC_INVARIANT2(q > 0);
    }

    if (q <= 65535u) {  // 65535 is 256*256-1 (we checked all factors < 256)
        if (q > 1)  // q is not composite since we checked all factors<256
            fca.push(q);
        return 1;
    }
    if (fca.size() == 0)
        return 0;  // we found no factors and we don't know if x is composite
    else {
        HPBC_ASSERT2(q > 1);
        return q;  // q is the quotient of x divided by all the factors found.
    }
}




template <class C>
typename std::enable_if<std::is_same<typename C::value_type, uint8_t>::value,
                        uint8_t>::type
small_trial_division256(C& fca, const uint8_t x)
{
    if (x < 4) {
        fca.push(x);
        return 1;
    }
    unsigned int q = x;  // use unsigned int to avoid integral promotion hassles

    // we only need to try primes < 16, since 16*16==256 covers all of uint8_t.
    while ((q & 1) == 0) {  // equivalent to (q % 2 == 0)
        fca.push(2);
        q = q >> 1;
    }
    while (q % 3 == 0) {
        fca.push(3);
        q /= 3;
    }
    while (q % 5 == 0) {
        fca.push(5);
        q /= 5;
    }
    while (q % 7 == 0) {
        fca.push(7);
        q /= 7;
    }
    while (q % 11 == 0) {
        fca.push(11);
        q /= 11;
    }
    while (q % 13 == 0) {
        fca.push(13);
        q /= 13;
    }
    HPBC_ASSERT2(q > 0);
    if (q > 1)
        fca.push(static_cast<uint8_t>(q));
    return 1;
}



template <typename T>
struct PrimeInfoPair {
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!std::numeric_limits<T>::is_signed, "");
    T inv_mod_R;   // R is shorthand for (2 ^ std::numeric_limits<T>::digits)
    T maxuint_div_prime;     // maxuint means std::numeric_limits<T>::max()
};


// Helper template function for the small_trial_division256() functions below.
template <class C>
typename C::value_type
perform_trial_divisions(C& fca, const typename C::value_type x, const uint8_t*
                      small_primes, const PrimeInfoPair<typename C::value_type>*
                      primes_info, size_t array_len)
{
    using T = typename C::value_type;
    namespace ma = montgomery_arithmetic;
    using P = typename ma::safely_promote_unsigned<T>::type;

    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!std::numeric_limits<T>::is_signed, "");
    HPBC_PRECONDITION2(small_primes != nullptr);
    HPBC_PRECONDITION2(primes_info != nullptr);
    HPBC_PRECONDITION2(array_len > 0);
    HPBC_PRECONDITION2(small_primes[0] == 3);  // must start primes at 3, not 2.
    // we must at least check all primes less than 256
    HPBC_PRECONDITION2(small_primes[array_len-1] >= 251);

    if (x < 4) {
        fca.push(x);
        return 1;  // x is not composite
    }

    T q = x;
    while ((q & 1) == 0) {  // equivalent to (q % 2 == 0)
        fca.push(2);
        q = static_cast<T>(q >> 1);
    }
    // See Hacker's Delight 2nd edition by Henry Warren, Section 10-17 "Test for
    // Zero Remainder after Division by a Constant", for a description of the
    // algorithm used below and for a proof of its correctness.

    HPBC_INVARIANT2(q > 0);
    for (size_t i=0; i<array_len; ++i) {
        const unsigned int prime = small_primes[i];
        // If no primes <= sqrt(q) are factors, q is prime or 0 or 1
        if (prime * prime > q) {
            if (q > 1)
                fca.push(q);
            return 1;
        }

        T inv = primes_info[i].inv_mod_R;
        T tmp = static_cast<T>(static_cast<P>(q) * inv);
        // Test  tmp <= UINT32_MAX / prime,  which is equivalent to testing
        // (q % prime == 0).  For details see Hacker's Delight.
        while (tmp <= primes_info[i].maxuint_div_prime) {
            fca.push(static_cast<T>(prime));
            // Since the prime divides q, tmp=q*inv equals q/prime.
            q = tmp;
            tmp = static_cast<T>(static_cast<P>(q) * inv);
        }
        HPBC_INVARIANT2(q > 0);
    }

    if (q <= 65535u) {  // 65535 is 256*256-1 (we checked all factors < 256)
        if (q > 1)  // q is not composite since we checked all factors<256
            fca.push(q);
        return 1;
    }
    if (fca.size() == 0)
        return 0;  // we found no factors and we don't know if x is composite
    else {
        HPBC_ASSERT2(q > 1);
        return q;  // q is the quotient of x divided by all the factors found.
    }
}



template <class C>
typename std::enable_if<std::is_same<typename C::value_type, uint16_t>::value,
                        uint16_t>::type
small_trial_division256(C& fca, const uint16_t x)
{
    // We'll populate small_primes with all primes (other than 2) less than 256.
    const uint8_t small_primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
        41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
        113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
        193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251 };
 #if 0
// I generated the uint16_t data for primes_info as follows:
    constexpr size_t array_len = sizeof(small_primes)/sizeof(small_primes[0]);
    for (size_t i=0; i<array_len; ++i) {
        namespace ma = hurchalla::montgomery_arithmetic;
        uint16_t inv = (uint16_t)0 -
                        ma::negative_inverse_mod_r((uint16_t)(small_primes[i]));
        std::cout << "{ UINT16_C(" << inv << "), UINT16_C(" <<
                                      UINT16_MAX / small_primes[i] << ") }, \n";
    }
    return 0;
#endif
    static const PrimeInfoPair<uint16_t> primes_info[] = {
        { UINT16_C(43691), UINT16_C(21845)}, // corresponding small_primes[]==3,
        { UINT16_C(52429), UINT16_C(13107)}, // 5,
        { UINT16_C(28087), UINT16_C(9362) }, // 7,
        { UINT16_C(35747), UINT16_C(5957) }, // 11,
        { UINT16_C(20165), UINT16_C(5041) }, // 13, etc...
        { UINT16_C(61681), UINT16_C(3855) },
        { UINT16_C(51739), UINT16_C(3449) }, 
        { UINT16_C(14247), UINT16_C(2849) }, 
        { UINT16_C(49717), UINT16_C(2259) }, 
        { UINT16_C(31711), UINT16_C(2114) }, 
        { UINT16_C(7085) , UINT16_C(1771) }, 
        { UINT16_C(39961), UINT16_C(1598) }, 
        { UINT16_C(48771), UINT16_C(1524) }, 
        { UINT16_C(18127), UINT16_C(1394) }, 
        { UINT16_C(21021), UINT16_C(1236) }, 
        { UINT16_C(55539), UINT16_C(1110) }, 
        { UINT16_C(38677), UINT16_C(1074) }, 
        { UINT16_C(19563), UINT16_C(978)  }, 
        { UINT16_C(43383), UINT16_C(923)  }, 
        { UINT16_C(61945), UINT16_C(897)  }, 
        { UINT16_C(5807) , UINT16_C(829)  }, 
        { UINT16_C(17371), UINT16_C(789)  }, 
        { UINT16_C(18409), UINT16_C(736)  }, 
        { UINT16_C(41889), UINT16_C(675)  }, 
        { UINT16_C(45421), UINT16_C(648)  }, 
        { UINT16_C(6999) , UINT16_C(636)  }, 
        { UINT16_C(44099), UINT16_C(612)  }, 
        { UINT16_C(2405) , UINT16_C(601)  }, 
        { UINT16_C(49297), UINT16_C(579)  }, 
        { UINT16_C(49023), UINT16_C(516)  }, 
        { UINT16_C(20011), UINT16_C(500)  }, 
        { UINT16_C(30137), UINT16_C(478)  }, 
        { UINT16_C(26403), UINT16_C(471)  }, 
        { UINT16_C(51901), UINT16_C(439)  }, 
        { UINT16_C(32551), UINT16_C(434)  }, 
        { UINT16_C(34229), UINT16_C(417)  }, 
        { UINT16_C(45835), UINT16_C(402)  }, 
        { UINT16_C(42775), UINT16_C(392)  }, 
        { UINT16_C(25381), UINT16_C(378)  }, 
        { UINT16_C(44667), UINT16_C(366)  }, 
        { UINT16_C(60829), UINT16_C(362)  }, 
        { UINT16_C(28479), UINT16_C(343)  }, 
        { UINT16_C(36673), UINT16_C(339)  }, 
        { UINT16_C(32269), UINT16_C(332)  }, 
        { UINT16_C(49399), UINT16_C(329)  }, 
        { UINT16_C(22363), UINT16_C(310)  }, 
        { UINT16_C(47903), UINT16_C(293)  }, 
        { UINT16_C(17611), UINT16_C(288)  }, 
        { UINT16_C(48365), UINT16_C(286)  }, 
        { UINT16_C(55129), UINT16_C(281)  }, 
        { UINT16_C(11791), UINT16_C(274)  }, 
        { UINT16_C(61457), UINT16_C(271)  }, 
        { UINT16_C(2611) , UINT16_C(261)  },
        };
    constexpr size_t arraylen = sizeof(small_primes)/sizeof(small_primes[0]);
    static_assert(arraylen == sizeof(primes_info)/sizeof(primes_info[0]), "");

    return perform_trial_divisions(fca, x, small_primes, primes_info, arraylen);
}


template <class C>
typename std::enable_if<std::is_same<typename C::value_type, uint32_t>::value,
                        uint32_t>::type
small_trial_division256(C& fca, const uint32_t x)
{
    // We'll populate small_primes with all primes (other than 2) less than 256.
    const uint8_t small_primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
        41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
        113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
        193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251 };
 #if 0
// I generated the uint32_t data for primes_info as follows:
    constexpr size_t array_len = sizeof(small_primes)/sizeof(small_primes[0]);
    for (size_t i=0; i<array_len; ++i) {
        namespace ma = hurchalla::montgomery_arithmetic;
        uint32_t inv = (uint32_t)0 -
                        ma::negative_inverse_mod_r((uint32_t)(small_primes[i]));
        std::cout << "{ UINT32_C(" << inv << "), UINT32_C(" <<
                                      UINT32_MAX / small_primes[i] << ") }, \n";
    }
    return 0;
#endif
    static const PrimeInfoPair<uint32_t> primes_info[] = {
        { UINT32_C(2863311531), UINT32_C(1431655765)}, // first small_prime is 3
        { UINT32_C(3435973837), UINT32_C(858993459) },
        { UINT32_C(3067833783), UINT32_C(613566756) },
        { UINT32_C(3123612579), UINT32_C(390451572) }, // 11 here
        { UINT32_C(3303820997), UINT32_C(330382099) },
        { UINT32_C(4042322161), UINT32_C(252645135) },
        { UINT32_C(678152731) , UINT32_C(226050910) }, // 19 here
        { UINT32_C(3921491879), UINT32_C(186737708) },
        { UINT32_C(1332920885), UINT32_C(148102320) },
        { UINT32_C(3186588639), UINT32_C(138547332) }, // 31
        { UINT32_C(2437684141), UINT32_C(116080197) },
        { UINT32_C(3247414297), UINT32_C(104755299) },
        { UINT32_C(799063683) , UINT32_C(99882960)  }, // 43
        { UINT32_C(1736263375), UINT32_C(91382282)  },
        { UINT32_C(2350076445), UINT32_C(81037118)  },
        { UINT32_C(2693454067), UINT32_C(72796055)  }, // 59
        { UINT32_C(3238827797), UINT32_C(70409299)  },
        { UINT32_C(128207979) , UINT32_C(64103989)  },
        { UINT32_C(3811027319), UINT32_C(60492497)  }, // 71
        { UINT32_C(3353604601), UINT32_C(58835168)  },
        { UINT32_C(1631000239), UINT32_C(54366674)  },
        { UINT32_C(724452315) , UINT32_C(51746593)  }, // 83
        { UINT32_C(4198451177), UINT32_C(48258059)  },
        { UINT32_C(1594008481), UINT32_C(44278013)  },
        { UINT32_C(2083697005), UINT32_C(42524428)  }, // 101
        { UINT32_C(3544390487), UINT32_C(41698711)  },
        { UINT32_C(2368252995), UINT32_C(40139881)  },
        { UINT32_C(3664513381), UINT32_C(39403369)  }, // 109
        { UINT32_C(266059921) , UINT32_C(38008560)  },
        { UINT32_C(4024418175), UINT32_C(33818640)  },
        { UINT32_C(3376959019), UINT32_C(32786009)  }, // 131
        { UINT32_C(125400505) , UINT32_C(31350126)  },
        { UINT32_C(1884841763), UINT32_C(30899045)  },
        { UINT32_C(2363673277), UINT32_C(28825283)  }, // 149
        { UINT32_C(3214114599), UINT32_C(28443492)  },
        { UINT32_C(738624949) , UINT32_C(27356479)  },
        { UINT32_C(1159377675), UINT32_C(26349492)  }, // 163
        { UINT32_C(3677726487), UINT32_C(25718367)  },
        { UINT32_C(223437605) , UINT32_C(24826400)  },
        { UINT32_C(3647123067), UINT32_C(23994230)  }, // 179
        { UINT32_C(284749213) , UINT32_C(23729101)  },
        { UINT32_C(4002639679), UINT32_C(22486739)  },
        { UINT32_C(2425655105), UINT32_C(22253716)  }, // 193
        { UINT32_C(1024687629), UINT32_C(21801864)  },
        { UINT32_C(4014391543), UINT32_C(21582750)  },
        { UINT32_C(1852331867), UINT32_C(20355295)  }, // 211
        { UINT32_C(3678649119), UINT32_C(19259943)  },
        { UINT32_C(2611037387), UINT32_C(18920560)  },
        { UINT32_C(1200340205), UINT32_C(18755315)  }, // 229
        { UINT32_C(534566745),  UINT32_C(18433336)  },
        { UINT32_C(1132146191), UINT32_C(17970574)  },
        { UINT32_C(285143057),  UINT32_C(17821441)  }, // 241
        { UINT32_C(2583824947), UINT32_C(17111423)  }
        };
    constexpr size_t arraylen = sizeof(small_primes)/sizeof(small_primes[0]);
    static_assert(arraylen == sizeof(primes_info)/sizeof(primes_info[0]), "");

    return perform_trial_divisions(fca, x, small_primes, primes_info, arraylen);
}


template <class C>
typename std::enable_if<std::is_same<typename C::value_type, uint64_t>::value,
                        uint64_t>::type
small_trial_division256(C& fca, const uint64_t x)
{
    // We'll populate small_primes with all primes (other than 2) less than 256.
    const uint8_t small_primes[] = { 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
        41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109,
        113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191,
        193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251 };
 #if 0
// I generated the uint64_t data for primes_info as follows:
    constexpr size_t array_len = sizeof(small_primes)/sizeof(small_primes[0]);
    for (size_t i=0; i<array_len; ++i) {
        namespace ma = hurchalla::montgomery_arithmetic;
        uint64_t inv = (uint64_t)0 -
                        ma::negative_inverse_mod_r((uint64_t)(small_primes[i]));
        std::cout << "{ UINT64_C(" << inv << "), UINT64_C(" <<
                                      UINT64_MAX / small_primes[i] << ") }, \n";
    }
    return 0;
#endif
    static const PrimeInfoPair<uint64_t> primes_info[] = {
        { UINT64_C(12297829382473034411), UINT64_C(6148914691236517205)},
        { UINT64_C(14757395258967641293), UINT64_C(3689348814741910323)},
        { UINT64_C(7905747460161236407) , UINT64_C(2635249153387078802)},
        { UINT64_C(3353953467947191203) , UINT64_C(1676976733973595601)},
        { UINT64_C(5675921253449092805) , UINT64_C(1418980313362273201)},
        { UINT64_C(17361641481138401521), UINT64_C(1085102592571150095)},
        { UINT64_C(9708812670373448219) , UINT64_C(970881267037344821) },
        { UINT64_C(15238614669586151335), UINT64_C(802032351030850070) },
        { UINT64_C(3816567739388183093) , UINT64_C(636094623231363848) },
        { UINT64_C(17256631552825064415), UINT64_C(595056260442243600) },
        { UINT64_C(1495681951922396077) , UINT64_C(498560650640798692) },
        { UINT64_C(10348173504763894809), UINT64_C(449920587163647600) },
        { UINT64_C(9437869060967677571) , UINT64_C(428994048225803525) },
        { UINT64_C(5887258746928580303) , UINT64_C(392483916461905353) },
        { UINT64_C(2436362424829563421) , UINT64_C(348051774975651917) },
        { UINT64_C(14694863923124558067), UINT64_C(312656679215416129) },
        { UINT64_C(5745707170499696405) , UINT64_C(302405640552615600) },
        { UINT64_C(17345445920055250027), UINT64_C(275324538413575397) },
        { UINT64_C(1818693077689674103) , UINT64_C(259813296812810586) },
        { UINT64_C(9097024474706080249) , UINT64_C(252695124297391118) },
        { UINT64_C(11208148297950107311), UINT64_C(233503089540627235) },
        { UINT64_C(11779246215742243803), UINT64_C(222249928598910260) },
        { UINT64_C(17617676924329347049), UINT64_C(207266787345051141) },
        { UINT64_C(11790702397628785569), UINT64_C(190172619316593315) },
        { UINT64_C(4200743699953660269) , UINT64_C(182641030432767837) },
        { UINT64_C(15760325033848937303), UINT64_C(179094602657374287) },
        { UINT64_C(8619973866219416643) , UINT64_C(172399477324388332) },
        { UINT64_C(12015769075535579493), UINT64_C(169236184162472950) },
        { UINT64_C(10447713457676206225), UINT64_C(163245522776190722) },
        { UINT64_C(9150747060186627967) , UINT64_C(145249953336295682) },
        { UINT64_C(281629680514649643)  , UINT64_C(140814840257324821) },
        { UINT64_C(16292379802327414201), UINT64_C(134647766961383588) },
        { UINT64_C(4246732448623781667) , UINT64_C(132710389019493177) },
        { UINT64_C(16094474695182830269), UINT64_C(123803651501406386) },
        { UINT64_C(8062815290495565607) , UINT64_C(122163868037811600) },
        { UINT64_C(6579730370240349621) , UINT64_C(117495185182863386) },
        { UINT64_C(2263404180823257867) , UINT64_C(113170209041162893) },
        { UINT64_C(10162278172342986519), UINT64_C(110459545351554201) },
        { UINT64_C(9809829218388894501) , UINT64_C(106628578460748853) },
        { UINT64_C(17107036403551874683), UINT64_C(103054436165975148) },
        { UINT64_C(3770881385233444253) , UINT64_C(101915713114417412) },
        { UINT64_C(2124755861893246783) , UINT64_C(96579811904238490)  },
        { UINT64_C(8124213711219232577) , UINT64_C(95578984837873324)  },
        { UINT64_C(14513935692512591373), UINT64_C(93638294790403815)  },
        { UINT64_C(2780916192016515319) , UINT64_C(92697206400550510)  },
        { UINT64_C(13900627050804827995), UINT64_C(87425327363552377)  },
        { UINT64_C(7527595115280579359) , UINT64_C(82720825442643729)  },
        { UINT64_C(1950316554048586955) , UINT64_C(81263189752024456)  },
        { UINT64_C(2094390156840385773) , UINT64_C(80553467570784068)  },
        { UINT64_C(7204522363551799129) , UINT64_C(79170575423646144)  },
        { UINT64_C(7255204782128442895) , UINT64_C(77183029597111094)  },
        { UINT64_C(17298606475760824337), UINT64_C(76542506529915151)  },
        { UINT64_C(2939720171109091891) , UINT64_C(73493004277727297)  }
        };
    constexpr size_t arraylen = sizeof(small_primes)/sizeof(small_primes[0]);
    static_assert(arraylen == sizeof(primes_info)/sizeof(primes_info[0]), "");

    return perform_trial_divisions(fca, x, small_primes, primes_info, arraylen);
}


}}

#endif
