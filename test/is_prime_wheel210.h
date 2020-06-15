
#ifndef HURCHALLA_FACTORING_IS_PRIME_WHEEL210_H_INCLUDED
#define HURCHALLA_FACTORING_IS_PRIME_WHEEL210_H_INCLUDED


#include "integer_sqrt.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>

namespace hurchalla { namespace factoring {


// This is a slight optimization compared to is_prime_bruteforce(), but it's
// just a constant time improvement and remains mostly a brute force approach.
// The algorithm is an adaptation of Wheel factorization, where we return false
// at the first factor found, or return true if no factors are found less than
// or equal to the sqrt of the input parameter.  For more info see:
// https://en.wikipedia.org/wiki/Wheel_factorization
template <typename T>
bool is_prime_wheel210(T x)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    HPBC_PRECONDITION2(x >= 0);

    if (x < 2) return false;
    // We test divisors to 13, to cover all possible factors for types uint8_t
    // and int8_t, meaning these types never need to use the wheel (int8_t would
    // overflow if it used the wheel).  Note: for larger types we could in
    // principle have stopped testing at 7 here and started the wheel at 11, but
    // even ignoring the incompatibility with 8 bit types, we are probably
    // getting slightly better performance on average by testing up to 13.
    if (x%2 == 0) return (x==2);
    if (x%3 == 0) return (x==3);
    if (x%5 == 0) return (x==5);
    if (x%7 == 0) return (x==7);
    if (x%11 == 0) return (x==11);
    if (x%13 == 0) return (x==13);

    T s = integer_sqrt(x);

    // The wheel spans the 210 number range [17, 227), skipping all multiples
    // of 2,3,5,7
    static constexpr uint8_t wheel[] = { 17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
        59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 121, 127,
        131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187,
        191, 193, 197, 199, 209, 211, 221, 223 };
    static constexpr size_t wheel_len = sizeof(wheel)/sizeof(wheel[0]);
    // Use the wheel to skip division by any multiple of 2,3,5,7 in the loop
    // below - we already tested 2,3,5,7 and found they were not factors.
    // Note the wheel pattern cycles every 2*3*5*7 == 210 numbers.
    static const uint8_t cycle_len = 210;
    for (T start=0; start + wheel[0] <= s;
                                  start = static_cast<T>(start + cycle_len)) {
        for (size_t i=0; i < wheel_len; ++i) {
            T maybe_factor = static_cast<T>(start + wheel[i]);
            if (x % maybe_factor == 0)
                return (x == maybe_factor);
        }
    }
    // we went past sqrt(x) without finding a factor, so x must be prime
    return true;
}


}}  // end namespace

#endif
