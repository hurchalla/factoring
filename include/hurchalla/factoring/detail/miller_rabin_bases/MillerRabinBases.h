// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---
// Author: Jeffrey Hurchalla

#ifndef HURCHALLA_FACTORING_MILLER_RABIN_BASES_H_INCLUDED
#define HURCHALLA_FACTORING_MILLER_RABIN_BASES_H_INCLUDED


#include <cstddef>

namespace hurchalla { namespace detail {


// Implementation note: this primary template uses a DUMMY parameter (defaulted
// to void), so that its specializations can be partial specializations
// templated on DUMMY rather than needing to be full explicit specializations.
// This allows the specializations to out-of-line define their static member
// variables within their headers without any danger of ODR violations.
// Unfortunately out-of-line definitions are necessary up until C++17 (17 allows
// the inline keyword to be used on member variables, which is much cleaner).
// For more information, see
// https://stackoverflow.com/questions/8016780/undefined-reference-to-static-constexpr-char
// and also
// https://stackoverflow.com/questions/45210631/odr-violation-with-template-specializations
// https://stackoverflow.com/questions/34552380/why-cs-vector-templated-class-doesnt-break-one-definition-rule
// https://en.cppreference.com/w/cpp/language/static
// https://en.cppreference.com/w/cpp/language/storage_duration
// https://en.cppreference.com/w/cpp/language/definition
// https://stackoverflow.com/questions/8016780/undefined-reference-to-static-constexpr-char
//
// Additionally, having partial specializations (via DUMMY) ensures that if one
// of the partially specialized structs is never used, then that struct's static
// variables (e.g. a large constexpr array table) will not exist and will thus
// use no memory.  A partial specialization could also require DUMMY to always
// be void via static_assert, which would ensure that only one instantiation is
// possible, which would in turn ensure that there could never be more than one
// copy of its static variables in memory.

template <int LOG2_MODULUS_LIMIT, std::size_t TOTAL_BASES, typename DUMMY=void>
struct MillerRabinBases {};


}}  // end namespace

#endif
