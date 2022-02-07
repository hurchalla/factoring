// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef SAMPLE_COMPOSITE_NUMBERS_H_INCLUDED
#define SAMPLE_COMPOSITE_NUMBERS_H_INCLUDED


#include "hurchalla/util/compiler_macros.h"
#include <cstdint>

static std::uint64_t prime_numbers64[] = {
    UINT64_C(2),
    UINT64_C(53),
    UINT64_C(127),
    UINT64_C(67967),
    UINT64_C(67979),
    UINT64_C(40000000003),
    UINT64_C(40000000031),
    UINT64_C(18446744073709551557),
    UINT64_C(18446744073709551533),
    UINT64_C(18446744073709551521)
};
static std::uint64_t nonprime_numbers64[] = {
    UINT64_C(0),
    UINT64_C(1),
    UINT64_C(49),
    UINT64_C(54),
    UINT64_C(55),
    UINT64_C(141),
    UINT64_C(140),
    UINT64_C(256),
    UINT64_C(67968),
    UINT64_C(67969),
    UINT64_C(67981),
    UINT64_C(67982),
    UINT64_C(40000000001),
    UINT64_C(40000000002),
    UINT64_C(40000000005),
    UINT64_C(40000000007),
    UINT64_C(40000000029),
    UINT64_C(40000000027),
    UINT64_C(40000000025),
    UINT64_C(40000000024),
    UINT64_C(8589934592),   // 2^33
    UINT64_C(18446744073709551558),
    UINT64_C(18446744073709551555),
    UINT64_C(18446744073709551554),
    UINT64_C(18446744073709551553),
    UINT64_C(18446744073709551551),
    UINT64_C(18446744073709551549),
    UINT64_C(18446744073709551523),
    UINT64_C(18446744073709551525),
    UINT64_C(18446744073709551527),
    UINT64_C(18446744073709551529),
    UINT64_C(18446744073709551530),
    UINT64_C(18446744073709551531)
};

#if HURCHALLA_COMPILER_HAS_UINT128_T()
// rely on 2^128 wrap-around in the subtractions below
static __uint128_t prime_numbers128[] = {
    static_cast<__uint128_t>(0) - 159,
    static_cast<__uint128_t>(0) - 173,
    static_cast<__uint128_t>(0) - 233
};
static __uint128_t nonprime_numbers128[] = {
    static_cast<__uint128_t>(0) - 160,
    static_cast<__uint128_t>(0) - 161,
    static_cast<__uint128_t>(0) - 163,
    static_cast<__uint128_t>(0) - 165,
    static_cast<__uint128_t>(0) - 167
};
#endif


#endif
