// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_FACTORING_FACTOR_BY_TABLE32_H_INCLUDED
#define HURCHALLA_FACTORING_FACTOR_BY_TABLE32_H_INCLUDED

#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/factoring/detail/ImplFactorByTable32.h"
#include <cstdint>
#include <array>

namespace hurchalla {


// Most commonly you would use FactorByTable32 (defined at bottom),
// rather than declaring a variable with this class template.

// Instantiating this template class will (indirectly) create a factor table
// with 2^(INPUT_BIT_LIMIT) entries.  If you specify an INPUT_BIT_LIMIT lower
// than the default 32, it will save memory and improve performance.  You are
// always strictly limited to factoring numbers below 2^(INPUT_BIT_LIMIT).
//
// Instantiating with FAVOR_SMALL_SIZE = false (the default) will usually
// result in faster factoring than if you choose true.  Initial tests suggest
// it tends to be about 10% faster on x86.  However, if you set FAVOR_SMALL_SIZE
// to true, you save ~15% memory space for the table.  You will probably want
// to use writeTableToFile() at some point, and if so, FAVOR_SMALL_SIZE would
// similarly affect your serialized data file's size.
// (Assuming INPUT_BIT_LIMIT = 32, with FAVOR_SMALL_SIZE = true you can expect
// memory/disk usage of ~1.4GB, and with FAVOR_SMALL_SIZE = false you can expect
// memory/disk usage of ~1.6GB.)
template <int INPUT_BIT_LIMIT = 32, bool FAVOR_SMALL_SIZE = false>
class FactorByTable
{
    static_assert(0 < INPUT_BIT_LIMIT && INPUT_BIT_LIMIT <= 32);
    detail::ImplFactorByTable32<INPUT_BIT_LIMIT, FAVOR_SMALL_SIZE> impl;

public:
    FactorByTable(const FactorByTable&) = delete;
    FactorByTable(FactorByTable&& other) : impl(std::move(other.impl)) {}

    // This might take a few minutes to construct, since it creates a ~1.5GB
    // table in memory from scratch.
    FactorByTable() : impl() {}

    // If you call this next constructor with an argument of false for
    // createTableIfCantOpen, then the constructor will throw if it is unable
    // to open table_filepath.  Otherwise (and by default), it will create the
    // table, which will very likely take a few minutes to complete.
    //
    // Can throw from a file open failure, a read failure, or mismatch in
    // file values read vs values expected.
    FactorByTable(const char* table_filepath, bool createTableIfCantOpen = true)
        : impl(table_filepath, createTableIfCantOpen) {}

    // Can throw from a file open failure or write failure.
    void writeTableToFile(const char* table_filepath) const
    {
        impl.writeTableToFile(table_filepath);
    }

    // the factorize function
    std::array<std::uint32_t, 32>
    operator()(std::uint32_t x, unsigned int& num_factors) const
    {
        // 0 and 1 do not have prime factorizations
        HPBC_PRECONDITION2(x >= 2);
        // we can only factor numbers below 2^(INPUT_BIT_LIMIT)
        HPBC_PRECONDITION2(x < (static_cast<uint64_t>(1) << INPUT_BIT_LIMIT));

        return impl.factorize(x, num_factors);
    }
};

using FactorByTable32 = FactorByTable<32, false>;


} // end namespace

#endif
