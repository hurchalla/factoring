// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_FACTORING_FACTORS_CONTAINER_ADAPTER_H_INCLUDED
#define HURCHALLA_FACTORING_FACTORS_CONTAINER_ADAPTER_H_INCLUDED


#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstddef>
#include <vector>
#include <array>
#include <limits>

namespace hurchalla { namespace factoring {


// This #if is set to 1, at least for now, since most likely there is no benefit
// to FactorsContainerAdapter supporting any containers beyond std::vector and
// std::array.  It can be set to 0 (and should work), if gaining compatibility
// with almost any container becomes more important than the resultant potential
// for misuse with some unknown and subtly incompatible container choice T.
#if 1
// Primary template
template <typename T> class FactorsContainerAdapter {
    static_assert(sizeof(T) < 0, "T must be std::vector or std::array");
};
#else
namespace fca_detail {
    // SFINAE: this function is enabled if reserve() exists for the container.
    template <class C>
    auto reserve_size(C& container, std::size_t num_elements)
        -> decltype(container.reserve(num_elements))
    {
        container.reserve(num_elements);
    }

    // The param "..." is a weaker match than size_t (above) for matching, which
    // makes this function a default fall-back when reserve() does not exist.
    template <class C>
    void reserve_size(C& container, ...) {}
}

// This primary template is a fall-back class that should work for most
// containers - though it's questionable if this compatibility has any benefit.
template <typename T>
class FactorsContainerAdapter {
    T& c_;
public:
    using value_type = typename T::value_type;
    explicit FactorsContainerAdapter(T& c) : c_(c) {};
    void reserve(std::size_t num_elements)
    {
        fca_detail::reserve_size(c_, num_elements);
    }
    std::size_t size() const { return c_.size(); }
    void push(value_type val) { c_.insert(c_.end(), val); }
};
#endif


template <typename T, std::size_t N>
class FactorsContainerAdapter<std::array<T, N>> {
    // std::numeric_limits<T>::digits, in general, is the minimum array length
    // needed to guarantee all factors of a type T number will fit in the array.
    // Note that the longest possible factor sequence for type T that can occur
    // is when the number to factor is the largest possible power of 2 that can
    // fit in T, which results in all 2s for factors.
    static_assert(N == std::numeric_limits<T>::digits, "");

    std::array<T, N>& c_;
    std::size_t num_factors_ = 0;
public:
    using value_type = T;
    explicit FactorsContainerAdapter(std::array<T, N>& c) : c_(c) {}
    void reserve(std::size_t) {}
    std::size_t size() const
    {
        return num_factors_;
    }
    void push(value_type val)
    {
        // given this class's static_assert, it should be impossible for the
        // assert here to fail, so long as we push only factors.
        HPBC_PRECONDITION2(num_factors_ < N);
        c_[num_factors_] = val;
        ++num_factors_;
    }
};


template <typename... Types>
class FactorsContainerAdapter<std::vector<Types...>> {
    std::vector<Types...>& c_;
public:
    using value_type = typename std::vector<Types...>::value_type;
    explicit FactorsContainerAdapter(std::vector<Types...>& c) : c_(c) {}
    void reserve(std::size_t num_elements)
    {
        c_.reserve(num_elements);
    }
    std::size_t size() const
    {
        return c_.size();
    }
    void push(value_type val)
    {
        c_.push_back(val);
    }
};


}}  // end namespace

#endif
