// Copyright (c) 2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_DUAL_MONTGOMERY_FORM_H_INCLUDED
#define HURCHALLA_DUAL_MONTGOMERY_FORM_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>


namespace hurchalla {


template <class MF>
class DualMontgomeryForm {
    const MF mf;
    using T = typename MF::IntegerType;

    class V {
        friend DualMontgomeryForm<MF>;
        struct OpenMV : public MF::MontgomeryValue {
            using MF::MontgomeryValue::get;
            OpenMV(typename MF::MontgomeryValue mfmv) : MF::MontgomeryValue(mfmv) {}
            OpenMV(T a) : MF::MontgomeryValue(a) {}
            OpenMV() {}
        };
        OpenMV v1;
        OpenMV v2;
    protected:
        HURCHALLA_FORCE_INLINE
        explicit V(T a, T b) : v1(a), v2(b) {}
    public:
        HURCHALLA_FORCE_INLINE V() = default;
    };
    class C {
        friend DualMontgomeryForm<MF>;
        struct OpenCV : public MF::CanonicalValue {
            using MF::CanonicalValue::get;
            OpenCV(typename MF::CanonicalValue mfcv) : MF::CanonicalValue(mfcv) {}
            OpenCV(T a) : MF::CanonicalValue(a) {}
            OpenCV() {}
        };
        OpenCV c1;
        OpenCV c2;
    protected:
        HURCHALLA_FORCE_INLINE
        explicit C(T a, T b) : c1(a), c2(b) {}
    public:
        HURCHALLA_FORCE_INLINE C() = default;
        HURCHALLA_FORCE_INLINE operator V() const
        {
            return V(c1.get(), c2.get());
        }
    };
public:
    using IntegerType = T;
    using CanonicalValue = C;
    using MontgomeryValue = V;

    DualMontgomeryForm(T modulus) : mf(modulus) {}
    DualMontgomeryForm(const MF& mfc) : mf(mfc) {}

    HURCHALLA_FORCE_INLINE T getModulus() const
    {
        return mf.getModulus();
    }
    HURCHALLA_FORCE_INLINE V convertIn(T a, T b) const
    {
        V x;
        x.v1 = mf.convertIn(a);
        x.v2 = mf.convertIn(b);
        return x;
    }
    HURCHALLA_FORCE_INLINE void convertOut(T& a, T& b, V x) const
    {
        a = mf.convertOut(x.v1);
        b = mf.convertOut(x.v2);
    }

    HURCHALLA_FORCE_INLINE C getUnityValue() const
    {
        C c;
        c.c1 = mf.getUnityValue();
        c.c2 = mf.getUnityValue();
        return c;
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        C c;
        c.c1 = mf.getCanonicalValue(x.v1);
        c.c2 = mf.getCanonicalValue(x.v2);
        return c;
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        V z;
        z.v1 = mf.add(x.v1, y.v1);
        z.v2 = mf.add(x.v2, y.v2);
        return z;
    }
    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        V z;
        z.v1 = mf.subtract(x.v1, y.v1);
        z.v2 = mf.subtract(x.v2, y.v2);
        return z;
    }

    HURCHALLA_FORCE_INLINE C add(C x, C y) const
    {
        C z;
        z.c1 = mf.add(x.c1, y.c1);
        z.c2 = mf.add(x.c2, y.c2);
        return z;
    }
    HURCHALLA_FORCE_INLINE C subtract(C x, C y) const
    {
        C z;
        z.c1 = mf.subtract(x.c1, y.c1);
        z.c2 = mf.subtract(x.c2, y.c2);
        return z;
    }

    template <class PTAG = hurchalla::LowlatencyTag>
    HURCHALLA_FORCE_INLINE V square(V x) const
    {
        V z;
        z.v1 = mf.template square<PTAG>(x.v1);
        z.v2 = mf.template square<PTAG>(x.v2);
        return z;
    }

    template <class PTAG = hurchalla::LowlatencyTag>
    HURCHALLA_FORCE_INLINE V multiply(V x, V y) const
    {
        bool isZero;
        V z;
        z.v1 = mf.template multiply<PTAG>(x.v1, y.v1, isZero);
        z.v2 = mf.template multiply<PTAG>(x.v2, y.v2, isZero);
        return z;
    }
    template <class PTAG = hurchalla::LowlatencyTag>
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, bool& isZero) const
    {
        bool isZero1, isZero2;
        V z;
        z.v1 = mf.template multiply<PTAG>(x.v1, y.v1, isZero1);
        z.v2 = mf.template multiply<PTAG>(x.v2, y.v2, isZero2);
        // this isZero assignment is a slight hack, but will work
        // out ok for ECM.  The problem is the multiply() API
        // sets exactly one isZero variable, and we really want
        // to reuse the API without change.
        isZero = isZero1 || isZero2;
        return z;
    }

    template <class PTAG = hurchalla::LowlatencyTag>
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, C z) const
    {
        V p;
        p.v1 = mf.template fmadd<PTAG>(x.v1, y.v1, z.c1);
        p.v2 = mf.template fmadd<PTAG>(x.v2, y.v2, z.c2);
        return p;
    }

    template <class PTAG = hurchalla::LowlatencyTag>
    HURCHALLA_FORCE_INLINE V fusedSquareSub(V x, C z) const
    {
        V p;
        p.v1 = mf.template fusedSquareSub<PTAG>(x.v1, z.c1);
        p.v2 = mf.template fusedSquareSub<PTAG>(x.v2, z.c2);
        return p;
    }
    template <class F> HURCHALLA_FORCE_INLINE
    T gcd_with_modulus(V x, const F& gcd_functor) const
    {
        // this is designed to extract a factor, and stretches
        // the meaning of finding the gcd (since we combine
        // two different values before taking the gcd with the
        // modulus).
        bool isZero;
        auto tmp = mf.template multiply(x.v1, x.v2, isZero);
        if (isZero) {
            if (mf.getCanonicalValue(x.v1) == mf.getZeroValue())
                tmp = x.v2;
            else
                tmp = x.v1;
        }
        return static_cast<T>(mf.gcd_with_modulus(tmp, gcd_functor));
    }

  // support functions intended for ECM, to help calculate the inverse
    template <class PTAG = hurchalla::LowlatencyTag>
    HURCHALLA_FORCE_INLINE T crossMultiplyAndConvertOut(V x) const
    {
        bool isZero;
        auto y = mf.template multiply<PTAG>(x.v1, x.v2, isZero);
        T a = mf.convertOut(y);
        return a;
    }
    HURCHALLA_FORCE_INLINE V convertInAndCopy(T a) const
    {
        V x;
        x.v1 = mf.convertIn(a);
        x.v2 = x.v1;
        return x;
    }
    HURCHALLA_FORCE_INLINE V swapChannels(V x) const
    {
        V y;
        y.v1 = x.v2;
        y.v2 = x.v1;
        return y;
    }
};


} // end namespace

#endif
