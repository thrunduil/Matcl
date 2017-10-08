/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#pragma once

#include "matcl-mp/matcl_mp.h"
#include "matcl-scalar/matcl_scalar.h"

namespace matcl { namespace test
{

#define DEFINE_TEST_FUNC(name)                              \
struct Func_##name                                          \
{                                                           \
    static std::string func_name() { return #name; };       \
                                                            \
    template<class T1>                                      \
    static auto eval_matcl(const T1& a1)                    \
        -> decltype(matcl::name(std::declval<T1>()))        \
    {                                                       \
        return matcl::name(a1);                             \
    };                                                      \
                                                            \
    template<class T1>                                      \
    static auto eval_std(const T1& a1)                      \
        -> decltype(std::name(std::declval<T1>()))          \
    {                                                       \
        return std::name(a1);                               \
    };                                                      \
                                                            \
    template<class T1>                                      \
    static auto eval_mp(const T1& a1, precision p)          \
        -> decltype(name(std::declval<T1>(), p))            \
    {                                                       \
        return name(a1, p);                                 \
    };                                                      \
};                                                          \

DEFINE_TEST_FUNC(sqrt)
DEFINE_TEST_FUNC(sqrt1pm1)
DEFINE_TEST_FUNC(cbrt)
DEFINE_TEST_FUNC(exp)
DEFINE_TEST_FUNC(expi)
DEFINE_TEST_FUNC(exp2)
DEFINE_TEST_FUNC(exp10)
DEFINE_TEST_FUNC(expm1)
DEFINE_TEST_FUNC(log)
DEFINE_TEST_FUNC(log2)
DEFINE_TEST_FUNC(log10)
DEFINE_TEST_FUNC(log1p)

DEFINE_TEST_FUNC(sin)
DEFINE_TEST_FUNC(cos)
DEFINE_TEST_FUNC(tan)
DEFINE_TEST_FUNC(cot)
DEFINE_TEST_FUNC(sec)
DEFINE_TEST_FUNC(csc)

DEFINE_TEST_FUNC(sinh)
DEFINE_TEST_FUNC(cosh)
DEFINE_TEST_FUNC(tanh)
DEFINE_TEST_FUNC(coth)
DEFINE_TEST_FUNC(sech)
DEFINE_TEST_FUNC(csch)

DEFINE_TEST_FUNC(asin)
DEFINE_TEST_FUNC(acos)
DEFINE_TEST_FUNC(atan)
DEFINE_TEST_FUNC(acot)
DEFINE_TEST_FUNC(asec)
DEFINE_TEST_FUNC(acsc)

DEFINE_TEST_FUNC(asinh)
DEFINE_TEST_FUNC(acosh)
DEFINE_TEST_FUNC(atanh)
DEFINE_TEST_FUNC(acoth)
DEFINE_TEST_FUNC(asech)
DEFINE_TEST_FUNC(acsch)

DEFINE_TEST_FUNC(uminus)
DEFINE_TEST_FUNC(abs)
DEFINE_TEST_FUNC(abs2)
DEFINE_TEST_FUNC(arg)
DEFINE_TEST_FUNC(inv)

DEFINE_TEST_FUNC(real)
DEFINE_TEST_FUNC(imag)
DEFINE_TEST_FUNC(conj)
DEFINE_TEST_FUNC(sign)

}};

namespace matcl
{
    // 0.5 ulp error implementations of these functions
    // are not available; just calculate then according to definition
    // with increased precision

    template<class MP_ty>
    inline MP_ty sqrt1pm1(const MP_ty& x, precision p)
    {
        precision p2    = p + 50;

        MP_ty x2(x, p2);        
        return MP_ty(x2 / (sqrt(1.0 + x2) + 1.0), p);
    };

    template<class MP_ty>
    inline MP_ty acot(const MP_ty& x, precision p)
    {
        precision p2    = p + 50;

        MP_ty x2(x, p2);
        return atan(1.0/x2, p);
    };

    template<class MP_ty>
    inline MP_ty asec(const MP_ty& x, precision p)
    {
        precision p2    = p + 50;

        MP_ty x2(x, p2);
        return acos(1.0/x2, p);
    };

    template<class MP_ty>
    inline MP_ty acsc(const MP_ty& x, precision p)
    {
        precision p2    = p + 50;

        MP_ty x2(x, p2);
        return asin(1.0/x2, p);
    };

    template<class MP_ty>
    inline MP_ty acoth(const MP_ty& x, precision p)
    {
        precision p2    = p + 50;

        MP_ty x2(x, p2);
        return atanh(1.0/x2, p);
    };

    template<class MP_ty>
    inline MP_ty asech(const MP_ty& x, precision p)
    {
        precision p2    = p + 50;

        MP_ty x2(x, p2);
        return acosh(1.0/x2, p);
    };

    template<class MP_ty>
    inline MP_ty acsch(const MP_ty& x, precision p)
    {
        precision p2    = p + 50;

        MP_ty x2(x, p2);
        return asinh(1.0/x2, p);
    };

    inline mp_float conj(const mp_float& x, precision)
    {
        return conj(x);
    }

    inline mp_float real(const mp_float& x, precision)
    {
        return real(x);
    }

    inline mp_float imag(const mp_float& x, precision)
    {
        return imag(x);
    }

    inline mp_complex conj(const mp_complex& x, precision)
    {
        return conj(x);
    }

    inline mp_float real(const mp_complex& x, precision)
    {
        return real(x);
    }

    inline mp_float imag(const mp_complex& x, precision)
    {
        return imag(x);
    }
};

//add missing functions to std
namespace std
{    
    template<class Ty>
    inline Ty sqrt1pm1(Ty val)  { return matcl::sqrt1pm1(val); };

    template<class Ty>
    inline Ty exp10(Ty val)     { return std::pow(Ty(10.0), val); };

    template<class Ty>
    inline Ty cot(Ty val)       { return matcl::cot(val); };

    template<class Ty>
    inline Ty sec(Ty val)       { return matcl::sec(val); };

    template<class Ty>
    inline Ty csc(Ty val)       { return matcl::csc(val); };

    template<class Ty>
    inline Ty coth(Ty val)      { return matcl::coth(val); };

    template<class Ty>
    inline Ty sech(Ty val)      { return matcl::sech(val); };

    template<class Ty>
    inline Ty csch(Ty val)      { return matcl::csch(val); };

    template<class Ty>
    inline Ty acot(Ty val)      { return matcl::acot(val); };

    template<class Ty>
    inline Ty asec(Ty val)      { return matcl::asec(val); };

    template<class Ty>
    inline Ty acsc(Ty val)      { return matcl::acsc(val); };

    template<class Ty>
    inline Ty acoth(Ty val)     { return matcl::acoth(val); };

    template<class Ty>
    inline Ty asech(Ty val)     { return matcl::asech(val); };

    template<class Ty>
    inline Ty acsch(Ty val)     { return matcl::acsch(val); };

    template<class Ty>
    inline Ty inv(Ty val)       { return Ty(1.0) / val; };

    template<class Ty>
    inline Ty abs2(Ty val)      { return val * val; };

    template<class Ty>
    inline Ty uminus(Ty val)    { return -val; };

    // std::arg is broken
    inline float arg(float val)     { return matcl::arg(val); };
    inline double arg(double val)   { return matcl::arg(val); };

    template<class Ty>
    inline std::complex<Ty> log1p(const std::complex<Ty>& val)
    {
        return matcl::log1p(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> log2(const std::complex<Ty>& val)
    {
        return matcl::log2(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> expm1(const std::complex<Ty>& val)
    {
        return matcl::expm1(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> exp2(const std::complex<Ty>& val)
    {
        return matcl::exp2(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> sign(const std::complex<Ty>& val)
    {
        return matcl::sign(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> sqrt1pm1(const std::complex<Ty>& val)
    {
        return matcl::sqrt1pm1(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> exp10(const std::complex<Ty>& val)
    {
        return matcl::exp10(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> cot(const std::complex<Ty>& val)
    {
        return matcl::cot(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> sec(const std::complex<Ty>& val)
    {
        return matcl::sec(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> csc(const std::complex<Ty>& val)
    {
        return matcl::csc(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> coth(const std::complex<Ty>& val)
    {
        return matcl::coth(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> sech(const std::complex<Ty>& val)
    {
        return matcl::sech(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> csch(const std::complex<Ty>& val)
    {
        return matcl::csch(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> acot(const std::complex<Ty>& val)
    {
        return matcl::acot(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> asec(const std::complex<Ty>& val)
    {
        return matcl::asec(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> acsc(const std::complex<Ty>& val)
    {
        return matcl::acsc(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> acoth(const std::complex<Ty>& val)
    {
        return matcl::acoth(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> asech(const std::complex<Ty>& val)
    {
        return matcl::asech(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> acsch(const std::complex<Ty>& val)
    {
        return matcl::acsch(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline std::complex<Ty> inv(const std::complex<Ty>& val)
    {
        return Ty(1.0)/val;
    }

    template<class Ty>
    inline std::complex<Ty> abs2(const std::complex<Ty>& val)
    {
        return val*val;
    }

    template<class Ty>
    inline std::complex<Ty> uminus(const std::complex<Ty>& val)
    {
        return -val;
    }

    template<class Ty>
    inline std::complex<Ty> expi(const std::complex<Ty>& val)
    {
        return matcl::expi(matcl::complex<Ty>(val)).value;
    }

    template<class Ty>
    inline Ty sign(Ty val)      { return matcl::sign(val); };
}