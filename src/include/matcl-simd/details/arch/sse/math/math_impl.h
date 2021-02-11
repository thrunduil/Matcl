/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/math/simd_math_func_def.h"

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                          exp
//-----------------------------------------------------------------------
template<>
struct simd_exp<double, 256, sse_tag>
{
    using simd_type     = simd<double, 256, sse_tag>;
    using simd_half     = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        simd_half v1    = exp(a.extract_low());
        simd_half v2    = exp(a.extract_high());

        return simd_type(v1, v2);
    };
};

template<>
struct simd_exp<float, 256, sse_tag>
{
    using simd_type     = simd<float, 256, sse_tag>;
    using simd_half     = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        simd_half v1    = exp(a.extract_low());
        simd_half v2    = exp(a.extract_high());

        return simd_type(v1, v2);
    };
};

//-----------------------------------------------------------------------
//                          log
//-----------------------------------------------------------------------
template<>
struct simd_log<double, 256, sse_tag>
{
    using simd_type     = simd<double, 256, sse_tag>;
    using simd_half     = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        simd_half v1    = log(a.extract_low());
        simd_half v2    = log(a.extract_high());

        return simd_type(v1, v2);
    };
};

template<>
struct simd_log<float, 256, sse_tag>
{
    using simd_type     = simd<float, 256, sse_tag>;
    using simd_half     = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        simd_half v1    = log(a.extract_low());
        simd_half v2    = log(a.extract_high());

        return simd_type(v1, v2);
    };
};

//-----------------------------------------------------------------------
//                          pow2k
//-----------------------------------------------------------------------
template<class T>
struct simd_pow2k<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;
    using int_type  = typename details::integer_type<T>::type;
    using simd_int  = simd<int_type, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_type(pow2k(k.extract_low()), pow2k(k.extract_high()));
    };

    force_inline
    static simd_type eval_i(const simd_int& k)
    {
        return simd_type(pow2ki(k.extract_low()), pow2ki(k.extract_high()));
    };
};

template<class T>
struct simd_pow2k<T, 128, sse_tag>
{
    using simd_type = simd<T, 128, sse_tag>;
    using int_type  = typename details::integer_type<T>::type;
    using simd_int  = simd<int_type, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_pow2k_impl<T, 128, sse_tag>::eval(k);
    };

    force_inline
    static simd_type eval_i(const simd_int& k)
    {
        return simd_pow2k_impl<T, 128, sse_tag>::eval_i(k);
    };
};

template<class T>
struct simd_pow2k<T, 128, scalar_sse_tag>
{
    using simd_type     = simd<T, 128, scalar_sse_tag>;
    using int_type      = typename details::integer_type<T>::type;
    using simd_int      = simd<int_type, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_type(pow2k(k.as_vector()));
    };

    force_inline
    static simd_type eval_i(const simd_int& k)
    {
        return simd_type(pow2ki(k.as_vector()));
    };
};

//-----------------------------------------------------------------------
//                          exponent
//-----------------------------------------------------------------------
template<class T>
struct simd_exponent<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;
    using int_type  = typename details::integer_type<T>::type;
    using simd_int  = simd<int_type, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_type(exponent(k.extract_low()), exponent(k.extract_high()));
    };

    force_inline
    static simd_int eval_i(const simd_type& k)
    {
        return simd_int(iexponent(k.extract_low()), iexponent(k.extract_high()));
    };
};

template<class T>
struct simd_exponent<T, 128, sse_tag>
{
    using simd_type = simd<T, 128, sse_tag>;
    using int_type  = typename details::integer_type<T>::type;
    using simd_int  = simd<int_type, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_exponent_impl<T, 128, sse_tag>::eval(k);
    };

    force_inline
    static simd_int eval_i(const simd_type& k)
    {
        return simd_exponent_impl<T, 128, sse_tag>::eval_i(k);
    };
};

template<class T>
struct simd_exponent<T, 128, scalar_sse_tag>
{
    using simd_type     = simd<T, 128, scalar_sse_tag>;
    using int_type      = typename details::integer_type<T>::type;
    using simd_int      = simd<int_type, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_type(exponent(k.as_vector()));
    };

    force_inline
    static simd_int eval_i(const simd_type& k)
    {
        return simd_int(iexponent(k.as_vector()));
    };
};

//-----------------------------------------------------------------------
//                          fraction
//-----------------------------------------------------------------------
template<class T>
struct simd_fraction<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;
    using int_type  = typename details::integer_type<T>::type;
    using simd_int  = simd<int_type, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_type(fraction(k.extract_low()), fraction(k.extract_high()));
    };
};

template<class T>
struct simd_fraction<T, 128, sse_tag>
{
    using simd_type = simd<T, 128, sse_tag>;
    using int_type  = typename details::integer_type<T>::type;
    using simd_int  = simd<int_type, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_fraction_impl<T, 128, sse_tag>::eval(k);
    };
};

template<class T>
struct simd_fraction<T, 128, scalar_sse_tag>
{
    using simd_type     = simd<T, 128, scalar_sse_tag>;
    using int_type      = typename details::integer_type<T>::type;
    using simd_int      = simd<int_type, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_type(fraction(k.as_vector()));
    };
};

//-----------------------------------------------------------------------
//                          copysign
//-----------------------------------------------------------------------
template<class T>
struct simd_copysign<T, 256, sse_tag>
{
    using simd_type = simd<T, 256, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(copysign(x.extract_low(), y.extract_low()), 
                        copysign(x.extract_high(), y.extract_high()));
    };
};

template<class T>
struct simd_copysign<T, 128, sse_tag>
{
    using simd_type = simd<T, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_copysign_impl<T, 128, sse_tag>::eval(x, y);
    };
};

template<class T>
struct simd_copysign<T, 128, scalar_sse_tag>
{
    using simd_type     = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(copysign(x.as_vector(), y.as_vector()));
    };
};

//-----------------------------------------------------------------------
//                          sin
//-----------------------------------------------------------------------
template<>
struct simd_sincos<double, 256, sse_tag>
{
    using simd_type     = simd<double, 256, sse_tag>;
    using simd_half     = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval_sin(const simd_type& a)
    {
        simd_half v1    = sin(a.extract_low());
        simd_half v2    = sin(a.extract_high());

        return simd_type(v1, v2);
    };

    force_inline
    static simd_type eval_cos(const simd_type& a)
    {
        simd_half v1    = cos(a.extract_low());
        simd_half v2    = cos(a.extract_high());

        return simd_type(v1, v2);
    };
};

template<>
struct simd_sincos<float, 256, sse_tag>
{
    using simd_type     = simd<float, 256, sse_tag>;
    using simd_half     = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval_sin(const simd_type& a)
    {
        simd_half v1    = sin(a.extract_low());
        simd_half v2    = sin(a.extract_high());

        return simd_type(v1, v2);
    };

    force_inline
    static simd_type eval_cos(const simd_type& a)
    {
        simd_half v1    = cos(a.extract_low());
        simd_half v2    = cos(a.extract_high());

        return simd_type(v1, v2);
    };
};

//-----------------------------------------------------------------------
//                          tan/cot
//-----------------------------------------------------------------------
template<>
struct simd_tancot<double, 256, sse_tag>
{
    using simd_type     = simd<double, 256, sse_tag>;
    using simd_half     = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval_tan(const simd_type& a)
    {
        simd_half v1    = tan(a.extract_low());
        simd_half v2    = tan(a.extract_high());

        return simd_type(v1, v2);
    };

    force_inline
    static simd_type eval_cot(const simd_type& a)
    {
        simd_half v1    = cot(a.extract_low());
        simd_half v2    = cot(a.extract_high());

        return simd_type(v1, v2);
    };

};

template<>
struct simd_tancot<float, 256, sse_tag>
{
    using simd_type     = simd<float, 256, sse_tag>;
    using simd_half     = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval_tan(const simd_type& a)
    {
        simd_half v1    = tan(a.extract_low());
        simd_half v2    = tan(a.extract_high());

        return simd_type(v1, v2);
    };

    force_inline
    static simd_type eval_cot(const simd_type& a)
    {
        simd_half v1    = cot(a.extract_low());
        simd_half v2    = cot(a.extract_high());

        return simd_type(v1, v2);
    };
};

}}}
