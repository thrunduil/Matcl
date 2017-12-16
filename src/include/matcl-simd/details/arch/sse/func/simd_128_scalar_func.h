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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/func/simd_func_def.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace simd { namespace details
{

template<class T>
struct simd_reverse<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<class T>
struct simd_mult<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(mult(x.as_vector(), y.as_vector()));
    };
};

template<>
struct simd_mult<int32_t, 128, scalar_sse_tag>
{
    using simd_type = simd<int32_t, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_mullo_epi32(x.data, y.data);
        #else
            return simd_type(x.first() * y.first());
        #endif        
    };
};

template<>
struct simd_mult<int64_t, 128, scalar_sse_tag>
{
    using simd_type = simd<int64_t, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.first() * y.first());
    };
};

template<class T>
struct simd_div<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(div(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_plus<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(plus(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_minus<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(minus(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_uminus<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(uminus(x.as_vector()));
    };
};

template<class T>
struct simd_sum_all<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static T eval(const simd_type& x)
    {
        return x.first();
    };
};

template<class T>
struct simd_fma_f<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fma_f(x.as_vector(), y.as_vector(), z.as_vector()));
    };
};

template<class T>
struct simd_fms_f<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fms_f(x.as_vector(), y.as_vector(), z.as_vector()));
    };
};

//
template<class T>
struct simd_fnma_f<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fnma_f(x.as_vector(), y.as_vector(), z.as_vector()));
    };
};

template<class T>
struct simd_fnms_f<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fnms_f(x.as_vector(), y.as_vector(), z.as_vector()));
    };
};

template<class T>
struct simd_fma_a<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fma_a(x.as_vector(), y.as_vector(), z.as_vector()));
    };
};

template<class T>
struct simd_fms_a<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fms_a(x.as_vector(), y.as_vector(), z.as_vector()));
    };
};

//
template<class T>
struct simd_fnma_a<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fnma_a(x.as_vector(), y.as_vector(), z.as_vector()));
    };
};

template<class T>
struct simd_fnms_a<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(fnms_a(x.as_vector(), y.as_vector(), z.as_vector()));
    };
};

template<class T>
struct simd_abs<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(abs(x.as_vector()));
    };
};

template<class T>
struct simd_bitwise_or<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(bitwise_or(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_bitwise_xor<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(bitwise_xor(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_bitwise_and<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(bitwise_and(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_bitwise_andnot<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(bitwise_andnot(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_bitwise_not<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(bitwise_not(x.as_vector()));
    };
};

template<class T>
struct simd_shift_left<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return simd_type(shift_left(x.as_vector(), y));
    };
};

template<class T>
struct simd_shift_right<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return simd_type(shift_right(x.as_vector(), y));
    };
};

template<class T>
struct simd_shift_right_arithmetic<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return simd_type(shift_right_arithmetic(x.as_vector(), y));
    };
};

template<>
struct simd_shift_right_arithmetic<int64_t, 128, scalar_sse_tag>
{
    using simd_type = simd<int64_t, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return simd_type(scalar_func::shift_right_arithmetic(x.first(), y));
    };
};

template<class T>
struct simd_max<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(max(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_min<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(min(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_sqrt<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(sqrt(x.as_vector()));
    };
};

template<class T>
struct simd_round<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return simd_type(round(x.as_vector()));
        #else
            auto s1   = scalar_func::round(x.first());
            return simd_type(s1);
        #endif
    };
};

template<class T>
struct simd_floor<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return simd_type(floor(x.as_vector()));
        #else
            auto s1   = scalar_func::floor(x.first());
            return simd_type(s1);
        #endif
    };
};

template<class T>
struct simd_ceil<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return simd_type(ceil(x.as_vector()));
        #else
            auto s1   = scalar_func::ceil(x.first());
            return simd_type(s1);
        #endif
    };
};

template<class T>
struct simd_trunc<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return simd_type(trunc(x.as_vector()));
        #else
            auto s1   = scalar_func::trunc(x.first());
            return simd_type(s1);
        #endif
    };
};

template<class T>
struct simd_eeq<T, 128, scalar_sse_tag> : simd_bool_base<double>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(eeq(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_neq<T, 128, scalar_sse_tag> : simd_bool_base<double>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(neq(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_lt<T, 128, scalar_sse_tag> : simd_bool_base<double>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(lt(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_gt<T, 128, scalar_sse_tag> : simd_bool_base<double>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(gt(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_leq<T, 128, scalar_sse_tag> : simd_bool_base<double>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(leq(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_geq<T, 128, scalar_sse_tag> : simd_bool_base<double>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(geq(x.as_vector(), y.as_vector()));
    };
};

template<class T>
struct simd_any_nan<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        auto res = is_nan(x.as_vector());
        return res.first() != 0;
    };
};

template<class T>
struct simd_is_nan<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(is_nan(x.as_vector()));
    };
};

template<class T>
struct simd_any<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        return x.first() != 0;
    };
};

template<class T>
struct simd_all<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        return x.first() != 0;
    };
};

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<class T>
struct simd_if_then_else<T, 128, scalar_sse_tag>
{
    using simd_type = simd<T, 128, scalar_sse_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return simd_type(if_then_else(test.as_vector(), val_true.as_vector(), 
                                          val_false.as_vector()));
        #else
            return (test.first() == T()) ? val_false : val_true;
        #endif  
    };
};

}}}
