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

#include "matcl-simd/arch/simd_impl.h"
#include "matcl-simd/func/simd_func_def.h"
#include "matcl-simd/details/scalar_func.h"

namespace matcl { namespace simd
{

template<>
struct simd_reverse<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_shuffle_ps(x.data, x.data, _MM_SHUFFLE(0,1,2,3));
    };
};

template<>
struct simd_mult<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_mul_ps( x.data, y.data );
    };
};

template<>
struct simd_div<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_div_ps( x.data, y.data );
    };
};

template<>
struct simd_plus<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_add_ps( x.data, y.data );
    };
};

template<>
struct simd_minus<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_sub_ps( x.data, y.data );
    };
};

template<>
struct simd_uminus<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_xor_ps( x.data, _mm_set1_ps(-0.0f) );
    };
};

template<>
struct simd_sum_all<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static float eval(const simd_type& x)
    {
       #if MATCL_ARCHITECTURE_HAS_SSE3
            __m128 s    = _mm_hadd_ps(x.data, x.data);
            s           = _mm_hadd_ps(s, s);
            return s.m128_f32[0];
        #else
            float s1    = x.get<0>() + x.get<1>();
            float s2    = x.get<2>() + x.get<3>();
            return s1 + s2;
        #endif
    };
};

template<>
struct simd_sub_add<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE3
            __m128 s   = _mm_addsub_ps(x.data, y.data);
            return s;
        #else
            float s1    = x.get<0>() - y.get<0>();
            float s2    = x.get<1>() + y.get<1>();
            float s3    = x.get<2>() - y.get<2>();
            float s4    = x.get<3>() + y.get<3>();
            return simd_type(s1, s2, s3, s4);
        #endif
    };
};
template<>
struct simd_fma_f<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fmadd_ps( x.data, y.data, z.data);
        #else
            return x * y + z;
        #endif
    };
};

template<>
struct simd_fms_f<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fmsub_ps( x.data, y.data, z.data);
        #else
            return x * y - z;
        #endif
    };
};

template<>
struct simd_abs<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const __m128 sign_mask = _mm_set1_ps(-0.f); // -0.f = 1 << 31
        return _mm_andnot_ps(sign_mask, x.data);
    };
};

template<>
struct simd_max<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_max_ps(x.data, y.data);
    };
};

template<>
struct simd_min<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_min_ps(x.data, y.data);
    };
};

template<>
struct simd_sqrt<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_sqrt_ps(x.data);
    };
};

template<>
struct simd_round<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_round_ps(x.data, _MM_FROUND_NINT);
        #else
            float s1    = scalar_func::round(x.get<0>());
            float s2    = scalar_func::round(x.get<1>());
            float s3    = scalar_func::round(x.get<2>());
            float s4    = scalar_func::round(x.get<3>());

            return simd_type(s1, s2, s3, s4);
        #endif
    };
};

template<>
struct simd_floor<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_round_ps(x.data, _MM_FROUND_FLOOR);
        #else
            float s1    = scalar_func::floor(x.get<0>());
            float s2    = scalar_func::floor(x.get<1>());
            float s3    = scalar_func::floor(x.get<2>());
            float s4    = scalar_func::floor(x.get<3>());

            return simd_type(s1, s2, s3, s4);
        #endif
    };
};

template<>
struct simd_ceil<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_round_ps(x.data, _MM_FROUND_CEIL);
        #else
            float s1    = scalar_func::ceil(x.get<0>());
            float s2    = scalar_func::ceil(x.get<1>());
            float s3    = scalar_func::ceil(x.get<2>());
            float s4    = scalar_func::ceil(x.get<3>());

            return simd_type(s1, s2, s3, s4);
        #endif
    };
};

template<>
struct simd_trunc<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_round_ps(x.data, _MM_FROUND_TRUNC);
        #else
            float s1    = scalar_func::trunc(x.get<0>());
            float s2    = scalar_func::trunc(x.get<1>());
            float s3    = scalar_func::trunc(x.get<2>());
            float s4    = scalar_func::trunc(x.get<3>());

            return simd_type(s1, s2, s3, s4);
        #endif
    };
};

template<>
struct simd_eeq<float, 128, sse_tag> : simd_cmp_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpeq_ps(x.data, y.data);
    };
};

template<>
struct simd_neq<float, 128, sse_tag> : simd_cmp_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpneq_ps(x.data, y.data);
    };
};

template<>
struct simd_lt<float, 128, sse_tag> : simd_cmp_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmplt_ps(x.data, y.data);
    };
};

template<>
struct simd_gt<float, 128, sse_tag> : simd_cmp_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpgt_ps(x.data, y.data);
    };
};

template<>
struct simd_leq<float, 128, sse_tag> : simd_cmp_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmple_ps(x.data, y.data);
    };
};

template<>
struct simd_geq<float, 128, sse_tag> : simd_cmp_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpge_ps(x.data, y.data);
    };
};

template<>
struct simd_any_nan<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        __m128 nt   = _mm_cmp_ps(x.data, x.data, _CMP_NEQ_UQ);
        int res     = _mm_movemask_ps(nt);

        return res != 0;
    };
};

template<>
struct simd_any<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        int res     = _mm_movemask_ps(x.data);
        return res != 0;
    };
};

template<>
struct simd_all<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        int res     = _mm_movemask_ps(x.data);
        return res == 15;
    };
};

}}
