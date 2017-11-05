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

namespace matcl { namespace simd
{

template<>
struct simd_reverse<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_shuffle_pd(x.data, x.data, _MM_SHUFFLE2(0,1));
    };
};

template<>
struct simd_mult<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_mul_pd( x.data, y.data );
    };
};

template<>
struct simd_div<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_div_pd( x.data, y.data );
    };
};

template<>
struct simd_plus<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_add_pd( x.data, y.data );
    };
};

template<>
struct simd_minus<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_sub_pd( x.data, y.data );
    };
};

template<>
struct simd_uminus<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_xor_pd( x.data, _mm_set1_pd(-0.0) );
    };
};

template<>
struct simd_sum_all<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static double eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE3
            __m128d s   = _mm_hadd_pd(x.data,x.data);
            return s.m128d_f64[0];
        #else
            double s    = x.get<0>() + x.get<1>();
            return s;
        #endif
    };
};

template<>
struct simd_sub_add<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE3
            __m128d s   = _mm_addsub_pd(x.data, y.data);
            return s;
        #else
            double s1   = x.get<0>() + y.get<0>();
            double s2   = x.get<1>() - y.get<1>();
            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_fma<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fmadd_pd( x.data, y.data, z.data);
        #else
            return x * y + z;
        #endif
    };
};

template<>
struct simd_fms<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fmsub_pd( x.data, y.data, z.data);
        #else
            return x * y - z;
        #endif
    };
};

template<>
struct simd_abs<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const __m128d sign_mask = _mm_set1_pd(-0.); // -0. = 1 << 63
        return _mm_andnot_pd(sign_mask, x.data);
    };
};

template<>
struct simd_max<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_max_pd(x.data, y.data);
    };
};

template<>
struct simd_min<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_min_pd(x.data, y.data);
    };
};

template<>
struct simd_sqrt<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_sqrt_pd(x.data);
    };
};

template<>
struct simd_round<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_round_pd(x.data, _MM_FROUND_NINT);
        #else
            double s1   = scalar_func::round(x.get<0>());
            double s2   = scalar_func::round(x.get<1>());

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_floor<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_round_pd(x.data, _MM_FROUND_FLOOR);
        #else
            double s1   = scalar_func::floor(x.get<0>());
            double s2   = scalar_func::floor(x.get<1>());

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_ceil<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_round_pd(x.data, _MM_FROUND_CEIL);
        #else
            double s1   = scalar_func::ceil(x.get<0>());
            double s2   = scalar_func::ceil(x.get<1>());

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_trunc<double, 128, sse_tag>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_round_pd(x.data, _MM_FROUND_TRUNC);
        #else
            double s1   = scalar_func::trunc(x.get<0>());
            double s2   = scalar_func::trunc(x.get<1>());

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_eeq<double, 128, sse_tag> : simd_cmp_base<double>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX
            return _mm_cmp_pd(x.data, y.data, _CMP_EQ_OQ);
        #else
            double s1   = (x.get<0>() == y.get<0>()) ? val_true : val_false;
            double s2   = (x.get<1>() == y.get<1>()) ? val_true : val_false;

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_neq<double, 128, sse_tag> : simd_cmp_base<double>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX
            return _mm_cmp_pd(x.data, y.data, _CMP_NEQ_OQ);
        #else
            double s1   = (x.get<0>() != y.get<0>()) ? val_true : val_false;
            double s2   = (x.get<1>() != y.get<1>()) ? val_true : val_false;

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_lt<double, 128, sse_tag> : simd_cmp_base<double>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX
            return _mm_cmp_pd(x.data, y.data, _CMP_LT_OQ);
        #else
            double s1   = (x.get<0>() < y.get<0>()) ? val_true : val_false;
            double s2   = (x.get<1>() < y.get<1>()) ? val_true : val_false;

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_gt<double, 128, sse_tag> : simd_cmp_base<double>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX
            return _mm_cmp_pd(x.data, y.data, _CMP_GT_OQ);
        #else
            double s1   = (x.get<0>() > y.get<0>()) ? val_true : val_false;
            double s2   = (x.get<1>() > y.get<1>()) ? val_true : val_false;

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_leq<double, 128, sse_tag> : simd_cmp_base<double>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX
            return _mm_cmp_pd(x.data, y.data, _CMP_LE_OQ);
        #else
            double s1   = (x.get<0>() <= y.get<0>()) ? val_true : val_false;
            double s2   = (x.get<1>() <= y.get<1>()) ? val_true : val_false;

            return simd_type(s1, s2);
        #endif
    };
};

template<>
struct simd_geq<double, 128, sse_tag> : simd_cmp_base<double>
{
    using simd_type = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX
            return _mm_cmp_pd(x.data, y.data, _CMP_GE_OQ);
        #else
            double s1   = (x.get<0>() >= y.get<0>()) ? val_true : val_false;
            double s2   = (x.get<1>() >= y.get<1>()) ? val_true : val_false;

            return simd_type(s1, s2);
        #endif
    };
};

}}
