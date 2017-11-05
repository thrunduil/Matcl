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
#include "matcl-simd/arch/avx/helpers.h"

namespace matcl { namespace simd
{

template<>
struct simd_reverse<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            static const __m256i control = details::vector_8_int<7,6,5,4,3,2,1,0>();
            return _mm256_permutevar8x32_ps(x.data, control);
        #else
            return simd_type(reverse(x.extract_high()), reverse(x.extract_low()));
        #endif
    };
};

template<>
struct simd_mult<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_mul_ps( x.data, y.data );
    };
};

template<>
struct simd_div<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_div_ps( x.data, y.data );
    };
};

template<>
struct simd_plus<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_add_ps( x.data, y.data );
    };
};

template<>
struct simd_minus<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_sub_ps( x.data, y.data );
    };
};

template<>
struct simd_uminus<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        float Z = -0.0f;
        return _mm256_xor_ps(x.data, _mm256_broadcast_ss(&Z));
    };
};

template<>
struct simd_sum_all<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static float eval(const simd_type& x)
    {
        float s = sum_all(x.extract_low()) + sum_all(x.extract_high());
        return s;
    };
};

template<>
struct simd_sub_add<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        __m256 s   = _mm256_addsub_ps (x.data, y.data);
        return s;
    };
};

template<>
struct simd_fma<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fmadd_ps( x.data, y.data, z.data);
        #else
            return x * y + z;
        #endif
    };
};

template<>
struct simd_fms<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fmsub_ps( x.data, y.data, z.data);
        #else
            return x * y - z;
        #endif
    };
};

template<>
struct simd_abs<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const __m256 sign_mask = _mm256_set1_ps(-0.f); // -0.f = 1 << 31
        return _mm256_andnot_ps(sign_mask, x.data);
    };
};

template<>
struct simd_max<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_max_ps(x.data, y.data);
    };
};

template<>
struct simd_min<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_min_ps(x.data, y.data);
    };
};

template<>
struct simd_sqrt<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_sqrt_ps(x.data);
    };
};

template<>
struct simd_round<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_ps(x.data, _MM_FROUND_NINT);
    };
};

template<>
struct simd_floor<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_ps(x.data, _MM_FROUND_FLOOR);
    };
};

template<>
struct simd_ceil<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_ps(x.data, _MM_FROUND_CEIL);
    };
};

template<>
struct simd_trunc<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_ps(x.data, _MM_FROUND_TRUNC);
    };
};

template<>
struct simd_eeq<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_ps(x.data, y.data, _CMP_EQ_OQ);
    };
};

template<>
struct simd_neq<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_ps(x.data, y.data, _CMP_NEQ_OQ);
    };
};

template<>
struct simd_lt<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_ps(x.data, y.data, _CMP_LT_OQ);
    };
};

template<>
struct simd_gt<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_ps(x.data, y.data, _CMP_GT_OQ);
    };
};

template<>
struct simd_leq<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_ps(x.data, y.data, _CMP_LE_OQ);
    };
};

template<>
struct simd_geq<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_ps(x.data, y.data, _CMP_GE_OQ);
    };
};

}}
