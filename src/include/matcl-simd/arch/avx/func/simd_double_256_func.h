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
struct simd_reverse<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            static const int control = 3 + (2 << 2) + (1 << 4);
            return _mm256_permute4x64_pd(x.data, control);
        #else
            return simd_type(reverse(x.extract_high()), reverse(x.extract_low()));
        #endif
    };
};

template<>
struct simd_mult<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_mul_pd( x.data, y.data );
    };
};

template<>
struct simd_div<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_div_pd( x.data, y.data );
    };
};

template<>
struct simd_plus<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_add_pd( x.data, y.data );
    };
};

template<>
struct simd_minus<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_sub_pd( x.data, y.data );
    };
};

template<>
struct simd_uminus<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        double Z = -0.0;
        return _mm256_xor_pd(x.data, _mm256_broadcast_sd(&Z));
    };
};

template<>
struct simd_sum_all<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static double eval(const simd_type& x)
    {
        double s    = sum_all(x.extract_low()) + sum_all(x.extract_high());
        return s;
    };
};

template<>
struct simd_sub_add<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        __m256d s   = _mm256_addsub_pd (x.data, y.data);
        return s;
    };
};

template<>
struct simd_fma<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fmadd_pd( x.data, y.data, z.data);
        #else
            return x * y + z;
        #endif
    };
};

template<>
struct simd_fms<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fmsub_pd( x.data, y.data, z.data);
        #else
            return x * y - z;
        #endif
    };
};

template<>
struct simd_abs<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const __m256d sign_mask = _mm256_set1_pd(-0.); // -0. = 1 << 63
        return _mm256_andnot_pd(sign_mask, x.data);
    };
};

template<>
struct simd_max<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_max_pd(x.data, y.data);
    };
};

template<>
struct simd_min<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_min_pd(x.data, y.data);
    };
};

template<>
struct simd_sqrt<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_sqrt_pd(x.data);
    };
};

template<>
struct simd_round<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_pd(x.data, _MM_FROUND_NINT);
    };
};

template<>
struct simd_floor<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_pd(x.data, _MM_FROUND_FLOOR);
    };
};

template<>
struct simd_ceil<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_pd(x.data, _MM_FROUND_CEIL);
    };
};

template<>
struct simd_trunc<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_pd(x.data, _MM_FROUND_TRUNC);
    };
};

template<>
struct simd_eeq<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_EQ_OQ);
    };
};

template<>
struct simd_neq<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_NEQ_OQ);
    };
};

template<>
struct simd_lt<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_LT_OQ);
    };
};

template<>
struct simd_gt<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_GT_OQ);
    };
};

template<>
struct simd_leq<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_LE_OQ);
    };
};

template<>
struct simd_geq<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_GE_OQ);
    };
};

}}
