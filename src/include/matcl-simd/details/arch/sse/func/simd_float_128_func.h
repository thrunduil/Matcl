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
#include "matcl-simd/details/scalar_func.h"
#include "matcl-core/details/float/fma_dekker_simd.inl"

namespace matcl { namespace simd { namespace details
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
        const simd_type mzero = simd_type::minus_zero();
        return _mm_xor_ps( x.data, mzero.data );
    };
};

template<>
struct simd_sum_all<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static float eval(const simd_type& x)
    {
        __m128 xp       = _mm_shuffle_ps(x.data, x.data, _MM_SHUFFLE(2, 3, 0, 1));
        __m128 sums     = _mm_add_ps(x.data, xp);
        xp              = _mm_movehl_ps(xp, sums);
        sums            = _mm_add_ss(sums, xp);

        return _mm_cvtss_f32(sums);
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
            const float* xp = x.get_raw_ptr();
            const float* yp = y.get_raw_ptr();

            float s1    = xp[0] - yp[0];
            float s2    = xp[1] + yp[1];
            float s3    = xp[2] - yp[2];
            float s4    = xp[3] + yp[3];
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

//
template<>
struct simd_fnma_f<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fnmadd_ps( x.data, y.data, z.data);
        #else
            return z - x * y;
        #endif
    };
};

template<>
struct simd_fnms_f<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fnmsub_ps( x.data, y.data, z.data);
        #else
            return -(x * y + z);
        #endif
    };
};

template<>
struct simd_fma_a<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fmadd_ps( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(x, y, z);
        #endif
    };
};

template<>
struct simd_fms_a<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fmsub_ps( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(x, y, -z);
        #endif
    };
};

//
template<>
struct simd_fnma_a<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fnmadd_ps( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(-x, y, z);
        #endif
    };
};

template<>
struct simd_fnms_a<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm_fnmsub_ps( x.data, y.data, z.data);
        #else
            return -fma_dekker_simd(x, y, z);
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
        const simd_type mzero = simd_type::minus_zero();
        return _mm_andnot_ps(mzero.data, x.data);
    };
};

template<>
struct simd_bitwise_or<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_or_ps(x.data, y.data);
    };
};

template<>
struct simd_bitwise_xor<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_xor_ps(x.data, y.data);
    };
};

template<>
struct simd_bitwise_and<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_and_ps(x.data, y.data);
    };
};

template<>
struct simd_bitwise_andnot<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_andnot_ps(x.data, y.data);
    };
};

template<>
struct simd_bitwise_not<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return bitwise_not(x.reinterpret_as_int32()).reinterpret_as_float();
    };
};

template<>
struct simd_shift_left<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        // cast to integer
        __m128i val_i   = _mm_castps_si128(x.data);

        // shift packed 32-bit integers in a left 
        __m128i res_i   = _mm_slli_epi32(val_i, y);

        //cast to float
        __m128 res_d    = _mm_castsi128_ps(res_i);
        return res_d;
    };
};

template<>
struct simd_shift_right<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        // cast to integer
        __m128i val_i   = _mm_castps_si128(x.data);

        // shift packed 32-bit integers in a left 
        __m128i res_i   = _mm_srli_epi32(val_i, y);

        //cast to float
        __m128 res_d    = _mm_castsi128_ps(res_i);
        return res_d;
    };
};

template<>
struct simd_shift_right_arithmetic<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return shift_right_arithmetic(x.reinterpret_as_int32(), y).reinterpret_as_float();
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
            const float* xp = x.get_raw_ptr();

            float s1    = scalar_func::round(xp[0]);
            float s2    = scalar_func::round(xp[1]);
            float s3    = scalar_func::round(xp[2]);
            float s4    = scalar_func::round(xp[3]);

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
            const float* xp = x.get_raw_ptr();

            float s1    = scalar_func::floor(xp[0]);
            float s2    = scalar_func::floor(xp[1]);
            float s3    = scalar_func::floor(xp[2]);
            float s4    = scalar_func::floor(xp[3]);

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
            const float* xp = x.get_raw_ptr();

            float s1    = scalar_func::ceil(xp[0]);
            float s2    = scalar_func::ceil(xp[1]);
            float s3    = scalar_func::ceil(xp[2]);
            float s4    = scalar_func::ceil(xp[3]);

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
            const float* xp = x.get_raw_ptr();

            float s1    = scalar_func::trunc(xp[0]);
            float s2    = scalar_func::trunc(xp[1]);
            float s3    = scalar_func::trunc(xp[2]);
            float s4    = scalar_func::trunc(xp[3]);

            return simd_type(s1, s2, s3, s4);
        #endif
    };
};

template<>
struct simd_eeq<float, 128, sse_tag> : simd_bool_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpeq_ps(x.data, y.data);
    };
};

template<>
struct simd_neq<float, 128, sse_tag> : simd_bool_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpneq_ps(x.data, y.data);
    };
};

template<>
struct simd_lt<float, 128, sse_tag> : simd_bool_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmplt_ps(x.data, y.data);
    };
};

template<>
struct simd_gt<float, 128, sse_tag> : simd_bool_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpgt_ps(x.data, y.data);
    };
};

template<>
struct simd_leq<float, 128, sse_tag> : simd_bool_base<float>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmple_ps(x.data, y.data);
    };
};

template<>
struct simd_geq<float, 128, sse_tag> : simd_bool_base<float>
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
struct simd_is_nan<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        __m128 nt   = _mm_cmp_ps(x.data, x.data, _CMP_NEQ_UQ);
        return nt;
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

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<>
struct simd_if_then_else<float, 128, sse_tag>
{
    using simd_type = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_blendv_ps(val_false.data, val_true.data, test.data);
        #else
            return _mm_or_ps( _mm_and_ps(test.data, val_true.data), 
                              _mm_andnot_ps(test.data, val_false.data));
        #endif
    };
};

}}}
