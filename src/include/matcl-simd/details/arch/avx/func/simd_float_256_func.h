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
#include "matcl-simd/details/arch/avx/helpers.h"
#include "matcl-core/details/float/fma_dekker_simd.inl"

namespace matcl { namespace simd { namespace details
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
            __m256 xp   = _mm256_permute2f128_ps(x.data, x.data, 1);
            xp          = _mm256_permute_ps(xp, _MM_SHUFFLE(0,1,2,3));
            return xp;
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
        const simd_type mzero = simd_type::minus_zero();
        return _mm256_xor_ps(x.data, mzero.data);
    };
};

template<>
struct simd_sum_all<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static float eval(const simd_type& x)
    {
        float s = sum_all(x.extract_low() + x.extract_high());
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
struct simd_fma_f<float, 256, avx_tag>
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
struct simd_fms_f<float, 256, avx_tag>
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

//
template<>
struct simd_fnma_f<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fnmadd_ps( x.data, y.data, z.data);
        #else
            return z - x * y;
        #endif
    };
};

template<>
struct simd_fnms_f<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fnmsub_ps( x.data, y.data, z.data);
        #else
            return -(x * y + z);
        #endif
    };
};

template<>
struct simd_fma_a<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fmadd_ps( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(x, y, z);
        #endif
    };
};

template<>
struct simd_fms_a<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fmsub_ps( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(x, y, -z);
        #endif
    };
};

//
template<>
struct simd_fnma_a<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fnmadd_ps( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(-x, y, z);
        #endif
    };
};

template<>
struct simd_fnms_a<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fnmsub_ps( x.data, y.data, z.data);
        #else
            return -fma_dekker_simd(x, y, z);
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
        const simd_type mzero = simd_type::minus_zero();
        return _mm256_andnot_ps(mzero.data, x.data);
    };
};

template<>
struct simd_bitwise_or<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_or_ps(x.data, y.data);
    };
};

template<>
struct simd_bitwise_xor<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_xor_ps(x.data, y.data);
    };
};

template<>
struct simd_bitwise_and<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_and_ps(x.data, y.data);
    };
};

template<>
struct simd_bitwise_andnot<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_andnot_ps(x.data, y.data);
    };
};

template<>
struct simd_bitwise_not<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return bitwise_not(x.reinterpret_as_int32()).reinterpret_as_float();
    };
};

template<>
struct simd_shift_left<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            // cast to integer
            __m256i val_i   = _mm256_castps_si256(x.data);

            // shift packed 32-bit integers in a left 
            __m256i res_i   = _mm256_slli_epi32(val_i, y);

            //cast to float
            __m256 res_d    = _mm256_castsi256_ps(res_i);
            return res_d;
        #else
            using simd_half = simd_type::simd_half;

            simd_half hi    = shift_left(x.extract_high(), y);
            simd_half lo    = shift_left(x.extract_low(), y);

            return simd_type(lo, hi);
        #endif
    };
};

template<>
struct simd_shift_right<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            // cast to integer
            __m256i val_i   = _mm256_castps_si256(x.data);

            // shift packed 32-bit integers in a left 
            __m256i res_i   = _mm256_srli_epi32(val_i, y);

            //cast to float
            __m256 res_d    = _mm256_castsi256_ps(res_i);
            return res_d;
        #else
            using simd_half = simd_type::simd_half;

            simd_half hi    = shift_right(x.extract_high(), y);
            simd_half lo    = shift_right(x.extract_low(), y);

            return simd_type(lo, hi);
        #endif
    };
};

template<>
struct simd_shift_right_arithmetic<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return shift_right_arithmetic(x.reinterpret_as_int32(), y).reinterpret_as_float();
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
        return _mm256_cmp_ps(x.data, y.data, _CMP_NEQ_UQ);
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

template<>
struct simd_any_nan<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        __m256 nt   = _mm256_cmp_ps(x.data, x.data, _CMP_NEQ_UQ);
        int res     = _mm256_movemask_ps(nt);

        return res != 0;
    };
};

template<>
struct simd_is_nan<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        __m256 nt   = _mm256_cmp_ps(x.data, x.data, _CMP_NEQ_UQ);
        return nt;
    };
};

template<>
struct simd_any<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        int res     = _mm256_movemask_ps(x.data);
        return res != 0;
    };
};

template<>
struct simd_all<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        int res     = _mm256_movemask_ps(x.data);
        return res == 255;
    };
};

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<>
struct simd_if_then_else<float, 256, avx_tag>
{
    using simd_type = simd<float, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        return _mm256_blendv_ps(val_false.data, val_true.data, test.data);
    };
};

}}}
