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
#include "matcl-simd/arch/sse/func/missing_intrinsics.h"

namespace matcl { namespace simd
{

template<>
struct simd_reverse<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

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

/*
TODO
template<>
struct simd_mult<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_mul_pd( x.data, y.data );
    };
};

template<>
struct simd_div<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_div_pd( x.data, y.data );
    };
};

template<>
struct simd_plus<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_add_pd( x.data, y.data );
    };
};

template<>
struct simd_minus<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_sub_pd( x.data, y.data );
    };
};

template<>
struct simd_uminus<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const simd_type mzero = simd_type::minus_zero();
        return _mm256_xor_pd(x.data, mzero.data);
    };
};

template<>
struct simd_sum_all<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static int64_t eval(const simd_type& x)
    {
        int64_t s    = sum_all(x.extract_low() + x.extract_high());
        return s;
    };
};

template<>
struct simd_sub_add<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        __m256d s   = _mm256_addsub_pd (x.data, y.data);
        return s;
    };
};

template<>
struct simd_fma_f<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

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
struct simd_fms_f<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

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
struct simd_fnma_f<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fnmadd_pd( x.data, y.data, z.data);
        #else
            return z - x * y;
        #endif
    };
};

template<>
struct simd_fnms_f<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fnmsub_pd( x.data, y.data, z.data);
        #else
            return -(x * y + z);
        #endif
    };
};

template<>
struct simd_fma_a<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fmadd_pd( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(x, y, z);
        #endif
    };
};

template<>
struct simd_fms_a<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fmsub_pd( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(x, y, -z);
        #endif
    };
};

template<>
struct simd_fnma_a<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fnmadd_pd( x.data, y.data, z.data);
        #else
            return fma_dekker_simd(-x, y, z);
        #endif
    };
};

template<>
struct simd_fnms_a<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        #if MATCL_ARCHITECTURE_HAS_FMA
            return _mm256_fnmsub_pd( x.data, y.data, z.data);
        #else
            return -fma_dekker_simd(x, y, z);
        #endif
    };
};

template<>
struct simd_abs<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const simd_type mzero = simd_type::minus_zero();
        return _mm256_andnot_pd(mzero.data, x.data);
    };
};

template<>
struct simd_bitwise_or<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_or_pd(x.data, y.data);
    };
};

template<>
struct simd_bitwise_xor<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_xor_pd(x.data, y.data);
    };
};

template<>
struct simd_bitwise_and<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_and_pd(x.data, y.data);
    };
};

template<>
struct simd_bitwise_andnot<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_andnot_pd(x.data, y.data);
    };
};

template<>
struct simd_shift_left<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            // cast to integer
            __m256i val_i   = _mm256_castpd_si256(x.data);

            // shift packed 64-bit integers in a left 
            __m256i res_i   = _mm256_slli_epi64(val_i, y);

            //cast to int64_t
            __m256d res_d   = _mm256_castsi256_pd(res_i);
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
struct simd_shift_right<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            // cast to integer
            __m256i val_i   = _mm256_castpd_si256(x.data);

            // shift packed 64-bit integers in a right 
            __m256i res_i   = _mm256_srli_epi64(val_i, y);

            //cast to int64_t
            __m256d res_d   = _mm256_castsi256_pd(res_i);
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
struct simd_max<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_max_pd(x.data, y.data);
    };
};

template<>
struct simd_min<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_min_pd(x.data, y.data);
    };
};

template<>
struct simd_sqrt<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_sqrt_pd(x.data);
    };
};

template<>
struct simd_round<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_pd(x.data, _MM_FROUND_NINT);
    };
};

template<>
struct simd_floor<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_pd(x.data, _MM_FROUND_FLOOR);
    };
};

template<>
struct simd_ceil<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_pd(x.data, _MM_FROUND_CEIL);
    };
};

template<>
struct simd_trunc<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_round_pd(x.data, _MM_FROUND_TRUNC);
    };
};

template<>
struct simd_eeq<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_EQ_OQ);
    };
};

template<>
struct simd_neq<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_NEQ_UQ);
    };
};

template<>
struct simd_lt<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_LT_OQ);
    };
};

template<>
struct simd_gt<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_GT_OQ);
    };
};

template<>
struct simd_leq<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_LE_OQ);
    };
};

template<>
struct simd_geq<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_cmp_pd(x.data, y.data, _CMP_GE_OQ);
    };
};

template<>
struct simd_any_nan<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        __m256d nt  = _mm256_cmp_pd(x.data, x.data, _CMP_NEQ_UQ);
        int res     = _mm256_movemask_pd(nt);

        return res != 0;
    };
};

template<>
struct simd_is_nan<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        __m256d nt  = _mm256_cmp_pd(x.data, x.data, _CMP_NEQ_UQ);
        return nt;
    };
};

template<>
struct simd_any<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        int res     = _mm256_movemask_pd(x.data);
        return res != 0;
    };
};

template<>
struct simd_all<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        int res     = _mm256_movemask_pd(x.data);
        return res == 15;
    };
};

//-----------------------------------------------------------------------
//                   MATHEMATICAL FUNCTIONS
//-----------------------------------------------------------------------
template<>
struct simd_pow2k<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;
    using simd_int  = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        // 2^52
        const int64_t pow2_52    = 4503599627370496.0;

        // bias in exponent
        const int64_t bias       = 1023.0;

        // put k + bias in least significant bits
        simd_type k2            = k + simd_type(bias + pow2_52);

        // shift left 52 places to get into exponent field
        simd_type pow2k         = shift_left(k2, 52);

        return pow2k;
    };
};

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<>
struct simd_if_then_else<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        return _mm256_blendv_pd(val_false.data, val_true.data, test.data);
    };
};

*/
}}
