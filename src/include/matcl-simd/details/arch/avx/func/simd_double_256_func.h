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
#include "matcl-simd/details/func/simd_func_def.h"
#include "matcl-core/details/float/fma_dekker_simd.inl"

namespace matcl { namespace simd { namespace details
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
            // swap LO <-> HI
            __m256d xp  = _mm256_permute2f128_pd(x.data, x.data, 1);

            // swap within 128-bit lanes
            xp          = _mm256_permute_pd(xp, 5);
            return xp;
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
        const simd_type mzero = simd_type::minus_zero();
        return _mm256_xor_pd(x.data, mzero.data);
    };
};

template<>
struct simd_horizontal_sum<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static double eval(const simd_type& x)
    {
        double s    = horizontal_sum(x.extract_low() + x.extract_high());
        return s;
    };
};

template<>
struct simd_horizontal_min<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static double eval(const simd_type& x)
    {
        double s    = horizontal_min(min(x.extract_low(), x.extract_high()));
        return s;
    };
};

template<>
struct simd_horizontal_max<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static double eval(const simd_type& x)
    {
        double s    = horizontal_max(max(x.extract_low(), x.extract_high()));
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
struct simd_fma_f<double, 256, avx_tag>
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
struct simd_fms_f<double, 256, avx_tag>
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
struct simd_fnma_f<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

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
struct simd_fnms_f<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

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
struct simd_fma_a<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

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
struct simd_fms_a<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

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
struct simd_fnma_a<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

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
struct simd_fnms_a<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

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
struct simd_abs<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const simd_type mzero = simd_type::minus_zero();
        return _mm256_andnot_pd(mzero.data, x.data);
    };
};

template<>
struct simd_bitwise_or<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_or_pd(x.data, y.data);
    };
};

template<>
struct simd_bitwise_xor<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_xor_pd(x.data, y.data);
    };
};

template<>
struct simd_bitwise_and<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_and_pd(x.data, y.data);
    };
};

template<>
struct simd_bitwise_andnot<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_andnot_pd(x.data, y.data);
    };
};

template<>
struct simd_bitwise_not<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return bitwise_not(x.reinterpret_as_int64()).reinterpret_as_double();
    };
};

template<>
struct simd_shift_left<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            // cast to integer
            __m256i val_i   = _mm256_castpd_si256(x.data);

            // shift packed 64-bit integers in a left 
            __m256i res_i   = _mm256_slli_epi64(val_i, y);

            //cast to double
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
struct simd_shift_right<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            // cast to integer
            __m256i val_i   = _mm256_castpd_si256(x.data);

            // shift packed 64-bit integers in a right 
            __m256i res_i   = _mm256_srli_epi64(val_i, y);

            //cast to double
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
struct simd_shift_right_arithmetic<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return shift_right_arithmetic(x.reinterpret_as_int64(), y).reinterpret_as_double();
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
        return _mm256_cmp_pd(x.data, y.data, _CMP_NEQ_UQ);
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

template<>
struct simd_any_nan<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        __m256d nt  = _mm256_cmp_pd(x.data, x.data, _CMP_NEQ_UQ);
        int res     = _mm256_movemask_pd(nt);

        return res != 0;
    };
};

template<>
struct simd_is_nan<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        __m256d nt  = _mm256_cmp_pd(x.data, x.data, _CMP_NEQ_UQ);
        return nt;
    };
};

template<>
struct simd_is_finite<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        using simd_int  = simd<int64_t, 256, avx_tag>;

        simd_int xi     = x.reinterpret_as_int64();

        // mask selecting all bits in the exponent
        simd_int mask   = simd_int(0x7FF0000000000000ll);

        xi              = bitwise_and(xi, mask);

        // return true if at least one bit in the exponent is not set
        simd_int res    = neq(xi, mask);
        return res.reinterpret_as_double();
    };
};

template<>
struct simd_any<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        int res     = _mm256_movemask_pd(x.data);
        return res != 0;
    };
};

template<>
struct simd_all<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        int res     = _mm256_movemask_pd(x.data);
        return res == 15;
    };
};

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<>
struct simd_if_then_else<double, 256, avx_tag>
{
    using simd_type = simd<double, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        return _mm256_blendv_pd(val_false.data, val_true.data, test.data);
    };
};

}}}
