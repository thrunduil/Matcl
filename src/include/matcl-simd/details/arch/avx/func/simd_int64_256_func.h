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
#include "matcl-simd/details/arch/sse/func/missing_intrinsics.h"

namespace matcl { namespace simd { namespace details
{

template<>
struct simd_reverse<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return reverse(x.reinterpret_as_double()).reinterpret_as_int64();
    };
};

template<>
struct simd_mult<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            // split into 32-bit multiplies

            __m256i zero    = _mm256_setzero_si256();

            // swap Hi and Lo: (y0_H, y0_L, y1_H, y1_L)
            __m256i y_swap  = _mm256_shuffle_epi32(y.data, 0xB1);

            // 32 bit L*H products: (x0_L * y0_H, x0_H * y0_L, x1_L * y1_H, x1_H * y1_L)
            __m256i prod_lh = _mm256_mullo_epi32(x.data, y_swap);
                        
            // (x0_L * y0_H + x0_H * y0_L, x1_L * y1_H + x1_H * y1_L, 0, 0)
            __m256i prod_lh2 = _mm256_hadd_epi32(prod_lh, zero);

            // (0, x0_L * y0_H + x0_H * y0_L, 0, x1_L * y1_H + x1_H * y1_L)
            __m256i prod_lh3 = _mm256_shuffle_epi32(prod_lh2, 0x73);

            // 64 bit unsigned products: (x0_L * y0_L, x1L * y1_L)
            __m256i prod_ll = _mm256_mul_epu32(x.data, y.data);

            // (x0_L * y0_L + (x0_L * y0_H + x0_H * y0_L) << 32, 
            //  x1_L * y1_L + (x1_L * y1_H + x1_H * y1_L) << 32)
            __m256i prod    = _mm256_add_epi64(prod_ll, prod_lh3);
            return  prod;
        #else
            simd_type res;
            int64_t* res_ptr        = res.get_raw_ptr();
            const int64_t* ptr_x    = x.get_raw_ptr();
            const int64_t* ptr_y    = y.get_raw_ptr();

            static const int vector_size    = simd_type::vector_size;

            for (int i = 0; i < vector_size; ++i)
                res_ptr[i]  = ptr_x[i] * ptr_y[i];

            return res;
        #endif
    };
};

template<>
struct simd_plus<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_add_epi64( x.data, y.data );
    };
};

template<>
struct simd_minus<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_sub_epi64( x.data, y.data );
    };
};

template<>
struct simd_uminus<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_sub_epi64(_mm256_setzero_si256(), x.data);
    };
};

template<>
struct simd_abs<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {        
        #if MATCL_ARCHITECTURE_HAS_AVX2
            __m256i sign  = _mm256_cmpgt_epi64(_mm256_setzero_si256(), x.data);

            // invert bits if negative
            __m256i inv   = _mm256_xor_si256(x.data, sign);

            // add 1
            return _mm256_sub_epi64(inv, sign);
        #else
            return simd_type(abs(x.extract_low()), 
                             abs(x.extract_high()));
        #endif
    };
};

template<>
struct simd_horizontal_sum<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static int64_t eval(const simd_type& x)
    {
        int64_t s    = horizontal_sum(x.extract_low() + x.extract_high());
        return s;
    };
};

template<>
struct simd_horizontal_min<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static int64_t eval(const simd_type& x)
    {
        int64_t s    = horizontal_min(min(x.extract_low(), x.extract_high()));
        return s;
    };
};

template<>
struct simd_horizontal_max<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static int64_t eval(const simd_type& x)
    {
        int64_t s    = horizontal_max(max(x.extract_low(), x.extract_high()));
        return s;
    };
};

template<>
struct simd_bitwise_or<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_or_si256(x.data, y.data);
        #else   
            return bitwise_or(x.reinterpret_as_double(), y.reinterpret_as_double())
                        .reinterpret_as_int64();
        #endif
    };
};

template<>
struct simd_bitwise_xor<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_xor_si256(x.data, y.data);
        #else   
            return bitwise_xor(x.reinterpret_as_double(), y.reinterpret_as_double())
                        .reinterpret_as_int64();
        #endif
    };
};

template<>
struct simd_bitwise_and<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_and_si256(x.data, y.data);
        #else   
            return bitwise_and(x.reinterpret_as_double(), y.reinterpret_as_double())
                        .reinterpret_as_int64();
        #endif        
    };
};

template<>
struct simd_bitwise_andnot<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_andnot_si256(x.data, y.data);
        #else   
            return bitwise_andnot(x.reinterpret_as_double(), y.reinterpret_as_double())
                        .reinterpret_as_int64();
        #endif        
    };
};

template<>
struct simd_bitwise_not<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_xor_si256(x.data, _mm256_set1_epi32(-1));
        #else   
            return simd_type(bitwise_not(x.extract_low()), bitwise_not(x.extract_high()));
        #endif        
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
            return _mm256_slli_epi64(x.data, y);
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
            return _mm256_srli_epi64(x.data, y);
        #else
            using simd_half = simd_type::simd_half;

            simd_half hi    = shift_right(x.extract_high(), y);
            simd_half lo    = shift_right(x.extract_low(), y);

            return simd_type(lo, hi);
        #endif
    };
};

template<>
struct simd_shift_right_arithmetic<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {        
        using simd_32 = simd<int32_t, 256, avx_tag>;

        #if MATCL_ARCHITECTURE_HAS_AVX2
            // instruction does not exist; split into 32-bit shifts
            if (y <= 32) 
            {
                // x >> y signed 32-bit ints
                __m256i sra  = _mm256_srai_epi32(x.data, y);

                // x >> y unsigned 64-bit ints
                __m256i srl  = _mm256_srli_epi64(x.data, y);

                // mask for signed high part
                __m256i mask = _mm256_setr_epi32(0,-1,0,-1, 0,-1,0,-1);

                return if_then_else(simd_32(mask), simd_32(sra), simd_32(srl))
                        .reinterpret_as_int64();
            }
            else 
            {  
                // y > 32
                y           = y - 32;

                // sign of x
                __m256i sign = _mm256_srai_epi32(x.data, 31);

                // x >> (y-32) signed int32
                __m256i sra2 = _mm256_srai_epi32(x.data, y);

                // x >> (y-32) >> 32 (second shift unsigned int64)
                __m256i sra3 = _mm256_srli_epi64(sra2, 32);

                // mask for high part containing only sign
                __m256i mask = _mm256_setr_epi32(0,-1,0,-1, 0,-1,0,-1);

                return  if_then_else(simd_32(mask), simd_32(sign), simd_32(sra3))
                                .reinterpret_as_int64();
            }
        #else   
            return simd_type(shift_right_arithmetic(x.extract_low(), y), 
                             shift_right_arithmetic(x.extract_high(), y));
        #endif
    };
};

template<>
struct simd_round<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_floor<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_ceil<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_trunc<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_eeq<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_cmpeq_epi64(x.data, y.data);
        #else
            return simd_type(eeq(x.extract_low(), y.extract_low()),
                             eeq(x.extract_high(), y.extract_high()));            
        #endif        
    };
};

template<>
struct simd_gt<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_cmpgt_epi64(x.data, y.data);
        #else
            return simd_type(gt(x.extract_low(), y.extract_low()),
                             gt(x.extract_high(), y.extract_high()));            
        #endif                
    };
};

template<>
struct simd_neq<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return bitwise_not(eeq(x, y));
    };
};

template<>
struct simd_lt<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return gt(y, x);
    };
};

template<>
struct simd_leq<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return geq(y, x);
    };
};

template<>
struct simd_geq<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return bitwise_not(gt(y, x));
    };
};

template<>
struct simd_max<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return if_then_else(gt(x, y), x, y);
    };
};

template<>
struct simd_min<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return if_then_else(gt(x, y), y, x);
    };
};

template<>
struct simd_any<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        return any(x.reinterpret_as_double());
    };
};

template<>
struct simd_all<int64_t, 256, avx_tag>
{
    using simd_type = simd<int64_t, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        return all(x.reinterpret_as_double());
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
        return if_then_else(test.reinterpret_as_double(), val_true.reinterpret_as_double(),
                            val_false.reinterpret_as_double()).reinterpret_as_int64();
    };
};

}}}
