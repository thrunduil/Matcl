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
#include "matcl-simd/details/arch/sse/func/missing_intrinsics.h"

namespace matcl { namespace simd { namespace details
{

template<>
struct simd_reverse<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        static const int control = 2 + (3 << 2) + (0 << 4) + (1 << 6);
        return _mm_shuffle_epi32(x.data, control);
    };
};

template<>
struct simd_mult<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            // split into 32-bit multiplies

            __m128i zero    = _mm_setzero_si128();

            // swap Hi and Lo: (y0_H, y0_L, y1_H, y1_L)
            __m128i y_swap  = _mm_shuffle_epi32(y.data, 0xB1);

            // 32 bit L*H products: (x0_L * y0_H, x0_H * y0_L, x1_L * y1_H, x1_H * y1_L)
            __m128i prod_lh = _mm_mullo_epi32(x.data, y_swap);
                        
            // (x0_L * y0_H + x0_H * y0_L, x1_L * y1_H + x1_H * y1_L, 0, 0)
            __m128i prod_lh2 = _mm_hadd_epi32(prod_lh, zero);

            // (0, x0_L * y0_H + x0_H * y0_L, 0, x1_L * y1_H + x1_H * y1_L)
            __m128i prod_lh3 = _mm_shuffle_epi32(prod_lh2, 0x73);

            // 64 bit unsigned products: (x0_L * y0_L, x1L * y1_L)
            __m128i prod_ll = _mm_mul_epu32(x.data, y.data);

            // (x0_L * y0_L + (x0_L * y0_H + x0_H * y0_L) << 32, 
            //  x1_L * y1_L + (x1_L * y1_H + x1_H * y1_L) << 32)
            __m128i prod    = _mm_add_epi64(prod_ll, prod_lh3);
            return  prod;
        #else
            const int64_t* ptr_x    = x.get_raw_ptr();
            const int64_t* ptr_y    = y.get_raw_ptr();

            int64_t xy0             = ptr_x[0] * ptr_y[0];
            int64_t xy1             = ptr_x[1] * ptr_y[1];

            return simd_type(xy0, xy1);
        #endif
    };
};

template<>
struct simd_plus<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_add_epi64( x.data, y.data );
    };
};

template<>
struct simd_minus<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_sub_epi64( x.data, y.data );
    };
};

template<>
struct simd_uminus<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_sub_epi64(_mm_setzero_si128(), x.data);
    };
};

template<>
struct simd_abs<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE42
            // 0 > x
            __m128i sign  = _mm_cmpgt_epi64(_mm_setzero_si128(), x.data);
        #else
            // sign in high int64
            __m128i signh = _mm_srai_epi32(x.data, 31);
            
            // copy sign to low int64
            __m128i sign  = _mm_shuffle_epi32(signh, 0xF5);
        #endif

        // invert bits if negative
        __m128i inv   = _mm_xor_si128(x.data, sign);

        // add 1
        return _mm_sub_epi64(inv,sign);
    };
};

template<>
struct simd_horizontal_sum<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static int64_t eval(const simd_type& x)
    {
        __m128i shuf    = missing::mm_movehl_epi64(x.data, x.data);
        __m128i res     = _mm_add_epi64(x.data, shuf);

        return missing::mm_cvtsi128_si64(res);
    };
};

template<>
struct simd_horizontal_min<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static int64_t eval(const simd_type& x)
    {
        __m128i shuf    = missing::mm_movehl_epi64(x.data, x.data);
        __m128i res     = missing::mm_min_epi64(x.data, shuf);

        return missing::mm_cvtsi128_si64(res);
    };
};

template<>
struct simd_horizontal_max<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static int64_t eval(const simd_type& x)
    {
        __m128i shuf    = missing::mm_movehl_epi64(x.data, x.data);
        __m128i res     = missing::mm_max_epi64(x.data, shuf);

        return missing::mm_cvtsi128_si64(res);
    };
};

template<>
struct simd_bitwise_or<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_or_si128(x.data, y.data);
    };
};

template<>
struct simd_bitwise_xor<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_xor_si128(x.data, y.data);
    };
};

template<>
struct simd_bitwise_and<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_and_si128(x.data, y.data);
    };
};

template<>
struct simd_bitwise_andnot<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_andnot_si128(x.data, y.data);
    };
};

template<>
struct simd_bitwise_not<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_xor_si128(x.data, _mm_set1_epi32(-1));
    };
};

template<>
struct simd_shift_left<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return _mm_slli_epi64(x.data, y);
    };
};

template<>
struct simd_shift_right<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return _mm_srli_epi64(x.data, y);
    };
};

template<>
struct simd_shift_right_arithmetic<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        using simd_32 = simd<int32_t, 128, sse_tag>;

        // instruction does not exist; split into 32-bit shifts
        if (y <= 32) 
        {
            // x >> y signed 32-bit ints
            __m128i sra  = _mm_srai_epi32(x.data, y);

            // x >> y unsigned 64-bit ints
            __m128i srl  = _mm_srli_epi64(x.data, y);

            // mask for signed high part
            __m128i mask = _mm_setr_epi32(0,-1,0,-1);

            return if_then_else(simd_32(mask), simd_32(sra), simd_32(srl))
                            .reinterpret_as_int64();
        }
        else 
        {  
            // y > 32
            y           = y - 32;

            // sign of x
            __m128i sign = _mm_srai_epi32(x.data, 31);

            // x >> (y-32) signed int32
            __m128i sra2 = _mm_srai_epi32(x.data, y);

            // x >> (y-32) >> 32 (second shift unsigned int64)
            __m128i sra3 = _mm_srli_epi64(sra2, 32);

            // mask for high part containing only sign
            __m128i mask = _mm_setr_epi32(0,-1,0,-1);
            return  if_then_else(simd_32(mask), simd_32(sign), simd_32(sra3))
                            .reinterpret_as_int64();
        };
    };
};

template<>
struct simd_round<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_floor<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_ceil<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_trunc<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_eeq<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_cmpeq_epi64(x.data, y.data);
        #else
            return missing::mm_cmpeq_epi64_sse(x.data, y.data);
        #endif
    };
};

template<>
struct simd_gt<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE42
            return _mm_cmpgt_epi64(x.data, y.data);
        #else
            return missing::mm_cmpgt_epi64_sse(x.data, y.data);
        #endif
    };
};

template<>
struct simd_neq<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return bitwise_not(eeq(x, y));
    };
};

template<>
struct simd_lt<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return gt(y, x);
    };
};

template<>
struct simd_leq<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return geq(y, x);
    };
};

template<>
struct simd_geq<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return bitwise_not(gt(y, x));
    };
};

template<>
struct simd_max<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return missing::mm_max_epi64(x.data, y.data);        
    };
};

template<>
struct simd_min<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return missing::mm_min_epi64(x.data, y.data);
    };
};

template<>
struct simd_any<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        return any(x.reinterpret_as_double());
    };
};

template<>
struct simd_all<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

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
struct simd_if_then_else<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        return if_then_else(test.reinterpret_as_double(), val_true.reinterpret_as_double(),
                            val_false.reinterpret_as_double()).reinterpret_as_int64();
    };
};

}}}
