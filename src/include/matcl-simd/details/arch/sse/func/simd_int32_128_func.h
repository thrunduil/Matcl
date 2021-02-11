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
struct simd_reverse<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_shuffle_epi32(x.data, _MM_SHUFFLE(0,1,2,3));
    };
};

template<>
struct simd_mult<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_mullo_epi32(x.data, y.data);
        #else
            // (x1, x1, x3, x3)
            __m128i x13    = _mm_shuffle_epi32(x.data, 0xF5);

            // (y1, y1, y3, y3)
            __m128i y13    = _mm_shuffle_epi32(y.data, 0xF5);

            // (x0 * y0, xy0_hi, x2 * y2, xy2_hi)
            // signed integer multiplication with 2-complement representation
            // is equivalent to unsigned integer multiplication 
            __m128i prod02 = _mm_mul_epu32(x.data, y.data);

            // (x1 * y1, xy1_hi, x3 * y3, xy3_hi)
            __m128i prod13 = _mm_mul_epu32(x13, y13);

            // (x0 * y0, x1 * y1, xy0_hi, xy1_hi)
            __m128i prod01 = _mm_unpacklo_epi32(prod02, prod13);

            // (x2 * y2, x3 * y3, xy2_hi, xy3_hi)
            __m128i prod23 = _mm_unpackhi_epi32(prod02, prod13);

            // (x0 * y0, x1 * y1, x2 * y2, x3 * y3)
            return _mm_unpacklo_epi64(prod01, prod23);
        #endif
    };
};

template<>
struct simd_plus<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_add_epi32( x.data, y.data );
    };
};

template<>
struct simd_minus<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_sub_epi32( x.data, y.data );
    };
};

template<>
struct simd_uminus<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_sub_epi32(_mm_setzero_si128(), x.data);
    };
};

template<>
struct simd_abs<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE3
            return _mm_sign_epi32(x.data, x.data);
        #else
            // sign of x
            __m128i sign    = _mm_srai_epi32(x.data, 31);

            // invert bits if negative
            __m128i inv     = _mm_xor_si128(x.data, sign);

            // add 1
            return _mm_sub_epi32(inv,sign); 
        #endif
    };
};

template<>
struct simd_horizontal_sum<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static int32_t eval(const simd_type& x)
    {
        __m128i xp      = _mm_shuffle_epi32(x.data, _MM_SHUFFLE(2, 3, 0, 1));
        __m128i sums    = _mm_add_epi32(x.data, xp);
        xp              = missing::mm_movehl_epi32(xp, sums);
        sums            = _mm_add_epi32(sums, xp);

        return _mm_cvtsi128_si32(sums);
    };
};

template<>
struct simd_horizontal_min<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static int32_t eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            __m128i xp      = _mm_shuffle_epi32(x.data, _MM_SHUFFLE(2, 3, 0, 1));
            __m128i sums    = _mm_min_epi32(x.data, xp);
            xp              = missing::mm_movehl_epi32(xp, sums);
            sums            = _mm_min_epi32(sums, xp);

            return _mm_cvtsi128_si32(sums);
        #else
            __m128i xp      = _mm_shuffle_epi32(x.data, _MM_SHUFFLE(2, 3, 0, 1));
            __m128i sums    = missing::mm_min_epi32_sse(x.data, xp);
            xp              = missing::mm_movehl_epi32(xp, sums);
            sums            = missing::mm_min_epi32_sse(sums, xp);

            return _mm_cvtsi128_si32(sums);
        #endif
    };
};

template<>
struct simd_horizontal_max<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static int32_t eval(const simd_type& x)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            __m128i xp      = _mm_shuffle_epi32(x.data, _MM_SHUFFLE(2, 3, 0, 1));
            __m128i sums    = _mm_max_epi32(x.data, xp);
            xp              = missing::mm_movehl_epi32(xp, sums);
            sums            = _mm_max_epi32(sums, xp);

            return _mm_cvtsi128_si32(sums);
        #else
            __m128i xp      = _mm_shuffle_epi32(x.data, _MM_SHUFFLE(2, 3, 0, 1));
            __m128i sums    = missing::mm_max_epi32_sse(x.data, xp);
            xp              = missing::mm_movehl_epi32(xp, sums);
            sums            = missing::mm_max_epi32_sse(sums, xp);

            return _mm_cvtsi128_si32(sums);
        #endif
    };
};

//
template<>
struct simd_bitwise_or<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_or_si128(x.data, y.data);
    };
};

template<>
struct simd_bitwise_xor<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_xor_si128(x.data, y.data);
    };
};

template<>
struct simd_bitwise_and<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_and_si128(x.data, y.data);
    };
};

template<>
struct simd_bitwise_andnot<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_andnot_si128(x.data, y.data);
    };
};

template<>
struct simd_bitwise_not<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm_xor_si128(x.data, _mm_set1_epi32(-1));
    };
};

template<>
struct simd_shift_left<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return _mm_slli_epi32(x.data, y);
    };
};

template<>
struct simd_shift_right<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return _mm_srli_epi32(x.data, y);
    };
};

template<>
struct simd_shift_right_arithmetic<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return _mm_srai_epi32(x.data, y);
    };
};

template<>
struct simd_round<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_floor<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_ceil<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_trunc<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_eeq<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpeq_epi32(x.data, y.data);
    };
};

template<>
struct simd_gt<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm_cmpgt_epi32(x.data, y.data);
    };
};

template<>
struct simd_neq<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return bitwise_not(eeq(x, y));
    };
};

template<>
struct simd_lt<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return gt(y, x);
    };
};

template<>
struct simd_leq<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return geq(y, x);
    };
};

template<>
struct simd_geq<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return bitwise_not(gt(y, x));
    };
};

template<>
struct simd_max<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_max_epi32(x.data, y.data);
        #else
            return if_then_else(gt(x, y), x, y);
        #endif  
    };
};

template<>
struct simd_min<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_SSE41
            return _mm_min_epi32(x.data, y.data);
        #else
            return if_then_else(gt(x, y), y, x);
        #endif          
    };
};

template<>
struct simd_any<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        return any(x.reinterpret_as_float());
    };
};

template<>
struct simd_all<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        return all(x.reinterpret_as_float());
    };
};

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<>
struct simd_if_then_else<int32_t, 128, sse_tag>
{
    using simd_type = simd<int32_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        return if_then_else(test.reinterpret_as_float(), val_true.reinterpret_as_float(),
                            val_false.reinterpret_as_float()).reinterpret_as_int32();
    };
};

}}}
