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
struct simd_reverse<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return reverse(x.reinterpret_as_float()).reinterpret_as_int32();
    };
};

template<>
struct simd_mult<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {        
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_mullo_epi32(x.data, y.data);
    
        #else
            return simd_type( x.extract_low() * y.extract_low(), 
                              x.extract_high() * y.extract_high());
        #endif
    };
};

template<>
struct simd_plus<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_add_epi32( x.data, y.data );
    };
};

template<>
struct simd_minus<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return _mm256_sub_epi32( x.data, y.data );
    };
};

template<>
struct simd_uminus<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return _mm256_sub_epi32(_mm256_setzero_si256(), x.data);
    };
};

template<>
struct simd_abs<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {        
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_sign_epi32(x.data, x.data);
        #else   
            return simd_type(abs(x.extract_low()), 
                             abs(x.extract_high()));       
        #endif
    };
};

template<>
struct simd_horizontal_sum<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static int32_t eval(const simd_type& x)
    {
        int32_t s = horizontal_sum(x.extract_low() + x.extract_high());
        return s;
    };
};

template<>
struct simd_horizontal_min<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static int32_t eval(const simd_type& x)
    {
        int32_t s = horizontal_min(min(x.extract_low(), x.extract_high()));
        return s;
    };
};

template<>
struct simd_horizontal_max<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static int32_t eval(const simd_type& x)
    {
        int32_t s = horizontal_max(max(x.extract_low(), x.extract_high()));
        return s;
    };
};

template<>
struct simd_bitwise_or<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {        
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_or_si256(x.data, y.data);
        #else   
            return bitwise_or(x.reinterpret_as_float(), y.reinterpret_as_float())
                        .reinterpret_as_int32();
        #endif  
    };
};

template<>
struct simd_bitwise_xor<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_xor_si256(x.data, y.data);
        #else   
            return bitwise_xor(x.reinterpret_as_float(), y.reinterpret_as_float())
                        .reinterpret_as_int32();
        #endif        
    };
};

template<>
struct simd_bitwise_and<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_and_si256(x.data, y.data);
        #else   
            return bitwise_and(x.reinterpret_as_float(), y.reinterpret_as_float())
                        .reinterpret_as_int32();
        #endif        
    };
};

template<>
struct simd_bitwise_andnot<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_andnot_si256(x.data, y.data);
        #else   
            return bitwise_andnot(x.reinterpret_as_float(), y.reinterpret_as_float())
                        .reinterpret_as_int32();
        #endif        
    };
};

template<>
struct simd_bitwise_not<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

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
struct simd_shift_left<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_slli_epi32(x.data, y);
        #else
            using simd_half = simd_type::simd_half;

            simd_half hi    = shift_left(x.extract_high(), y);
            simd_half lo    = shift_left(x.extract_low(), y);

            return simd_type(lo, hi);
        #endif
    };
};

template<>
struct simd_shift_right<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_srli_epi32(x.data, y);
        #else
            using simd_half = simd_type::simd_half;

            simd_half hi    = shift_right(x.extract_high(), y);
            simd_half lo    = shift_right(x.extract_low(), y);

            return simd_type(lo, hi);
        #endif
    };
};

template<>
struct simd_shift_right_arithmetic<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_srai_epi32(x.data, y);
        #else
            using simd_half = simd_type::simd_half;

            simd_half hi    = shift_right_arithmetic(x.extract_high(), y);
            simd_half lo    = shift_right_arithmetic(x.extract_low(), y);

            return simd_type(lo, hi);
        #endif
    };
};

template<>
struct simd_round<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_floor<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_ceil<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};  

template<>
struct simd_trunc<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    };
};

template<>
struct simd_eeq<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_cmpeq_epi32(x.data, y.data);
        #else
            return simd_type(eeq(x.extract_low(), y.extract_low()),
                             eeq(x.extract_high(), y.extract_high()));            
        #endif
    };
};

template<>
struct simd_gt<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_cmpgt_epi32(x.data, y.data);
        #else
            return simd_type(gt(x.extract_low(), y.extract_low()),
                             gt(x.extract_high(), y.extract_high()));            
        #endif        
    };
};

template<>
struct simd_neq<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return bitwise_not(eeq(x, y));
    };
};

template<>
struct simd_lt<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return gt(y, x);
    };
};

template<>
struct simd_leq<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return geq(y, x);
    };
};

template<>
struct simd_geq<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return bitwise_not(gt(y, x));
    };
};

template<>
struct simd_max<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_max_epi32(x.data, y.data);
        #elif MATCL_ARCHITECTURE_HAS_SSE41
            return simd_type(max(x.extract_low(), y.extract_low()),
                             max(x.extract_high(), y.extract_high()));
        #else   
            return if_then_else(gt(x, y), x, y);
        #endif
    };
};

template<>
struct simd_min<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        #if MATCL_ARCHITECTURE_HAS_AVX2
            return _mm256_min_epi32(x.data, y.data);
        #elif MATCL_ARCHITECTURE_HAS_SSE41
            return simd_type(min(x.extract_low(), y.extract_low()),
                             min(x.extract_high(), y.extract_high()));
        #else
            return if_then_else(gt(x, y), y, x);
        #endif  
    };
};

template<>
struct simd_any<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static bool eval(const simd_type& x)
    {
        return any(x.reinterpret_as_float());
    };
};

template<>
struct simd_all<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

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
struct simd_if_then_else<int32_t, 256, avx_tag>
{
    using simd_type = simd<int32_t, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        return if_then_else(test.reinterpret_as_float(), val_true.reinterpret_as_float(),
                            val_false.reinterpret_as_float()).reinterpret_as_int32();
    };
};

}}}
