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

#if MATCL_ARCHITECTURE_HAS_SSE3
    #include "pmmintrin.h"
    #include "tmmintrin.h"
#endif

#if MATCL_ARCHITECTURE_HAS_SSE41
    #include "smmintrin.h"
#endif

#if MATCL_ARCHITECTURE_HAS_SSE42
    #include "nmmintrin.h"
#endif

namespace matcl { namespace simd { namespace missing
{

// move the upper double-precision (64-bit) floating-point element from b
// to the lower element of dst, and copy the upper element from a to the
// upper element of dst.
//
// Operation
//      dst[63:0]   := b[127:64]]
//      dst[127:64] := a[127:64]
force_inline
__m128d mm_movehl_pd(__m128d a, __m128d b)
{
    __m128 as   = _mm_castpd_ps(a);
    __m128 bs   = _mm_castpd_ps(b);
    __m128 rs   = _mm_movehl_ps(as, bs);
    __m128d r   = _mm_castps_pd(rs);

    return r;
};

// move the upper 64-bit integer element from b to the lower element of dst, 
// and copy the upper element from a to the upper element of dst.
//
// Operation
//      dst[63:0]   := b[127:64]]
//      dst[127:64] := a[127:64]
force_inline
__m128i mm_movehl_epi64(__m128i a, __m128i b)
{
    __m128 as   = _mm_castsi128_ps(a);
    __m128 bs   = _mm_castsi128_ps(b);
    __m128 rs   = _mm_movehl_ps(as, bs);
    __m128i r   = _mm_castps_si128(rs);

    return r;
};

// move the upper two 32-bit elements from b to the lower element of dst, 
// and copy the upper two 32-bits elements from a to the upper element of dst.
//
// Operation
//      dst[63:0]   := b[127:64]]
//      dst[127:64] := a[127:64]
force_inline
__m128i mm_movehl_epi32(__m128i a, __m128i b)
{
    __m128 as   = _mm_castsi128_ps(a);
    __m128 bs   = _mm_castsi128_ps(b);
    __m128 rs   = _mm_movehl_ps(as, bs);
    __m128i r   = _mm_castps_si128(rs);

    return r;
};

#if MATCL_ARCHITECTURE_HAS_SSE3
    // duplicate the low 64-bit integer element from a, and store the results in dst.
    //
    // Operation
    //      tmp[63:0] := a[63:0]
    //      tmp[127:64] := a[63:0]
    force_inline
    __m128i mm_movedup_epi64(__m128i a)
    {
        __m128d as  = _mm_castsi128_pd(a);
        __m128d r   = _mm_movedup_pd(as);

        return _mm_castpd_si128(r);
    };
#endif

// store the upper 64-bit integer element of a into memory.
//
// Operation
//      MEM[mem_addr+63:mem_addr] := a[127:64]
force_inline
void mm_storeh_epi64 (__m128i* mem_addr, __m128i a)
{
    _mm_storeh_pi((__m64 *)mem_addr, _mm_castsi128_ps(a));
}

// shuffle 64-bit integer elements using the control in imm8, and store the results in dst.
//
// Operation
//      dst[63:0]   := (imm8[0] == 0) ? a[63:0] : a[127:64]
//      dst[127:64] := (imm8[1] == 0) ? b[63:0] : b[127:64]
template<int imm8>
force_inline
__m128i mm_shuffle_epi64(__m128i a, __m128i b)
{
    __m128d ai  = _mm_castsi128_pd(a);
    __m128d bi  = _mm_castsi128_pd(b);
    __m128d res = _mm_shuffle_pd(ai, bi, imm8);

    return _mm_castpd_si128(res);
};

// compare packed 64-bit integers in a and b for equality, and store the results in dst
//
//  Operation:
//      FOR j := 0 to 1
//      	i := j*64
//      	dst[i+63:i] := ( a[i+63:i] == b[i+63:i] ) ? 0xFFFFFFFFFFFFFFFF : 0
//      ENDFOR
force_inline
__m128i mm_cmpeq_epi64_sse(__m128i x, __m128i y)
{
    // 32 bit compares
    __m128i cmp_32      = _mm_cmpeq_epi32(x, y);

    // swap low and high
    __m128i cmp_32_swap = _mm_shuffle_epi32(cmp_32, 0xB1);

    // both low and high must be -1
    __m128i test        = _mm_and_si128(cmp_32, cmp_32_swap);
    return test;
}

// compare packed 64-bit integers in x and y for greater-than, and store the results in dst
//
//  Operation:
//      FOR j := 0 to 1
//      	i := j*64
//      	dst[i+63:i] := ( x[i+63:i] > y[i+63:i] ) ? 0xFFFFFFFFFFFFFFFF : 0
//      ENDFOR
force_inline
__m128i mm_cmpgt_epi64_sse(__m128i x, __m128i y)
{
    // y - x
    __m128i s                   = _mm_sub_epi64(y, x);

    // y < x if y and x have same sign and s < 0 or (y < 0 and x >= 0)

    // check if x and y have different signs, i.e. x ^ y
    __m128i xy_sign_different   = _mm_xor_si128(x, y);

    // check if x >= 0 and y < 0, i.e. ~x & y
    __m128i x_pos_y_neg         = _mm_andnot_si128(x, y);

    // check s < 0 and sign(x) == sign(y), i.e. s & ~(x ^ y)
    __m128i cond_1              = _mm_andnot_si128(xy_sign_different,s);

    __m128i test                = _mm_or_si128(x_pos_y_neg, cond_1);

    // extend sign bit to 32 bits
    test                        = _mm_srai_epi32(test, 31);

    // extend sign bit to 64 bits
    test                        = _mm_shuffle_epi32(test, 0xF5);
    return  test;
};

// compare packed 32-bit integers in x and y, and store packed minimum values in dst.
force_inline
__m128i mm_min_epi32_sse(__m128i x, __m128i y)
{
    __m128i test        = _mm_cmpgt_epi32(x, y);

    // x > y ? y : x
    return _mm_or_si128( _mm_and_si128(test, y), _mm_andnot_si128(test, x));
}

// compare packed 32-bit integers in x and y, and store packed maximum values in dst
force_inline
__m128i mm_max_epi32_sse(__m128i x, __m128i y)
{
    __m128i test        = _mm_cmpgt_epi32(x, y);

    // x > y ? x : y
    return _mm_or_si128( _mm_and_si128(test, x), _mm_andnot_si128(test, y));
}


// compare packed 64-bit integer elements in x and y, and store packed minimum values in dst.
//
//  Operation
//      FOR j := 0 to 1
//      	i := j*64
//      	dst[i+63:i] := MIN(x[i+63:i], y[i+63:i])
//      ENDFOR
force_inline
__m128i mm_min_epi64(__m128i x, __m128i y)
{
    #if MATCL_ARCHITECTURE_HAS_SSE42
        __m128i test    = _mm_cmpgt_epi64(x, y);
    #else
        __m128i test    = missing::mm_cmpgt_epi64_sse(x, y);
    #endif

    #if MATCL_ARCHITECTURE_HAS_SSE41
        __m128d xd  = _mm_castsi128_pd(x);
        __m128d yd  = _mm_castsi128_pd(y);
        __m128d td  = _mm_castsi128_pd(test);
        __m128d res = _mm_blendv_pd(xd, yd, td);
        return _mm_castpd_si128(res);
    #else
        return _mm_or_si128( _mm_and_si128(test, y), 
                            _mm_andnot_si128(test, x));
    #endif
}

// compare packed 64-bit integer elements in x and y, and store packed maximum values in dst.
//
//  Operation
//      FOR j := 0 to 1
//      	i := j*64
//      	dst[i+63:i] := MAX(x[i+63:i], y[i+63:i])
//      ENDFOR
force_inline
__m128i mm_max_epi64(__m128i x, __m128i y)
{
    #if MATCL_ARCHITECTURE_HAS_SSE42
        __m128i test    = _mm_cmpgt_epi64(x, y);
    #else
        __m128i test    = missing::mm_cmpgt_epi64_sse(x, y);
    #endif

    #if MATCL_ARCHITECTURE_HAS_SSE41
        __m128d xd  = _mm_castsi128_pd(x);
        __m128d yd  = _mm_castsi128_pd(y);
        __m128d td  = _mm_castsi128_pd(test);
        __m128d res = _mm_blendv_pd(yd, xd, td);
        return _mm_castpd_si128(res);
    #else
        return _mm_or_si128( _mm_and_si128(test, x), 
                            _mm_andnot_si128(test, y));
    #endif
}

// copy 64-bit integer a to the lower element of dst, and zero the upper element.
//  Operation:
//      dst[63:0] := a[63:0]
//      dst[127:64] := 0
force_inline
__m128i mm_cvtsi64_si128(int64_t a)
{
    #if defined (_M_X64)
        return _mm_cvtsi64_si128(a);
    #else
        return _mm_set_epi64x(0, a);
    #endif
};

// copy the lower 64-bit integer in a to dst.
//  Operation:
//      dst[63:0] := a[63:0]
force_inline
int64_t mm_cvtsi128_si64(__m128i a)
{ 
    #if defined (_M_X64)
        return _mm_cvtsi128_si64(a);
    #else
        return reinterpret_cast<const int64_t*>(&a)[0];
    #endif
};

}}}
