/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017
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

}}}