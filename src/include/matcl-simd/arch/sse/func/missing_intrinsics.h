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

}}}
