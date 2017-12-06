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

//TODO
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
struct simd_shift_left<int64_t, 128, sse_tag>
{
    using simd_type = simd<int64_t, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return _mm_slli_epi64(x.data, y);
    };
};

}}
