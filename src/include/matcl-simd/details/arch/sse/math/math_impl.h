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
#include "matcl-simd/details/math/simd_math_func_def.h"

namespace matcl { namespace simd { namespace details
{

template<>
struct simd_exp<double, 256, sse_tag>
{
    using simd_type     = simd<double, 256, sse_tag>;
    using simd_half     = simd<double, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        simd_half v1    = exp(a.extract_low());
        simd_half v2    = exp(a.extract_high());

        return simd_type(v1, v2);
    };
};

template<>
struct simd_exp<float, 256, sse_tag>
{
    using simd_type     = simd<float, 256, sse_tag>;
    using simd_half     = simd<float, 128, sse_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        simd_half v1    = exp(a.extract_low());
        simd_half v2    = exp(a.extract_high());

        return simd_type(v1, v2);
    };
};

}}}
