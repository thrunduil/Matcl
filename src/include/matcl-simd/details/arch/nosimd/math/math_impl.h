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

template<int Bits>
struct simd_exp<double, Bits, nosimd_tag>
{
    using simd_type     = simd<double, Bits, nosimd_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        static const 
        int vec_size    = simd_type::vector_size;

        simd_type res;

        const double* ptr_a = a.get_raw_ptr();
        double* ptr_res     = res.get_raw_ptr();

        for (int i = 0; i < vec_size; ++i)
            ptr_res[i]      = ms::exp(ptr_a[i]);

        return res;
    };
};

template<int Bits>
struct simd_exp<float, Bits, nosimd_tag>
{
    using simd_type     = simd<float, Bits, nosimd_tag>;

    force_inline
    static simd_type eval(const simd_type& a)
    {
        static const 
        int vec_size    = simd_type::vector_size;

        simd_type res;

        const float* ptr_a  = a.get_raw_ptr();
        float* ptr_res      = res.get_raw_ptr();

        for (int i = 0; i < vec_size; ++i)
            ptr_res[i]      = ms::exp(ptr_a[i]);

        return res;
    };
};

}}}
