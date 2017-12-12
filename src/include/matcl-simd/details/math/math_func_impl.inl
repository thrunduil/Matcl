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

#include "matcl-simd/math_functions.h"

namespace matcl { namespace simd
{

namespace  ms = matcl::simd;

//-----------------------------------------------------------------------
//                   MATHEMATICAL FUNCTIONS
//-----------------------------------------------------------------------
template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::pow2k(const simd<Val, Bits, Simd_tag>& k)
{
    return details::simd_pow2k<Val, Bits, Simd_tag>::eval(k);
};

template<int Bits, class Simd_tag>
force_inline simd<float, Bits, Simd_tag> 
ms::pow2ki(const simd<int32_t, Bits, Simd_tag>& k)
{
    return details::simd_pow2k<float, Bits, Simd_tag>::eval_i(k);
};

template<int Bits, class Simd_tag>
force_inline simd<double, Bits, Simd_tag> 
ms::pow2ki(const simd<int64_t, Bits, Simd_tag>& k)
{
    return details::simd_pow2k<double, Bits, Simd_tag>::eval_i(k);
};

template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::exp(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_exp<Val, Bits, Simd_tag>::eval(x);
};

}}
