/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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

namespace matcl { namespace simd
{

template<>
struct default_scalar_simd_type<double>
{ 
    using type = simd<double, 128, scalar_nosimd_tag>; 
};

template<>
struct default_scalar_simd_type<float>
{ 
    using type = simd<float, 128, scalar_nosimd_tag>; 
};

template<>
struct default_scalar_simd_type<int32_t>
{ 
    using type = simd<int32_t, 128, scalar_nosimd_tag>; 
};

template<>
struct default_scalar_simd_type<int64_t>
{ 
    using type = simd<int64_t, 128, scalar_nosimd_tag>; 
};

}}