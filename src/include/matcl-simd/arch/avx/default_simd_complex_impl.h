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

#include "matcl-simd/arch/simd_impl.h"
#include "matcl-simd/complex/simd_complex.h"

namespace matcl { namespace simd
{

template<>
struct default_simd_type<simd_double_complex>
{ 
    using type = simd_compl<double, 256, avx_tag>; 
};

template<>
struct default_simd_type<simd_single_complex>
{ 
    using type = simd_compl<float, 256, avx_tag>; 
};


//
template<>
struct default_simd_type_size<simd_double_complex, 256>
{
    using type = simd_compl<double, 256, avx_tag>; 
};

template<>
struct default_simd_type_size<simd_single_complex, 256>
{
    using type = simd_compl<float, 256, avx_tag>; 
};

}}