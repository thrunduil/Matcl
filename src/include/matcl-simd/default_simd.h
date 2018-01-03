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

namespace matcl { namespace simd
{

// simd type for storing values of given type with maximum size
template<class V>
struct default_simd_type 
{
    static_assert(md::dependent_false<V>::value, "unsupported simd type");
};

// simd type for storing one value of given type
template<class V>
struct default_scalar_simd_type
{
    static_assert(md::dependent_false<V>::value, "unsupported simd type");
};

// simd type for storing values of given type with given register size
// (allowed values for Bits are 128 and 256)
template<class V, int Bits>
struct default_simd_bit_size
{
    static_assert(md::dependent_false<V>::value, "unsupported simd type");
};

// simd type for storing N values of given type
// (allowed values for N are 1, 2, 4, 8)
template<class V, int N>
struct default_simd_vector_size
{
    using type = typename default_simd_bit_size<V, N * sizeof(V) * 8>::type;
};

template<class V>
struct default_simd_vector_size<V, 1>
{
    using type = typename default_scalar_simd_type<V>::type;
};

}}

// specialize default_simd_type type for given type and given architecture
#include "matcl-simd/details/arch/default_simd.h"
#include "matcl-simd/details/arch/default_simd_scalar.h"

