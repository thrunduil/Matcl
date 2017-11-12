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

#include "matcl-simd/config.h"
#include "matcl-core/details/mpl.h"

namespace matcl { namespace simd
{

// no simd instructions tag
struct nosimd_tag;

// sse instructions tag
struct sse_tag;

// avx instructions tag
struct avx_tag;

// simd type storing elements of type Val of total size Bits in bits
// using instruction set determined by Simd_tag
template<class Val, int Bits, class Simd_tag>
class simd;

// simd type storing elements of complex type with real type Val of total size
// Bits in bits, using instruction set determined by Simd_tag
template<class Val, int Bits, class Simd_tag>
class simd_compl;

// number of elements stored in given simd type
template<class Simd_type>
struct vector_size;

// simd type for storing values of given type with maximum size
template<class V>
struct default_simd_type;

// simd type for storing values of given type with given register size
// (allowed values for Bits are 128 and 256)
template<class V, int Bits>
struct default_simd_type_size;

}}

