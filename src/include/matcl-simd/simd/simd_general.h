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

namespace md = matcl::details;

// 128 bit vector tag
struct reg_128{};

// 256 bit vector tag
struct reg_256{};

// simd type; this class must be specialized for given value type 'Val'
// and simd tag 'Simd_tag'
template<class Val, class Simd_tag>
struct simd
{
    static_assert(md::dependent_false<Val>::value, "unsupported simd type");
};

// simd type; this class must be specialized for given value type 'Val' 
// and simd tag 'Simd_tag'
template<class Val, class Simd_tag>
struct simd_compl
{
    static_assert(md::dependent_false<Val>::value, "unsupported simd_compl type");
};

// number of elements stored in given simd type 'Simd_type'
template<class Simd_type>
struct vector_size
{
    static const int value = 1;
};

template<class Val, class Simd_tag>
struct vector_size<simd<Val, Simd_tag>>
{
    static const int value = simd<Val, Simd_tag>::vector_size;
};

template<class Val, class Simd_tag>
struct vector_size<simd_compl<Val, Simd_tag>>
{
    static const int value = simd_compl<Val, Simd_tag>::vector_size;
};

}}

