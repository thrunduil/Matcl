/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-simd/simd.h"

namespace matcl { namespace simd
{

//TODO

// check if there exists simd type for given value type of given length
template<class V, int Size>
struct has_simd_size
{ 
    static const bool value = false;
};

template<>
struct has_simd_size<double,2>
{ 
    static const bool value = true;
};

template<>
struct has_simd_size<double,4>
{ 
    static const bool value = true;
};

template<>
struct has_simd_size<float,4>
{ 
    static const bool value = true;
};

template<>
struct has_simd_size<float,8>
{ 
    static const bool value = true;
};

template<>
struct has_simd_size<simd_double_complex, 1>
{ 
    static const bool value = true;
};

template<>
struct has_simd_size<simd_double_complex, 2>
{ 
    static const bool value = true;
};

template<>
struct has_simd_size<simd_single_complex, 2>
{ 
    static const bool value = true;
};

template<>
struct has_simd_size<simd_single_complex, 4>
{ 
    static const bool value = true;
};

// simd type for storing values of given type with given size
template<class V, int Size>
struct default_simd_type_elems 
{};

template<>
struct default_simd_type_elems<double,2>
{ 
    using type = simd<double, 128, sse_tag>; 
};

template<>
struct default_simd_type_elems<double,4>
{ 
    using type = simd<double, 256, avx_tag>; 
};

template<>
struct default_simd_type_elems<float,4>
{ 
    using type = simd<float, 128, sse_tag>; 
};

template<>
struct default_simd_type_elems<float,8>
{ 
    using type = simd<float, 256, avx_tag>; 
};

template<>
struct default_simd_type_elems<simd_double_complex, 1>
{ 
    using type = simd_compl<double, 128, sse_tag>; 
};

template<>
struct default_simd_type_elems<simd_double_complex, 2>
{ 
    using type = simd_compl<double, 256, avx_tag>; 
};

template<>
struct default_simd_type_elems<simd_single_complex, 2>
{ 
    using type = simd_compl<float, 128, sse_tag>; 
};

template<>
struct default_simd_type_elems<simd_single_complex, 4>
{ 
    using type = simd_compl<float, 256, avx_tag>; 
};

}}

