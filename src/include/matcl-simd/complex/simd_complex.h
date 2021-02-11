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

#include "matcl-simd/config.h"

#if MATCL_USE_MATCL_COMPLEX
    #include "matcl-core/matrix/scalar_types.h"
    #include "matcl-core/details/complex_details.h"
#endif

#include "matcl-simd/simd_general.h"

namespace matcl { namespace simd
{

#if MATCL_USE_MATCL_COMPLEX

    using simd_double_complex   = matcl::Complex;
    using simd_single_complex   = matcl::Float_complex;

#else

    //types simd_double_complex and simd_single_complex
    //must already be defined by a user in the namespace matcl::simd

#endif

}}

namespace matcl { namespace simd { namespace details
{

template<class T>   struct real_type                        {};
template<>          struct real_type<simd_single_complex>   { using type = float; };
template<>          struct real_type<simd_double_complex>   { using type = double; };

template<class T, int Bits, class Tag>
                    struct real_type<simd_compl<T,Bits,Tag>>{ using type = simd<T,Bits,Tag>; };

template<class T>   struct complex_type                     {};
template<>          struct complex_type<float>              { using type = simd_single_complex; };
template<>          struct complex_type<double>             { using type = simd_double_complex; };

template<class T, int Bits, class Tag>
                    struct complex_type<simd<T,Bits,Tag>>   { using type = simd_compl<T,Bits,Tag>; };

}}}
