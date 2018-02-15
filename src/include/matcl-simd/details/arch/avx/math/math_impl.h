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

//-----------------------------------------------------------------------
//                          pow2k
//-----------------------------------------------------------------------
template<class Val>
struct simd_pow2k<Val, 256, avx_tag>
{
    using simd_type = simd<Val, 256, avx_tag>;
    using int_type  = typename details::integer_type<Val>::type;
    using simd_int  = simd<int_type, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_pow2k_impl<Val, 256, avx_tag>::eval(k);
    };

    force_inline
    static simd_type eval_i(const simd_int& k)
    {
        return simd_pow2k_impl<Val, 256, avx_tag>::eval_i(k);
    };
};

//-----------------------------------------------------------------------
//                          exponent
//-----------------------------------------------------------------------
template<class Val>
struct simd_exponent<Val, 256, avx_tag>
{
    using simd_type = simd<Val, 256, avx_tag>;
    using int_type  = typename details::integer_type<Val>::type;
    using simd_int  = simd<int_type, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_exponent_impl<Val, 256, avx_tag>::eval(k);
    };

    force_inline
    static simd_int eval_i(const simd_type& k)
    {
        return simd_exponent_impl<Val, 256, avx_tag>::eval_i(k);
    };
};

//-----------------------------------------------------------------------
//                          copysign
//-----------------------------------------------------------------------
template<class Val>
struct simd_copysign<Val, 256, avx_tag>
{
    using simd_type = simd<Val, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_copysign_impl<Val, 256, avx_tag>::eval(x, y);
    };
};

//-----------------------------------------------------------------------
//                          fraction
//-----------------------------------------------------------------------
template<class Val>
struct simd_fraction<Val, 256, avx_tag>
{
    using simd_type = simd<Val, 256, avx_tag>;
    using int_type  = typename details::integer_type<Val>::type;
    using simd_int  = simd<int_type, 256, avx_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        return simd_fraction_impl<Val, 256, avx_tag>::eval(k);
    };
};

}}}
