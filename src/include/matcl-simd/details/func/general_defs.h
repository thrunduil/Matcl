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

#include "matcl-simd/simd_general.h"
#include "matcl-simd/details/func/simd_func_def.h"

namespace matcl { namespace simd { namespace details
{

template<int Bits, class Simd_tag>
struct simd_signbit_base<float, Bits, Simd_tag>
{
    using simd_type = simd<float, Bits, Simd_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return bitwise_and(x, simd_type::minus_zero());
    };
};

template<int Bits, class Simd_tag>
struct simd_signbit_base<double, Bits, Simd_tag>
{
    using simd_type = simd<double, Bits, Simd_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        return bitwise_and(x, simd_type::minus_zero());
    };
};

//-----------------------------------------------------------------------
//                  simd_reinterpret_as
//-----------------------------------------------------------------------
template<class Val, int Bits, class Tag>
struct simd_reinterpret_as<double, Val, Bits, Tag>
{
    using simd_in   = simd<Val, Bits, Tag>;
    using simd_ret  = simd<double, Bits, Tag>;

    force_inline
    static simd_ret eval(const simd_in& x)
    {
        return x.reinterpret_as_double();
    }
};

template<class Val, int Bits, class Tag>
struct simd_reinterpret_as<float, Val, Bits, Tag>
{
    using simd_in   = simd<Val, Bits, Tag>;
    using simd_ret  = simd<float, Bits, Tag>;

    force_inline
    static simd_ret eval(const simd_in& x)
    {
        return x.reinterpret_as_float();
    }
};

template<class Val, int Bits, class Tag>
struct simd_reinterpret_as<int32_t, Val, Bits, Tag>
{
    using simd_in   = simd<Val, Bits, Tag>;
    using simd_ret  = simd<int32_t, Bits, Tag>;

    force_inline
    static simd_ret eval(const simd_in& x)
    {
        return x.reinterpret_as_int32();
    }
};

template<class Val, int Bits, class Tag>
struct simd_reinterpret_as<int64_t, Val, Bits, Tag>
{
    using simd_in   = simd<Val, Bits, Tag>;
    using simd_ret  = simd<int64_t, Bits, Tag>;

    force_inline
    static simd_ret eval(const simd_in& x)
    {
        return x.reinterpret_as_int64();
    }
};

template<class Val, int Bits, class Tag>
struct simd_reinterpret_as<Val, Val, Bits, Tag>
{
    using simd_in   = simd<Val, Bits, Tag>;
    using simd_ret  = simd<Val, Bits, Tag>;

    force_inline
    static simd_ret eval(const simd_in& x)
    {
        return x;
    }
};

}}}
