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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/math/simd_math_func_def.h"

namespace matcl { namespace simd { namespace details
{

//-----------------------------------------------------------------------
//                          pow2k
//-----------------------------------------------------------------------
template<class Val, int Bits, class Simd_tag>
struct simd_pow2k_impl
{};

template<int Bits, class Simd_tag>
struct simd_pow2k_impl<double, Bits, Simd_tag>
{
    using simd_type = simd<double, Bits, Simd_tag>;
    using simd_int  = simd<int64_t, Bits, Simd_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        // 2^52
        const double pow2_52    = 4503599627370496.0;

        // bias in exponent
        const double bias       = 1023.0;

        // put k + bias in least significant bits
        simd_type k2            = k + simd_type(bias + pow2_52);

        // shift left 52 places to get into exponent field
        simd_type pow2k         = shift_left(k2, 52);

        return pow2k;
    };

    force_inline
    static simd_type eval_i(const simd_int& k)
    {
        simd_int ik         = k + simd_int(1023);
        ik                  = shift_left(ik, 52);
        simd_type res       = ik.reinterpret_as_double();
        return res;
    };
};

template<int Bits, class Simd_tag>
struct simd_pow2k_impl<float, Bits, Simd_tag>
{
    using simd_type = simd<float, Bits, Simd_tag>;
    using simd_int  = simd<int32_t, Bits, Simd_tag>;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        // 2^23
        const float pow2_23     = 8388608.0f;

        // bias in exponent
        const float bias        = 127.0f;

        // put k + bias in least significant bits
        simd_type k2            = k + simd_type(bias + pow2_23);

        // shift left 52 places to get into exponent field
        simd_type pow2k         = shift_left(k2, 23);

        return pow2k;
    };

    force_inline
    static simd_type eval_i(const simd_int& k)
    {
        simd_int ik         = k + simd_int(127);
        ik                  = shift_left(ik, 23);
        simd_type res       = ik.reinterpret_as_float();
        return res;
    };
};

//-----------------------------------------------------------------------
//                          fraction
//-----------------------------------------------------------------------
template<class Val, int Bits, class Simd_tag>
struct simd_fraction_impl
{};

template<int Bits, class Simd_tag>
struct simd_fraction_impl<double, Bits, Simd_tag>
{
    using simd_type = simd<double, Bits, Simd_tag>;
    using simd_int  = simd<int64_t, Bits, Simd_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_int mask_1     = simd_int(0x800FFFFFFFFFFFFFll);
        simd_int mask_2     = simd_int(0x3FE0000000000000ll);

        // set exponent to 0 + bias-1
        simd_type xs        = simd_type(x);
        simd_type x1        = ms::bitwise_and(mask_1.reinterpret_as_double(), xs);
        simd_type ret       = ms::bitwise_or(x1, mask_2.reinterpret_as_double());

        return ret;
    };
};

template<int Bits, class Simd_tag>
struct simd_fraction_impl<float, Bits, Simd_tag>
{
    using simd_type = simd<float, Bits, Simd_tag>;
    using simd_int  = simd<int32_t, Bits, Simd_tag>;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_int mask_1     = simd_int(0x807FFFFFl);
        simd_int mask_2     = simd_int(0x3F000000l);

        // set exponent to 0 + bias-1
        simd_type xs        = simd_type(x);
        simd_type x1        = ms::bitwise_and(mask_1.reinterpret_as_float(), xs);
        simd_type ret       = ms::bitwise_or(x1, mask_2.reinterpret_as_float());

        return ret;
    };
};

//-----------------------------------------------------------------------
//                          exponent
//-----------------------------------------------------------------------
template<class Val, int Bits, class Simd_tag>
struct simd_exponent_impl
{};

template<int Bits, class Simd_tag>
struct simd_exponent_impl<double, Bits, Simd_tag>
{
    using simd_type = simd<double, Bits, Simd_tag>;
    using simd_int  = simd<int64_t, Bits, Simd_tag>;

    force_inline
    static simd_int eval_i(const simd_type& x)
    {
        simd_int mask_1     = simd_int(0x7FF0000000000000ll);

        simd_int x1         = ms::bitwise_and(mask_1, x.reinterpret_as_int64());
        simd_int ret        = ms::shift_right(x1, 52) - simd_int(1022);

        return ret;
    };

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const double pow2_52    = 4503599627370496.0;   // 2^52
        const double bias       = 1022.0;               // bias - 1

        // shift exponent to low bits
        simd_type tmp1  = shift_left(x, 1);
        tmp1            = shift_right(tmp1, 53);

        // insert new exponent
        simd_type tmp2  = bitwise_or(tmp1, simd_type(pow2_52));

        // subtract magic number and bias
        simd_type res   = tmp2 - simd_type(pow2_52 + bias);
        return res;
    };
};

template<int Bits, class Simd_tag>
struct simd_exponent_impl<float, Bits, Simd_tag>
{
    using simd_type = simd<float, Bits, Simd_tag>;
    using simd_int  = simd<int32_t, Bits, Simd_tag>;

    force_inline
    static simd_int eval_i(const simd_type& x)
    {
        simd_int mask_1     = simd_int(0x7F800000l);

        simd_int x1         = ms::bitwise_and(mask_1, x.reinterpret_as_int32());
        simd_int ret        = ms::shift_right(x1, 23) - simd_int(126);

        return ret;    
    };

    force_inline
    static simd_type eval(const simd_type& x)
    {
        const float pow2_23    = 8388608.0f;    // 2^23
        const float bias       = 126.0f;        // bias - 1

        // shift exponent to low bits
        simd_type tmp1  = shift_left(x, 1);
        tmp1            = shift_right(tmp1, 24);

        // insert new exponent
        simd_type tmp2  = bitwise_or(tmp1, simd_type(pow2_23));

        // subtract magic number and bias
        simd_type res   = tmp2 - simd_type(pow2_23 + bias);
        return res;
    };
};

//-----------------------------------------------------------------------
//                          copysign
//-----------------------------------------------------------------------
template<class Val, int Bits, class Simd_tag>
struct simd_copysign_impl
{
    using simd_type = simd<Val, Bits, Simd_tag>;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return ms::bitwise_or(ms::signbit_base(y), ms::abs(x));
    };
};

}}}
