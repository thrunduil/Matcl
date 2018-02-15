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

#include "matcl-simd/details/scalar_func.h"

namespace matcl { namespace simd { namespace scalar_func
{

//-------------------------------------------------------------------
//                         pow2k
//-------------------------------------------------------------------
struct pow2k_impl_double
{
    force_inline
    static double eval(double k)
    {
        // 2^52
        const double pow2_52    = 4503599627370496.0;

        // bias in exponent
        const double bias       = 1023.0;

        // put k + bias in least significant bits
        double k2               = k + (bias + pow2_52);

        // shift left 52 places to get into exponent field
        int64_t pow2k_i         = reinterpret_cast<int64_t&>(k2) << 52;

        return reinterpret_cast<double&>(pow2k_i);
    };
};

struct pow2k_impl_float
{
    force_inline
    static float eval(float k)
    {
        // 2^23
        const float pow2_23     = 8388608.0f;

        // bias in exponent
        const float bias        = 127.0f;

        // put k + bias in least significant bits
        float k2                = k + (bias + pow2_23);

        // shift left 52 places to get into exponent field
        int32_t pow2k_i         = reinterpret_cast<int32_t&>(k2) << 23;

        return reinterpret_cast<float&>(pow2k_i);
    };
};

struct pow2ki_impl_double
{
    force_inline
    static double eval(int64_t k)
    {
        int64_t ik          = k + int64_t(1023);
        ik                  = ik << 52;

        return reinterpret_cast<double&>(ik);
    };
};

struct pow2ki_impl_float
{
    force_inline
    static float eval(int32_t k)
    {
        int32_t ik          = k + int32_t(127);
        ik                  = ik << 23;

        return reinterpret_cast<float&>(ik);
    };
};

//-------------------------------------------------------------------
//                         fraction
//-------------------------------------------------------------------
struct fraction_impl_double
{
    force_inline
    static double eval(double x)
    {
        int64_t mask_1      = 0x800FFFFFFFFFFFFFll;
        int64_t mask_2      = 0x3FE0000000000000ll;

        // set exponent to 0 + bias-1
        int64_t xs          = reinterpret_cast<int64_t&>(x);
        int64_t x1          = scalar_func::bitwise_and(mask_1, xs);
        int64_t iret        = scalar_func::bitwise_or(x1, mask_2);

        double ret          = reinterpret_cast<double&>(iret);
        return ret;
    };
};

struct fraction_impl_float
{
    force_inline
    static float eval(float x)
    {
        int32_t mask_1      = 0x807FFFFFl;
        int32_t mask_2      = 0x3F000000l;

        // set exponent to 0 + bias-1
        int32_t xs          = reinterpret_cast<int32_t&>(x);
        int32_t x1          = scalar_func::bitwise_and(mask_1, xs);
        int32_t iret        = scalar_func::bitwise_or(x1, mask_2);

        float ret           = reinterpret_cast<float&>(iret);
        return ret;
    }
};

//-------------------------------------------------------------------
//                         exponent
//-------------------------------------------------------------------
struct exponent_impl_double
{
    force_inline
    static int64_t eval_i(double x)
    {
        int64_t mask_1      = 0x7FF0000000000000ll;

        int64_t xs          = reinterpret_cast<int64_t&>(x);
        int64_t x1          = scalar_func::bitwise_and(mask_1, xs);
        int64_t iret        = scalar_func::shift_right(x1, 52) - int64_t(1022);
        return iret;
    };

    force_inline
    static double eval(double x)
    {
        const double pow2_52    = 4503599627370496.0;   // 2^52
        const double bias       = 1022.0;               // bias - 1

        // shift exponent to low bits
        double tmp1     = scalar_func::shift_left(x, 1);
        tmp1            = scalar_func::shift_right(tmp1, 53);

        // insert new exponent
        double tmp2     = scalar_func::bitwise_or(tmp1, double(pow2_52));

        // subtract magic number and bias
        double res      = tmp2 - double(pow2_52 + bias);
        return res;
    };
};

struct exponent_impl_float
{
    force_inline
    static int32_t eval_i(float x)
    {
        int32_t mask_1      = 0x7F800000l;

        int32_t xs          = reinterpret_cast<int32_t&>(x);
        int32_t x1          = scalar_func::bitwise_and(mask_1, xs);
        int32_t iret        = scalar_func::shift_right(x1, 23) - int32_t(126);
        return iret;
    }

    force_inline
    static float eval(float x)
    {
        const float pow2_23    = 8388608.0f;    // 2^23
        const float bias       = 126.0f;        // bias - 1

        // shift exponent to low bits
        float tmp1      = scalar_func::shift_left(x, 1);
        tmp1            = scalar_func::shift_right(tmp1, 24);

        // insert new exponent
        float tmp2      = scalar_func::bitwise_or(tmp1, pow2_23);

        // subtract magic number and bias
        float res       = tmp2 - float(pow2_23 + bias);
        return res;
    }
};

//-------------------------------------------------------------------
//                         copysign
//-------------------------------------------------------------------
template<class Val>
struct copysign_impl
{
    force_inline
    static Val eval(Val x, Val y)
    {
        Val sign    = scalar_func::bitwise_and(y, Val(-0.0));
        return scalar_func::bitwise_or(sign, std::abs(x));
    };
};

//-------------------------------------------------------------------
//                         missing scalar functions
//-------------------------------------------------------------------

force_inline double pow2k(double f)     { return pow2k_impl_double::eval(f); }
force_inline float pow2k(float f)       { return pow2k_impl_float::eval(f); }

force_inline double pow2ki(int64_t f)   { return pow2ki_impl_double::eval(f); }
force_inline float pow2ki(int32_t f)    { return pow2ki_impl_float::eval(f); }

force_inline double fraction(double f)  { return fraction_impl_double::eval(f); }
force_inline float fraction(float f)    { return fraction_impl_float::eval(f); }

force_inline int64_t iexponent(double f){ return exponent_impl_double::eval_i(f); }
force_inline int32_t iexponent(float f) { return exponent_impl_float::eval_i(f); }

force_inline double exponent(double f)  { return exponent_impl_double::eval(f); }
force_inline float exponent(float f)    { return exponent_impl_float::eval(f); }

force_inline double copysign(double x, double y)    { return copysign_impl<double>::eval(x, y); }
force_inline float copysign(float x, float y)       { return copysign_impl<float>::eval(x, y); }

}}}
