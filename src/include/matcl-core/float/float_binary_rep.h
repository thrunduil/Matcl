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

#include "matcl-core/config.h"
#include <cstdint>

namespace matcl 
{

// return single precision floating point number with binary representation
// given by 'bits' argument
float hex_float(const uint32_t& bits);

// return double precision floating point number with binary representation
// given by 'bits' argument
double hex_double(const uint64_t& bits);

// binary representation of IEEE 758 double precision number
struct double_binary_rep
{
    using base_type = uint64_t;
       
    static const int bias   = 1023;

    base_type    fraction   : 52;
    base_type    exponent   : 11;
    base_type    sign       : 1; 
};

// binary representation of IEEE 758 single precision number
struct float_binary_rep
{
    using base_type = uint32_t;
       
    static const int bias   = 127;

    base_type    fraction   : 23;
    base_type    exponent   : 8;
    base_type    sign       : 1; 
};

// decode pointers to IEEE 754 double value; denormal numbers are not
// supported;
union double_decoder
{
    public:
        // stored double precision value
        const double*       float_value;

        // binary representation of the stored value
        const double_binary_rep*
                            ieee_754;

        // stored value as unsigned integer
        const uint64_t*     int_value;

    public:
        // decode a value val (val cannot be denormal)
        double_decoder(const double* val);

        // get stored float value
        const double&        value() const;

        // get float representation
        const double_binary_rep&
                            get_representation() const;

        // get exponent (bias correction is performed, but raw exponent is not
        // corrected for denormal numbers)
        int                 get_exponent() const;

        // get exponent (bias correction not performed)
        size_t              get_raw_exponent() const;

        // get sign bit
        size_t              get_sign() const;
};

// decode pointers to IEEE 754 double value; denormal numbers are not
// supported;
union float_decoder
{
    public:
        // stored double precision value
        const float*        float_value;

        // binary representation of the stored value
        const float_binary_rep*
                            ieee_754;

        // stored value as unsigned integer
        const uint32_t*     int_value;

    public:
        // decode a value val (val cannot be denormal)
        float_decoder(const float* val);

        // get stored float value
        const float&        value() const;

        // get float representation
        const float_binary_rep&
                get_representation() const;

        // get exponent (bias correction is performed, but raw exponent is not
        // corrected for denormal numbers)
        int                 get_exponent() const;

        // get exponent (bias correction not performed)
        size_t              get_raw_exponent() const;

        // get sign bit
        size_t              get_sign() const;
};

}

#include "matcl-core/details/float/float_binary_rep.inl"
