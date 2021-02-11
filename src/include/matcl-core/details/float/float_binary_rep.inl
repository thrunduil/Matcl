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

#include "matcl-core/float/float_binary_rep.h"

namespace matcl { namespace details
{

union convert_int_float
{
    uint32_t    m_int;
    float       m_float;
};

union convert_int_double
{
    uint64_t    m_int;
    double      m_float;
};

}}

namespace matcl
{

force_inline 
float matcl::hex_float(const uint32_t& val)
{
    return reinterpret_cast<const details::convert_int_float*>(&val)->m_float;
};

force_inline double matcl::hex_double(const uint64_t& val)
{
    return reinterpret_cast<const details::convert_int_double*>(&val)->m_float;
};

//-----------------------------------------------------------------------
//                   DOUBLE PRECISION DECODER
//-----------------------------------------------------------------------
force_inline
double_decoder::double_decoder(const double* val)
    :float_value(val)
{};

force_inline
const double& double_decoder::value() const
{ 
    return *float_value;
}

force_inline
const double_binary_rep& double_decoder::get_representation() const
{
    return *ieee_754;
};

force_inline
int double_decoder::get_exponent() const
{ 
    return int(ieee_754->exponent) - double_binary_rep::bias; 
};

force_inline
size_t double_decoder::get_raw_exponent() const
{ 
    return ieee_754->exponent; 
};

force_inline
size_t double_decoder::get_sign() const
{ 
    return ieee_754->sign; 
};

//-----------------------------------------------------------------------
//                   SINGLE PRECISION DECODER
//-----------------------------------------------------------------------
force_inline
float_decoder::float_decoder(const float* val)
    :float_value(val)
{};

force_inline
const float& float_decoder::value() const
{ 
    return *float_value;
}

force_inline
const float_binary_rep& float_decoder::get_representation() const
{
    return *ieee_754;
};

force_inline
int float_decoder::get_exponent() const
{ 
    return int(ieee_754->exponent) - double_binary_rep::bias; 
};

force_inline
size_t float_decoder::get_raw_exponent() const
{ 
    return ieee_754->exponent; 
};

force_inline
size_t float_decoder::get_sign() const
{ 
    return ieee_754->sign; 
};

}
