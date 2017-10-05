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
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/matrix/enums.h"

namespace matcl { namespace constants
{

//-----------------------------------------------------------------
//                      integer constants
//-----------------------------------------------------------------
// maximum value of Integer
Integer                 max_int();

// minumum value of Integer
Integer                 min_int();

//-----------------------------------------------------------------
//                      float constants
//-----------------------------------------------------------------
// machine epsilon; distance from 1.0 to the next larger single precision number
Float                   f_eps();

// smallest positive normalized floating point number.
Float                   f_min_real();

// largest positive normalized floating point number.
Float                   f_max_real();

// infinite value
Float                   f_inf();

// not-a-number value
Float                   f_nan();

// e value; e = exp(1)
Float                   f_e();

// pi value
Float                   f_pi();

// pi/2 value
Float                   f_pi_2();

// the logarithm to base 2 of f_e
Float                   f_log2e();

// the logarithm to base 10 of f_e
Float                   f_log10e();

// the logarithm of 2
Float                   f_ln2();

// the logarithm of 10
Float                   f_ln10();

// the square root of 2
Float                   f_sqrt2();

// the square root of 1/2
Float                   f_sqrt1_2();

//-----------------------------------------------------------------
//                      double constants
//-----------------------------------------------------------------
// machine epsilon; distance from 1.0 to the next larger double precision number
Real                    eps();

// smallest positive normalized floating point number.
Real                    min_real();

// largest positive normalized floating point number.
Real                    max_real();

// infinite value
Real                    inf();

// not-a-number value
Real                    nan();

// e value; e = exp(1)
Real                    e();

// pi value
Real                    pi();

// pi/2 value
Real                    pi_2();

// the logarithm to base 2 of e
Real                    log2e();

// the logarithm to base 10 of e
Real                    log10e();

// the logarithm of 2
Real                    ln2();

// the logarithm of 10
Real                    ln10();

// the square root of 2
Real                    sqrt2();

// the square root of 1/2
Real                    sqrt1_2();

//-----------------------------------------------------------------
//                      complex constants
//-----------------------------------------------------------------
// imaginary unit
Float_complex           f_i();
// imaginary unit
Complex                 i();

//-----------------------------------------------------------------
//                      templated constants
//-----------------------------------------------------------------
// template argument Value is Float or Real

// machine epsilon; distance from 1.0 to the next larger number
template<class Value> 
Value                   eps();

// smallest positive normalized floating point number
template<class Value> 
Value                   min_real();

// largest positive normalized floating point number
template<class Value> 
Value                   max_real();

// infinite value
template<class Value>
Value                   inf();

// not-a-number value
template<class Value> 
Value                   nan();

// e value; e = exp(1)
template<class Value>
Value                   e();

// pi value
template<class Value>
Value                   pi();

// pi/2 value
template<class Value>
Value                   pi_2();

// the logarithm to base 2 of e
template<class Value>
Value                   log2e();

// the logarithm to base 10 of e
template<class Value> 
Value                   log10e();

// the logarithm of 2
template<class Value>
Value                   ln2();

// the logarithm of 10
template<class Value> 
Value                   ln10();

// the square root of 2
template<class Value>
Value                   sqrt2();

// the square root of 1/2
template<class Value>
Value                   sqrt1_2();

// imaginary unit
template<class Value>
complex<Value>          i();

//-----------------------------------------------------------------
//      constants for given value code converted to Real
//-----------------------------------------------------------------
// machine epsilon; distance from 1.0 to the next larger number
MATCL_CORE_EXPORT Real  eps(value_code vc);

// smallest positive normalized floating point number
MATCL_CORE_EXPORT Real  min_real(value_code vc);

// largest positive normalized floating point number
MATCL_CORE_EXPORT Real  max_real(value_code vc);

};};

#include "matcl-core/details/constants.inl"
