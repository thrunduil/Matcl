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

#include "matcl-dynamic/function_name.h"
#include "matcl-dynamic/predefined_functions_names.h"

namespace matcl { namespace dynamic
{

//--------------------------------------------------------------------
//               predefined functions for object
//--------------------------------------------------------------------
// these functions are registered for predefined object types in
// matcl-dynamic

// explicitly cast object obj to object of type new_type; a registered
// promotion, equivalent, decay, explicit, or cast converter is called
MATCL_DYN_EXPORT object cast(Type new_type, const object& obj);

// explicitly cast object obj to object of type new_type; a registered
// promotion, equivalent, decay or explicit converter is called
MATCL_DYN_EXPORT object convert(Type new_type, const object& obj);

// cast to boolean value
MATCL_DYN_EXPORT bool   cast_bool(const object& x);

// logical negation; should be equivalent to !cast_bool(a)
MATCL_DYN_EXPORT bool   operator!(const object& x);

// check if value is zero
MATCL_DYN_EXPORT bool   is_zero(const object& x);

// check if value is one
MATCL_DYN_EXPORT bool   is_one(const object& x);

// real part of a complex number
MATCL_DYN_EXPORT object real(const object& x);

// imaginary part of a complex number
MATCL_DYN_EXPORT object imag(const object& x);

// comparison functions
MATCL_DYN_EXPORT object operator==(const object& x, const object& y);
MATCL_DYN_EXPORT object operator!=(const object& x, const object& y);
MATCL_DYN_EXPORT object operator<(const object& x, const object& y);
MATCL_DYN_EXPORT object operator<=(const object& x, const object& y);
MATCL_DYN_EXPORT object operator>=(const object& x, const object& y);
MATCL_DYN_EXPORT object operator>(const object& x, const object& y);

// arithmetic operators
MATCL_DYN_EXPORT object operator+(const object& x, const object& y);
MATCL_DYN_EXPORT object operator-(const object& x, const object& y);
MATCL_DYN_EXPORT object operator*(const object& x, const object& y);

// unary minus
MATCL_DYN_EXPORT object operator-(const object& x);

// for objects representing integer values this should give values as 
// close as possible to value when integers are conveted to floats;
// for integer division one should use idiv function
MATCL_DYN_EXPORT object operator/(const object& x, const object& y);

// division, equivalent to operator/ for objects not representing 
// integers;
MATCL_DYN_EXPORT object idiv(const object& x, const object& y);

};};

