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

namespace matcl { namespace dynamic { namespace functions
{

//--------------------------------------------------------------------
//                 predefined function names
//--------------------------------------------------------------------

// (bool)(object); function with this name must take one argument 
// and return OBool
struct MATCL_DYN_EXPORT op_bool     { static const function_name& eval(); };

// OBool operator!(object); function with this name must take one argument 
// and return OBool
struct MATCL_DYN_EXPORT op_not      { static const function_name& eval(); };

// check if value is zero; this function is generated automatically and
// should not be registered explicitly
struct MATCL_DYN_EXPORT is_zero     { static const function_name& eval(); };

// check if value is one; this function is generated automatically and
// should not be registered explicitly
struct MATCL_DYN_EXPORT is_one      { static const function_name& eval(); };

// real part of complex number
struct MATCL_DYN_EXPORT real        { static const function_name& eval(); };

// imaginary part of complex number
struct MATCL_DYN_EXPORT imag        { static const function_name& eval(); };

// OBool operator==(object,object)
struct MATCL_DYN_EXPORT op_eeq      { static const function_name& eval(); };

// operator!=(object,object)
struct MATCL_DYN_EXPORT op_neq      { static const function_name& eval(); };

// OBool operator<(object,object)
struct MATCL_DYN_EXPORT op_lt       { static const function_name& eval(); };

// OBool operator<=(object,object)
struct MATCL_DYN_EXPORT op_leq      { static const function_name& eval(); };

// OBool operator>=(object,object)
struct MATCL_DYN_EXPORT op_geq      { static const function_name& eval(); };

// OBool operator>(object,object)
struct MATCL_DYN_EXPORT op_gt       { static const function_name& eval(); };

// operator-(object)
struct MATCL_DYN_EXPORT op_uminus   { static const function_name& eval(); };

// operator+(object,object)
struct MATCL_DYN_EXPORT op_plus     { static const function_name& eval(); };

// operator-(object,object)
struct MATCL_DYN_EXPORT op_minus    { static const function_name& eval(); };

// operator*(object,object)
struct MATCL_DYN_EXPORT op_mul      { static const function_name& eval(); };

// name of function operator/(object,object), for objects
// representing integer values this should give values as close
// as possible to value when integers are conveted to floats;
// integer division is represented by function op_idiv; 
// function with this name must take two arguments
struct MATCL_DYN_EXPORT op_div      { static const function_name& eval(); };

// division, equivalent to op_div for objects not representing integers
struct MATCL_DYN_EXPORT idiv        { static const function_name& eval(); };

};};};
