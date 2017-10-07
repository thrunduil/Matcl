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

#include "matcl-scalar/config.h"

#include "matcl-core/matrix/scalar_types.h"
#include "matcl-core/matrix/enums.h"
#include "matcl-dynamic/details/object.inl"

namespace matcl 
{

// type of object
using Object            = dynamic::object;

// string object
using String            = dynamic::OString;

// object storing an element of type T
using dynamic::object_type;

// bool object
using dynamic::OBool;

// integer object
using dynamic::OInteger;

// double precision floating point value object
using dynamic::OReal;

// single precision floating point value object
using dynamic::OFloat;

// double precision complex floating point value object
using dynamic::OComplex;

// single precision complex floating point value object
using dynamic::OFloat_complex;

// object of string type
using dynamic::OString;

// cast object to Integer
MATCL_SCALAR_EXPORT
Integer                 cast_integer(const Object& obj);

// cast object to Float
MATCL_SCALAR_EXPORT
Float                   cast_float(const Object& obj);

// cast object to Real
MATCL_SCALAR_EXPORT
Real                    cast_real(const Object& obj);

// cast object to Complex
MATCL_SCALAR_EXPORT
Complex                 cast_complex(const Object& obj);

// cast object to Float_complex
MATCL_SCALAR_EXPORT
Float_complex           cast_float_complex(const Object& obj);

// convert Integer to Object of type ty
MATCL_SCALAR_EXPORT
Object                  convert_to_object(const dynamic::Type& ty, Integer v);

// convert Float to Object of type ty
MATCL_SCALAR_EXPORT
Object                  convert_to_object(const dynamic::Type& ty, Float v);

// convert Real to Object of type ty
MATCL_SCALAR_EXPORT
Object                  convert_to_object(const dynamic::Type& ty, Real v);

// convert Complex to Object of type ty
MATCL_SCALAR_EXPORT
Object                  convert_to_object(const dynamic::Type& ty, const Complex& v);

// convert Float_complex to Object of type ty
MATCL_SCALAR_EXPORT
Object                  convert_to_object(const dynamic::Type& ty, const Float_complex& v);

// convert Object to Object of type ty
MATCL_SCALAR_EXPORT
Object                  convert_to_object(const dynamic::Type& ty, const Object& v);

};
