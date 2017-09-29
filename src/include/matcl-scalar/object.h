/* 
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2017
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
#include "matcl-dynamic/object.h"

namespace matcl { namespace details
{

class printer;

}};

namespace matcl 
{

// type of object
using Object            = dynamic::object;

// object storing an element of type T
template<class T>
using object_type       = dynamic::object_type<T>;

// string object
using String            = dynamic::object_type<std::string>;

// bool object
using OBool             = dynamic::object_type<bool>;

// integer object
using OInteger          = dynamic::object_type<Integer>;

// double precision floating point value object
using OReal             = dynamic::object_type<Real>;

// single precision floating point value object
using OFloat            = dynamic::object_type<Float>;

// double precision complex floating point value object
using OComplex          = dynamic::object_type<Complex>;

// single precision complex floating point value object
using OFloat_complex    = dynamic::object_type<Float_complex>;

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

// implement disp function for objects; internal use only
MATCL_SCALAR_EXPORT
void                    disp_object(const Object& obj, details::printer& pr, Integer elem_width, 
                            matcl::align_type at, Integer value_pos);

// get string representation of object v
MATCL_SCALAR_EXPORT
std::string             to_string(const Object& v);

};
