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

#include "matcl-scalar/object.h"
#include "matcl-scalar/objects/object_functions.h"
#include "matcl-dynamic/predefined_functions.h"
#include "matcl-dynamic/object_type.h"
#include "matcl-core/details/IO/printer.h"
#include "matcl-core/options/options_disp.h"
#include "matcl-scalar/IO/scalar_io.h"

namespace matcl 
{

namespace md = matcl::details;

Integer matcl::cast_integer(const Object& v)
{
    dynamic::Type t = OInteger::get_static_type();
    return OInteger(cast(t, v), dynamic::from_object()).get();
};

Real matcl::cast_real(const Object& v)
{
    dynamic::Type t = OReal::get_static_type();
    return OReal(cast(t, v), dynamic::from_object()).get();
};

Float matcl::cast_float(const Object& v)
{
    dynamic::Type t = OFloat::get_static_type();
    return OFloat(cast(t, v), dynamic::from_object()).get();
};

Complex matcl::cast_complex(const Object& v)
{
    dynamic::Type t = OComplex::get_static_type();
    return OComplex(cast(t, v), dynamic::from_object()).get();
};

Float_complex matcl::cast_float_complex(const Object& v)
{
    dynamic::Type t = OFloat_complex::get_static_type();
    return OFloat_complex(cast(t, v), dynamic::from_object()).get();
};

Object matcl::convert_to_object(const dynamic::Type& ty, Integer v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, Float v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, Real v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, const Complex& v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, const Float_complex& v)
{
    return Object(matcl::dynamic::convert(ty, Object(v)));
}

Object matcl::convert_to_object(const dynamic::Type& ty, const Object& v)
{
    return Object(matcl::dynamic::convert(ty, v));
}

};
