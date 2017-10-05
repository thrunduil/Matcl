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

#include "matcl-core/memory/global_objects.h"
#include "matcl-core/details/utils.h"
#include <string>

namespace matcl
{

// these templates should be specialized for scalar types defined in other
// libraries, in order to enable some of functions defined in matcl

// mark type T as a scalar type
template<class T>
struct is_external_scalar
{ 
    static const bool value = false; 
};

template<>
struct is_external_scalar<std::string>
{ 
    static const bool value = true; 
};

// get complex type associated with given type; it T is already a complex
// type, then set type to T
template<class T>
struct make_complex_type
{ 
    using type = typename details::complex_type<T>::type; 
};

};