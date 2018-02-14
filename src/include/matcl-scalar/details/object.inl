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

#include "matcl-scalar/object.h"
#include "matcl-core/general/type_traits.h"

namespace matcl { namespace details
{

template<class Type, 
    bool Is_extern_scalar = is_external_scalar<typename std::decay<Type>::type>::value>
struct make_object_impl
{
    static Object eval(Type&& t)
    {
        using base_type     = typename std::decay<Type>::type;
        return Object(dynamic::object_type<base_type>(std::forward<Type>(t)));
    }
};

template<class Type>
struct make_object_impl<Type, false>
{
    static Object eval(Type&& t)
    {
        return Object(std::forward<Type>(t));
    }
};

}}

namespace matcl 
{

template<class T>
Object matcl::make_object(T&& val)
{
    return details::make_object_impl<T>::eval(std::forward<T>(val));
}

};
