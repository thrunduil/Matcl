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

#include "matcl-dynamic/details/data_constructor.h"
#include "matcl-dynamic/dynamic_function/details/type_list.h"

namespace matcl { namespace dynamic { namespace details
{

template<class FuncTraits>
template<int N, class T>
typename enable_const<T>::type
    data_constructor<FuncTraits>::get() const
{
    using in_type       = typename get_elem<typename FuncTraits::input_type,N>::type;
    using value_type    = typename std::decay<T>::type;

    static_assert(std::is_same<T,in_type>::value,"incompatible types");

    return *reinterpret_cast<const value_type*>(m_args[N]);
};

template<class FuncTraits>
template<int N, class T>
typename enable_nconst<T>::type
    data_constructor<FuncTraits>::get() const
{
    using in_type       = typename get_elem<typename FuncTraits::input_type,N>::type;
    using value_type    = typename std::decay<T>::type;
    using data_type     = object_data<value_type>;

    static_assert(std::is_same<T,in_type>::value,"incompatible types");

    const_cast<object*>(m_args[N])->make_unique();

    return *reinterpret_cast<value_type*>(const_cast<object*>(m_args[N]));
};

};};};
