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

#include "matcl-dynamic/special_types.h"

namespace matcl { namespace dynamic
{

inline const object& any_type::get_stored() const
{
    return m_stored;
};

template<class T>
any_type::any_type(const object_type<T>& other)
    :m_stored(other)
{};

template<class T>
any_type::any_type(object_type<T>&& other)
    :m_stored(std::move(other))
{};

template<class S, class Enable>
any_type::any_type(S&& other)
    :m_stored(object_type<typename std::decay<S>::type>(std::forward<S>(other)))
{};

template<class T>
any_type& any_type::operator=(const object_type<T>& other) &
{
    m_stored.reset(object(other));
    return *this;
};

template<class T>
any_type& any_type::operator=(object_type<T>&& other) &
{
    m_stored.reset(object(std::move(other)));
    return *this;
};

template<class S, class Enable>
any_type& any_type::operator=(S&& other) &
{
    using S0 = typename std::decay<S>::type;
    m_stored.reset(object(object_type<S0>(std::forward<S>(other))));
    return *this;
};

//---------------------------------------------------------------
//                  reference type
//---------------------------------------------------------------
template<class T>
Type object_reference<T>::get_static_type()
{
    return dynamic::details::mark_reference_type<T>().get();
};

};};
