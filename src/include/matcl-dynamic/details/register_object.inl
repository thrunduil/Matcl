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

#include "matcl-core/general/fwd_decls.h"
#include "matcl-dynamic/details/global_objects.h"
#include "matcl-dynamic/details/type_object.inl"

namespace matcl { namespace dynamic { namespace details
{

template<class T> 
std::string mark_type<T>::name() const
{ 
    return details::register_object_helper::get_name(typeid(T).name()); 
};

template<class T>
details::register_object_helper 
register_object<T>::g_hook(typeid(T).name(), details::register_object_helper::create_data<T>);

template<class T> 
Type mark_type<T>::get() const
{    
    //static Type ti( register_object<T>::g_hook.get_object(typeid(T).name()) );
    return register_object<T>::g_hook.get_type();
};

template<class T>
Type details::register_object_helper::create_data(pool_type*& pool)
{     
    return Type(global_objects::initialize_type<T>(pool)); 
};

inline Type details::register_object_helper::get_type() const
{
    return m_type;
};

inline details::register_object_helper::pool_type* 
details::register_object_helper::get_pool() const
{
    return m_pool;
};

force_inline
Type mark_type<null_type>::get() const
{    
    return Type();
};

inline std::string mark_type<null_type>::name() const
{ 
    return "null"; 
};

template<class T> 
Type mark_reference_type<T>::get() const
{
    static Type g_ti = operations::make_reference_type(T::get_static_type());
    return g_ti;
};

};};};
