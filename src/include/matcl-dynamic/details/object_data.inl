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

#include "matcl-dynamic/details/object_data.h"
#include "matcl-dynamic/details/object_data_pool.h"
#include "matcl-dynamic/details/register_object.inl"

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 
#pragma warning(disable:4505)   //unreferenced local function has been removed

namespace matcl { namespace dynamic { namespace details
{

force_inline
void object_data_base::increase_refcount() const
{
    ++m_refcount;
};

force_inline
void object_data_base::decrease_refcount()
{
    if((--m_refcount ) == 0)
        delete this;
};

template<class T>
force_inline
const T& object_data_base::get_value() const
{ 
    return static_cast<const object_data<T>*>(this)->get();
};

template<class T>
force_inline
T& object_data_base::get_value()
{ 
    return static_cast<object_data<T>*>(this)->get();
};

force_inline
object_data_base::~object_data_base()
{};

force_inline
bool object_data_base::is_unique() const
{ 
    return m_refcount == 1; 
};

force_inline
object_data_base::object_data_base()
    :m_refcount(0)
{};


template<class T>
force_inline
object_data<T>* object_data<T>::create(const T& val)
{
    return new object_data(val);
}

template<class T>
force_inline
object_data<T>* object_data<T>::create(T&& val)
{
    return new object_data(std::move(val));
}

template<class T>
force_inline
void* object_data<T>::operator new(size_t)
{
    return object_data_pool<T>::malloc();
};

template<class T>
force_inline
void object_data<T>::operator delete(void* ptr)
{
    object_data_pool<T>::free(ptr);
};

template<class T>
force_inline
object_data<T>::object_data()
    : m_data(T())
{};

template<class T>
force_inline
object_data<T>::object_data(const T& val)
    : m_data(val)
{};

template<class T>
force_inline
object_data<T>::object_data(T&& val)
    : m_data(std::move(val))
{};

template<class T>
force_inline
const typename object_data<T>::value_type& 
object_data<T>::get() const
{ 
    return m_data; 
};

template<class T>
force_inline
typename  object_data<T>::value_type& object_data<T>::get()
{ 
    return m_data; 
};

};};};

#pragma warning(pop)
