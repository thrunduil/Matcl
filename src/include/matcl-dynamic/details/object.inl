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

#include "matcl-dynamic/object.h"
#include "matcl-dynamic/details/object_data.inl"

namespace matcl { namespace dynamic
{

//------------------------------------------------------------
//                      object
//------------------------------------------------------------
force_inline object::object()
    :m_data(nullptr)
{};

force_inline object::object(const object& other)    
    :m_type(other.m_type), m_data(other.m_data)
{
    if (m_data)
        m_data->increase_refcount();
};

force_inline object::object(object&& other)
    :m_type(std::move(other.m_type)), m_data(other.m_data)
{
    other.m_data = nullptr;
};

force_inline object::object(Type ti, data_type* data, not_null)
    :m_type(ti), m_data(data)
{};

force_inline object::object(Type ti, null)
    :m_type(ti), m_data(nullptr)
{
};

force_inline object::~object()
{
    if (m_data)
        m_data->decrease_refcount();
};

force_inline void object::reset(const object& other) &
{
    data_type* ptr = other.m_data;

    if (ptr)
        ptr->increase_refcount();

    if(m_data != nullptr)
        m_data->decrease_refcount();

    m_data  = ptr;
    m_type  = other.m_type;
    return;
};

force_inline void object::reset(object&& other) &
{
    object tmp(std::move(other));
    swap(*this, tmp);
    return;
};

// conversion from typed object
template<class T>
force_inline
object::object(const object_type<T>& other)
    :object(other.m_data)
{};

// conversion from temporary typed object
template<class T>
force_inline
object::object(object_type<T>&& other)
    :object(std::move(other.m_data))
{};

force_inline bool object::is_unique() const
{
    if (m_data == nullptr || m_data->is_unique())
        return true;
    else 
        return false;
}

force_inline void swap(object& o1, object& o2)
{
    std::swap(o1.m_type,o2.m_type);
    std::swap(o1.m_data,o2.m_data);
};

force_inline
object& object::operator=(const object& other) &
{    
    Type in1                = this->get_type();
    Type in2                = other.get_type();

    if (operations::has_trivial_assignment(in1, in2) == true)
        reset(other);        
    else
        make_assignment(other, in1, in2);
    
    return *this;
};

force_inline
object& object::operator=(object&& other) &
{
    Type in1                = this->get_type();
    Type in2                = other.get_type();

    if (operations::has_trivial_assignment(in1, in2) == true)
        reset(std::move(other));        
    else
        make_assignment(std::move(other), in1, in2);
    
    return *this;
};

};};
