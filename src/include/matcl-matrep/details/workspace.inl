/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/utils/workspace.h"

namespace matcl
{

template<class T>
pod_workspace<T>::pod_workspace()
    :m_ptr(nullptr)
{};

template<class T>
pod_workspace<T>::pod_workspace(size_t size)
    : base_type(size * sizeof(T)) 
{
    m_ptr = static_cast<T*>(base_type::get_ptr());
};

template<class T>
pod_workspace<T>::pod_workspace(size_t size, const T& val)
    : base_type(size * sizeof(T)) 
{
    m_ptr = static_cast<T*>(base_type::get_ptr());
    set_value(val);
}

template<class T>
pod_workspace<T>::pod_workspace(pod_workspace&& other)
    : base_type(std::move(other))
{
    m_ptr = static_cast<T*>(base_type::get_ptr());
};

template<class T>
void pod_workspace<T>::resize(size_t new_size)
{
    base_type::resize(new_size*sizeof(T));
    m_ptr = static_cast<T*>(base_type::get_ptr());
};

template<class T>
void pod_workspace<T>::resize(size_t new_size, const T& val)
{
    base_type::resize(new_size*sizeof(T));
    m_ptr = static_cast<T*>(base_type::get_ptr());
    set_value(val);
};

template<class T>
void pod_workspace<T>::set_value(const T& val)
{
    size_t s = this->size();

    for (size_t i = 0; i < s; ++i)
        ptr()[i] = val;
};

};