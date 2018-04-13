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

#include "matcl-blas-ext/blas_concurrency.h"

namespace matcl
{

workspace_pool::workspace_pool()
{
    m_mutex_vec.v_.clear();
};

workspace_pool::~workspace_pool()
{};

void workspace_pool::add(void* work)
{
    m_work.push_back(work);
};

workspace_pool::workspace_ptr workspace_pool::get()
{
    void* work = nullptr;

    while(work == nullptr)
    {
        lock_type lock(m_mutex_vec);

        if (this->m_work.size() > 0)
        {
            work = m_work.back();
            m_work.pop_back();
        };
    };

    return workspace_ptr(work,this);
};

void workspace_pool::release(void* work)
{
    lock_type lock(m_mutex_vec);
    m_work.push_back(work);
};

workspace_pool::workspace_ptr::workspace_ptr(void* ptr, workspace_pool* owner)
    :m_ptr(ptr), m_owner(owner)
{};
workspace_pool::workspace_ptr::workspace_ptr()
    :m_ptr(nullptr), m_owner(nullptr)
{};

workspace_pool::workspace_ptr::~workspace_ptr()
{
    if (m_ptr)
        m_owner->release(m_ptr);
};

workspace_pool::workspace_ptr::workspace_ptr(workspace_ptr&& other)
    :m_ptr(other.m_ptr), m_owner(other.m_owner)
{
    other.m_ptr = nullptr;
};
workspace_pool::workspace_ptr& 
workspace_pool::workspace_ptr::operator=(workspace_ptr&& other)
{
    std::swap(m_ptr,other.m_ptr);
    std::swap(m_owner,other.m_owner);
    return *this;
};
void* workspace_pool::workspace_ptr::get() &
{
    return m_ptr;
};
void workspace_pool::workspace_ptr::release()
{
    if (m_ptr)
    {
        m_owner->release(m_ptr);
        m_ptr = nullptr;
    };
};

};