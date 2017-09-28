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

#include "matcl-core/details/hash_table/object_table.h"

namespace matcl { namespace details
{

template<class Main_allocator>
inline pool_allocator<Main_allocator>::pool_allocator(size_t size)
    :m_pool(size)
{};

template<class Main_allocator>
inline void pool_allocator<Main_allocator>::free(void* ptr,size_t)
{
    #ifdef _DEBUG
        assertion(m_pool.is_from(ptr) == true, "invalid free");
    #endif
    m_pool.free(ptr);
};

template<class Main_allocator>
inline void* pool_allocator<Main_allocator>::malloc_impl()
{
    return m_pool.malloc();
};

template<class Main_allocator>
inline void pool_allocator<Main_allocator>::purge_memory()
{
    m_pool.purge_memory();
};

template<class V, class H, class E, class A>
object_table<V,H,E,A>::object_table(size_t capacity)
:base_type(sizeof(VT)), m_table(capacity)
{};

template<class V, class H, class E, class A>
object_table<V,H,E,A>::~object_table()
{
    close(true);
};

template<class V, class H, class E, class A>
void object_table<V,H,E,A>::close(bool call_destructors)
{
    m_table.close(call_destructors);
    base_type::purge_memory();
};

template<class V, class H, class E, class A>
template<class Y>
V object_table<V,H,E,A>::get(const Y& str)
{
    using entry = hash_table::entry;	
    entry ptr   = m_table.get(str);

    if (!ptr.empty())
    {
        V::type_traits::copy(*ptr);
        return V::make(*ptr);
    }
    else
    {
        VT* m_str	= register_obj(str);
        ptr.assign(m_str);
        return V::make(*ptr);
    };
};

template<class V, class H, class E, class A>
template<class Y>
V object_table<V,H,E,A>::get_existing(const Y& str) const
{
    using const_entry = hash_table::const_entry;	
    const_entry ptr   = m_table.get(str);

    if (!ptr.empty())
    {
        V::type_traits::copy(*ptr);
        return V::make(*ptr);
    }
    else
    {
        return V::make(nullptr);
    };
};

template<class V, class H, class E, class A>
template<class Y>
typename object_table<V,H,E,A>::VT* object_table<V,H,E,A>::register_obj(Y&& val)
{   
    size_t b = allocator::get_bytes_info(val);
    void* ptr  = base_type::malloc(b);
    new(ptr) VT(std::forward<Y>(val));
    return reinterpret_cast<VT*>(ptr);
};

template<class V, class H, class E, class A>
void object_table<V,H,E,A>::unregister_obj(VT* ptr)
{
    m_table.remove(ptr);
    size_t b = allocator::get_bytes(ptr);
    ptr->~VT();	
    base_type::free(const_cast<void*>(static_cast<const void*>(ptr)),b);
};

template<class V, class H, class E, class A>
template<class stack_type>
void object_table<V,H,E,A>::unregister_obj(VT* ptr, stack_type& st)
{
    using VT_nc = std::remove_const<VT>::type;

    m_table.remove(ptr);
    size_t b = allocator::get_bytes(ptr);
    const_cast<VT_nc*>(ptr)->destroy(st);	
    base_type::free(const_cast<void*>(static_cast<const void*>(ptr)),b);
};

template<class V, class H, class E, class A>
double object_table<V,H,E,A>::reuse_stats() const
{
    size_t N = m_table.capacity();
    const VT* const* ptr = m_table.get_entries();

    double M = 0;
    double K = 0;
    for (size_t i = 0; i < N; ++i)
    {
        const VT* elem = ptr[i];
        if (elem > (const VT*)1)
        {
            K += elem->refcount();
            M += 1;
        };
    };
    return K/(M+1e-5);
};

};};