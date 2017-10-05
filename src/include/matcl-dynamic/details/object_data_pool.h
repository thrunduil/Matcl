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
#include "matcl-core/details/hash_table/object_table.h"
#include "matcl-core/memory/alloc.h"

namespace matcl { namespace dynamic { namespace details
{

namespace md = matcl :: details;

// object_data pool must be shared between dlls
class object_data_pool_impl
    : protected md::object_allocator<md::default_allocator_simple<true, true, char>>
{
    private:
        using allocator     = md::default_allocator_simple<true, true, char>;
        using base          = md::object_allocator<allocator>;
        using this_type     = object_data_pool_impl;
        using mutex_type    = matcl::default_spinlock_mutex;
        using lock_type     = std::unique_lock<mutex_type>;

        mutex_type          m_mutex;

    private:
        object_data_pool_impl(bool is_global, size_t size)
            :base(size, is_global)
        {};

        object_data_pool_impl(const object_data_pool_impl&) = delete;
        object_data_pool_impl& operator=(const object_data_pool_impl&) = delete;

        ~object_data_pool_impl()
        {};

        friend class global_objects;

    public:
        void* malloc()
        {
            lock_type lock(m_mutex);
            return base::malloc();
        };

        void free(void* ptr)
        {
            lock_type lock(m_mutex);
            base::free(ptr);
        };
};

template<class T>
struct object_data_pool
{
    private:
        using allocator     = object_data_pool_impl;

    private:
        
        static object_data_pool_impl* get()
        {            
            return register_object<T>::g_hook.get_pool();
        }

    public:
        static void* malloc()
        {
            return get()->malloc();
        };

        static void free(void* ptr)
        {
            return get()->free(ptr);
        };
};

};};};
