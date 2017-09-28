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

namespace boost
{

template <typename UserAllocator>
class pool;

struct default_user_allocator_malloc_free;
};

#pragma warning(push)
#pragma warning(disable:4251)	//needs to have dll-interface 

namespace matcl { namespace dynamic { namespace details
{

struct item_counter
{
    private:
        using ATOMIC_LONG = matcl::atomic<long>;

    private:
        static MATCL_DYN_EXPORT ATOMIC_LONG n_item;

    public:
        static long get()		{ return n_item; };
        static void inc()		{ ++n_item;	};
        static void dec()		{ --n_item;	};
};

class MATCL_DYN_EXPORT object_data_pool_base
{
    private:
        using pool_type         = boost::pool<boost::default_user_allocator_malloc_free>;
        pool_type*				m_pool;

        object_data_pool_base(const object_data_pool_base&) = delete;
        object_data_pool_base& operator=(const object_data_pool_base&) = delete;

    protected:
        object_data_pool_base(size_t size);

        void	free(void* ptr);
        void*	malloc();

        ~object_data_pool_base();
};

template<class T>
struct object_data_pool_impl : protected object_data_pool_base
{
    private:
        using base          = object_data_pool_base;
        using mutex_type    = matcl::default_spinlock_mutex;
        using lock_type     = std::unique_lock<mutex_type>;

        mutex_type          m_mutex;

        object_data_pool_impl()
            :base(sizeof(object_data<T>))
        {};

        object_data_pool_impl(const object_data_pool_impl&) = delete;
        object_data_pool_impl& operator=(const object_data_pool_impl&) = delete;

        ~object_data_pool_impl();

        void* malloc_impl()
        {
            lock_type lock(m_mutex);
            item_counter::inc();
            return base::malloc();
        };
        void free_impl(void* ptr)
        {
            lock_type lock(m_mutex);
            item_counter::dec();
            base::free(ptr);
        };

        static object_data_pool_impl* get()
        {
            static object_data_pool_impl* instance = new object_data_pool_impl<T>();
            return instance;
        }

    public:
        long no_existing_objects()
        {
            return item_counter::n_item;
        };

        static void* malloc()
        {
            return get()->malloc_impl();
        };

        static void free(void* ptr)
        {
            return get()->free_impl(ptr);
        };
};

};};};

#pragma warning(pop)