/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/details/refcount.h"
#include <boost/pool/pool.hpp>
#include "thread_local_cache.h"
#include "matcl-core/error/exception_classes.h"

namespace matcl { namespace details
{

template<class refcount_pool>
struct refstruct_pool_impl : protected boost::pool<boost::default_user_allocator_new_delete>
{
    private:
        static const int cache_cap      = refcount_pool::cache_capacity;

    private:
        using refcount_type = refcount_str<refcount_pool>;
        using cache_type    = cache<refcount_pool, cache_cap>;
        using item_counter  = typename refcount_pool::item_counter;

        using base          = boost::pool<boost::default_user_allocator_new_delete>;
        using spin_lock     = matcl::default_spinlock_mutex;
        using scoped_lock   = std::unique_lock<spin_lock>;

        spin_lock       m_mutex;
        item_counter    n_item;

    public:
        static refstruct_pool_impl*     m_pool;

    public:
        refstruct_pool_impl()
            :base(sizeof(refcount_type))
        {
            n_item.reset();
        };

        void* malloc()
        {
            n_item.increase();
            void* ptr = cache_type::get().malloc(this);
            return ptr;
        };

        void* hard_malloc()
        {
            scoped_lock lock(m_mutex);            
            return base::malloc();
        };

        void free(void* ptr)
        {
            n_item.decrease();
            cache_type::get().free(ptr, this);
        }

        void hard_free(void* ptr)
        {
            scoped_lock lock(m_mutex);            
            base::free(ptr);
        };

        long no_existing_objects()
        {
            return n_item.count();
        };
};

//refstruct_pool_impl must be alive while destructors of static objects are called
//does not create important memleaks

template<class refcount_pool>
refstruct_pool_impl<refcount_pool>* refstruct_pool_impl<refcount_pool>::m_pool 
                = new refstruct_pool_impl<refcount_pool>();

template<class refcount_pool>
refcount_str<refcount_pool>* refstruct_pool<refcount_pool>::create(long val)
{
    using ptr_type  = refcount_str<refcount_pool>;
    ptr_type* ptr   = (ptr_type*) refstruct_pool_impl<refcount_pool>
                                        ::m_pool->malloc();

    new(ptr) ptr_type(val);
    return ptr;
};

template<class refcount_pool>
void refstruct_pool<refcount_pool>::destroy(refcount_str<refcount_pool>* ptr)
{
    return refstruct_pool_impl<refcount_pool>::m_pool->free(ptr);
};

long no_existing_objects()
{
    return refstruct_pool_impl<default_refcount_pool>::m_pool->no_existing_objects();
};

void close_refcount()
{
    delete refstruct_pool_impl<default_refcount_pool>::m_pool;
    refstruct_pool_impl<default_refcount_pool>::m_pool = nullptr;
};

template struct refstruct_pool_impl<default_refcount_pool>;
template struct refstruct_pool_impl<scatter_refcount_pool>;

template class refstruct_pool<default_refcount_pool>;
template class refstruct_pool<scatter_refcount_pool>;

};};
