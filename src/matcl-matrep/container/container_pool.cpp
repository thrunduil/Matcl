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

#include "matcl-matrep/container/container_pool.h"
#include "matcl-matrep/container/matrix_container.h"
#include "matcl-matrep/base/thread_local_cache.h"
#include <boost/pool/pool.hpp>
#include "matcl-matrep/details/extract_type_switch.h"

namespace matcl { namespace details
{

#if (DEBUG_MEMORY == 1)
  using item_counter = default_item_counter;
#else
  using item_counter = empty_item_counter;
#endif

template<int code, int size>
class pool_impl : private boost::pool<boost::default_user_allocator_new_delete>
{
    private:
        static const int cache_cap      = 30;

    private:
        using base_type     = boost::pool<boost::default_user_allocator_new_delete>;
        using cache_type    = cache<pool_impl, cache_cap>;
        using spin_lock     = matcl::default_spinlock_mutex;
        using scoped_lock   = std::unique_lock<spin_lock>;

        spin_lock           m_mutex;
        item_counter        n_item;
        static pool_impl    g_instance;

    public:
        pool_impl()
            :base_type(size)
        {
            n_item.reset();
        };

        static pool_impl& get()  { return g_instance; };

        void destroy()
        {
            this->~pool_impl();
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
            return base_type::malloc();
        };
        
        void free(void* ptr)
        {
            n_item.decrease();
            cache_type::get().free(ptr, this);
        };
        
        void hard_free(void*ptr)
        {
            scoped_lock lock(m_mutex);            
            base_type::free(ptr);
        };
};

template<int code, int size>
pool_impl<code,size> pool_impl<code,size>::g_instance;

template<class T>
struct pool_size_impl
{
    static const size_t value = 4;
};

template<class V, class S>
struct pool_size_impl<raw::Matrix<V,S>>
{
    using container_type = matrix_container<V,S>;
    static const size_t value = sizeof(container_type);
};

template<int code>
struct pool_size
{
    using matrix_type = typename code_to_type<(matcl::mat_code)code>::type;
    static const size_t value = pool_size_impl<matrix_type>::value;
};

template<Integer start, Integer end, template<Integer code> class functor>
struct For
{
    template<class ... Arg>
    static void eval(Arg&& ... args)
    {
        functor<start>::eval(std::forward<Arg>(args)...);
        For<start+1,end,functor>::eval<Arg...>(std::forward<Arg>(args)...);
    };
};

template<Integer end, template<Integer code> class functor>
struct For<end,end,functor>
{
    template<class ... Arg>
    static void eval(Arg&& ... ){};
};

template<Integer code>
struct close_pool
{
    static void eval()
    {
        static const int size = pool_size<code>::value;
        return pool_impl<code,size>::get().destroy();
    };
};

void container_pool::close()
{
    For<macro_first_matrix_type_code, macro_last_matrix_type_code, close_pool>::eval();
};

template<int code>
void* container_pool::malloc()
{
    static const int size = pool_size<code>::value;
    return pool_impl<code,size>::get().malloc();
};

struct pool_visitor : public extract_type_from_code<void,pool_visitor>
{
    template<class Mat_type>
    static void eval(void* ptr)
    {
        static const int code = (int)type_to_code<Mat_type>::value;
        static const int size = pool_size<code>::value;
        return pool_impl<code,size>::get().free(ptr);
    };

    template<class Mat_type>
    static void eval_scalar(void*)
    {};
};

void container_pool::free(void* ptr, int code)
{
    return pool_visitor::make((mat_code)code, ptr);
}

template<class Matrix_Type>
struct inst_pool
{
    static const int code = (int)type_to_code<Matrix_Type>::value;
    
    void* malloc_inst()
    {
        return container_pool::malloc<code>();
    };
};

#define MACRO_INST_POOLS(arg,d1,d2) template struct inst_pool<arg>;

MACRO_FOREACH_CODE_ALL(MACRO_INST_POOLS,0,0)    

};};
