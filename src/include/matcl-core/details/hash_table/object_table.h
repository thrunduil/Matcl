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

#include <boost/pool/pool.hpp>
#include "matcl-core/details/hash_table/hash_table.h"

namespace matcl { namespace details
{

//TODO: take from sym_arrow

#pragma warning (push)
#pragma warning (disable:4251)

class default_allocator
{
    public:
        using size_type         = std::size_t;
        using difference_type   = std::ptrdiff_t;
  
        static char*    malloc(size_type bytes) { return static_cast<char *>(std::malloc(bytes)); }
        static void     free(void* block)       { std::free(block); }
};

template<class Main_allocator>
class pool_allocator
{
    private:
        using pool          = boost::pool<Main_allocator>;

    private:
        pool                m_pool;

    public:
        pool_allocator(size_t size);

        void                free(void* ptr, size_t b);
        void*               malloc(size_t)              { return malloc_impl(); };
        void                purge_memory();

        static size_t       get_bytes(const void*)      { return size_t(); }
        template<class Y>
        static size_t       get_bytes_info(const Y&)    { return size_t(); }

    private:
        void*               malloc_impl();
};

#pragma warning (pop)


template<class object_ptr, class hasher, class equaler, class main_alloc>
class object_table : public pool_allocator<main_alloc>
{
    private:
        using VT0           = typename object_ptr::value_type;
        using VT            = typename std::remove_pointer<VT0>::type;
        using track_value   = default_track_value<VT>;
        using hash_table    = hash_table<VT, hasher, equaler, track_value, main_alloc>;     
        using allocator     = pool_allocator<main_alloc>;
        using base_type     = allocator;

    private:
        hash_table          m_table;        

        object_table(const object_table&) = delete;
        object_table& operator=(const object_table&) = delete;

        template<class Y>
        VT*                 register_obj(Y&& str);        

    public:
        object_table(size_t capacity = 0);
        ~object_table();

        template<class Y>
        object_ptr          get(const Y& str);
        template<class Y>
        object_ptr          get_existing(const Y& str) const;
        void                unregister_obj(VT* ptr);
        template<class stack_type>
        void                unregister_obj(VT* ptr, stack_type& st);
        size_t              size() const                    { return m_table.size(); };
        size_t              capacity() const                { return m_table.capacity(); };
        double              reuse_stats() const;
        double              collisions() const              { return m_table.collisions();};
        void                close(bool call_destructors);
};

};};