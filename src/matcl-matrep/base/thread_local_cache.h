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

#include "matcl-matrep/details/fwd_decls.h"

namespace matcl { namespace details
{

//for now this class stores small number of heavily used small objects
//(refcounts and containers) therefore not registered as a cache
template<class tag, int capacity>
class cache
{
    //internal use
    public:
        int             m_pos;
        void*           m_array[capacity];

        static MATCL_THREAD_LOCAL
        cache           g_instance;

    public:
        static cache&   get()   { return g_instance; };

        template<class Allocator>
        void* malloc(Allocator* alloc)
        {
            if (m_pos == 0)
                return alloc->hard_malloc();
            return m_array[--m_pos];
        };

        template<class Allocator>
        void free(void* ptr, Allocator* alloc)
        {
            if (m_pos == capacity)
                return alloc->hard_free(ptr);

            m_array[m_pos++] = ptr;
        };
};

template<class tag>
class cache<tag,0>
{
    public:
        static cache get()  { return cache(); };

        template<class Allocator>
        void* malloc(Allocator* alloc)
        {
            return alloc->hard_malloc();
        };

        template<class Allocator>
        void free(void* ptr, Allocator* alloc)
        {
            return alloc->hard_free(ptr);
        };
};

template<class tag, int capacity>
MATCL_THREAD_LOCAL cache<tag, capacity> cache<tag, capacity>::g_instance;

};};
