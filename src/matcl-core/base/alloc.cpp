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

#include "matcl-core/general/alloc.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/general/exception.h"
#include "matcl-core/general/memory.h"

#include <algorithm>

//TODO: move to config
//#define USE_DLMALLOC

#ifdef USE_DLMALLOC
    #include "matcl-core/base/malloc.h"
#endif

namespace matcl { namespace details
{

#ifdef USE_DLMALLOC
   
    matcl::default_spinlock_mutex m_mutex;

    struct allocator_impl
    {
        static void* malloc(size_t n)
        {
            std::unique_lock<default_spinlock_mutex> lock(m_mutex);
            return ::dlmalloc(n);
        }

        static void* aligned_malloc(size_t n, size_t align)
        {
            std::unique_lock<default_spinlock_mutex> lock(m_mutex);
            return dlmemalign(align, n);
        };

        static void* realloc(void* ptr, size_t n)
        {
            std::unique_lock<default_spinlock_mutex> lock(m_mutex);
            return ::dlrealloc(ptr,n);
        }

        static void* aligned_realloc(void* ptr, size_t old_size, size_t n, size_t align)
        {
            std::unique_lock<default_spinlock_mutex> lock(m_mutex);

            void* ret_ptr = dlrealloc_in_place(ptr, n);

            if (ret_ptr != nullptr)
                return ret_ptr;

            ret_ptr     = dlmemalign(align, n);

            if (!ret_ptr)
                return ret_ptr;

            memcpy(ret_ptr, ptr, std::min(old_size, n));
            return ret_ptr;
        }

        static void free(void* ptr)
        {
            std::unique_lock<default_spinlock_mutex> lock(m_mutex);
            return ::dlfree(ptr);
        }

        static void aligned_free(void* ptr)
        {
            std::unique_lock<default_spinlock_mutex> lock(m_mutex);
            return ::dlfree(ptr);
        }
    };

#else

    struct allocator_impl
    {
        static void* malloc(size_t n)
        {
            return ::malloc(n);
        }

        static void* aligned_malloc(size_t n, size_t align)
        {
            return ::_aligned_malloc(n, align);
        };

        static void* realloc(void* ptr, size_t n)
        {
            return ::realloc(ptr, n);
        }

        static void* aligned_realloc(void* ptr, size_t old_size, size_t n, size_t align)
        {
            (void)old_size;
            return _aligned_realloc(ptr, n, align);
        }

        static void free(void* ptr)
        {
            return ::free(ptr);
        }

        static void aligned_free(void* ptr)
        {
            return _aligned_free(ptr);
        }
    };

#endif

static const Integer MAGIC_NUMBER   = 537597246;

static void set_magic_number(void* ptr0, size_t n)
{
    char* ptr       = (char*)ptr0;
    Integer* last   = (Integer*)(ptr + n);
    *last           = MAGIC_NUMBER;
};

static void check_magic_number(void* ptr0, size_t bytes)
{
    char* ptr       = (char*)ptr0;
    Integer* last   = (Integer*)(ptr + bytes);

    if (*last != MAGIC_NUMBER)
        error::memory_corrupted();
};

void* default_allocator::malloc(size_t n)
{
    void* ptr = allocator_impl::malloc(n + sizeof(Integer));

    if (!ptr)
    {
        matcl::free_caches();
        ptr = allocator_impl::malloc(n + sizeof(Integer));
    }
        
    if (ptr)
        set_magic_number(ptr, n);

    return ptr;
};
void* default_allocator::aligned_malloc(size_t n)
{
    void* ptr   = allocator_impl::aligned_malloc(n + sizeof(Integer), MATCL_CACHE_LINE_SIZE);

    if (!ptr)
    {
        matcl::free_caches();
        ptr     = allocator_impl::aligned_malloc(n + sizeof(Integer), MATCL_CACHE_LINE_SIZE);
    }

    if (ptr)
        set_magic_number(ptr, n);

    return ptr;
};

void* default_allocator::realloc(void* ptr, size_t old_size, size_t n)
{
    if (ptr)
        check_magic_number(ptr, old_size);

    if (n > 0)
    {
        void* new_ptr = allocator_impl::realloc(ptr,n + sizeof(Integer));

        if (!new_ptr)
        {
            matcl::free_caches();
            new_ptr     = allocator_impl::realloc(ptr,n + sizeof(Integer));
        }

        if (new_ptr)
            set_magic_number(new_ptr, n);

        return new_ptr;
    }
    else
    {
        void* new_ptr = allocator_impl::realloc(ptr,n);
        return new_ptr;
    };
}
void* default_allocator::aligned_realloc(void* ptr, size_t old_size, size_t n)
{
    if (ptr)
        check_magic_number(ptr, old_size);

    old_size    += sizeof(Integer);

    if (n > 0)
    {
        void* new_ptr = allocator_impl::aligned_realloc(ptr, old_size,
                                 n + sizeof(Integer), MATCL_CACHE_LINE_SIZE);

        if (!new_ptr)
        {
            matcl::free_caches();
            new_ptr     = allocator_impl::aligned_realloc(ptr, old_size, 
                                n + sizeof(Integer), MATCL_CACHE_LINE_SIZE);
        }

        if (new_ptr)
            set_magic_number(new_ptr, n);

        return new_ptr;
    }
    else
    {
        return allocator_impl::aligned_realloc(ptr, old_size, n, 
                                            MATCL_CACHE_LINE_SIZE);
    };
}

void default_allocator::free(void* ptr, size_t bytes)
{
    if (ptr)
        check_magic_number(ptr, bytes);

    return allocator_impl::free(ptr);
}

void default_allocator::aligned_free(void* ptr, size_t bytes)
{
    if (ptr)
        check_magic_number(ptr, bytes);

    return allocator_impl::aligned_free(ptr);
}

};};