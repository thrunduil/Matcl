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

#include "matcl-core/memory/alloc.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/error/exception_classes.h"
#include "matcl-core/memory/memory.h"

#if MATCL_DEBUG_MEMORY
    #include "matcl-core/details/leak_detector.h"
#endif

#include <algorithm>

#if MATCL_USE_DLMALLOC
    #include "matcl-core/base/malloc.h"
#endif

#pragma warning(push)
#pragma warning(disable : 4127) // conditional expression is constant

namespace matcl { namespace details
{

#if MATCL_USE_DLMALLOC
   
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

template<bool Throw>
struct check_ptr;

template<>
struct check_ptr<false>
{
    static void eval(void*, size_t)
    {};
};

template<>
struct check_ptr<true>
{
    static void eval(void* ptr, size_t size)
    {
        if (!ptr && size > 0)
            throw error::alloc(size);
    }
};

}}

namespace matcl
{

template<bool Throw_bad_alloc>
void* default_allocator<Throw_bad_alloc>::malloc(size_t n)
{
    void* ptr = md::allocator_impl::malloc(n + sizeof(Integer));

    if (!ptr)
    {
        matcl::free_caches();
        ptr = md::allocator_impl::malloc(n + sizeof(Integer));
    }
        
    if (ptr)
    {
        md::set_magic_number(ptr, n);

        #if MATCL_DEBUG_MEMORY
            details::leak_detector::report_malloc(ptr);
        #endif
    }

    md::check_ptr<Throw_bad_alloc>::eval(ptr, n);
    return ptr;
};

template<bool Throw_bad_alloc>
void* default_allocator<Throw_bad_alloc>::aligned_malloc(size_t n)
{
    void* ptr   = md::allocator_impl::aligned_malloc(n + sizeof(Integer), MATCL_CACHE_LINE_SIZE);

    if (!ptr)
    {
        matcl::free_caches();
        ptr     = md::allocator_impl::aligned_malloc(n + sizeof(Integer), MATCL_CACHE_LINE_SIZE);
    }

    if (ptr)
    {
        md::set_magic_number(ptr, n);

        #if MATCL_DEBUG_MEMORY
            details::leak_detector::report_malloc(ptr);
        #endif
    }

    md::check_ptr<Throw_bad_alloc>::eval(ptr, n);
    return ptr;
};

template<bool Throw_bad_alloc>
void* default_allocator<Throw_bad_alloc>::simple_malloc(size_t n)
{
    void* ptr = md::allocator_impl::malloc(n);

    if (!ptr)
    {
        matcl::free_caches();
        ptr = md::allocator_impl::malloc(n);
    }
    
    #if MATCL_DEBUG_MEMORY
        details::leak_detector::report_malloc(ptr);
    #endif

    md::check_ptr<Throw_bad_alloc>::eval(ptr, n);
    return ptr;
};

template<bool Throw_bad_alloc>
void* default_allocator<Throw_bad_alloc>
            ::realloc(void* ptr, size_t old_size, size_t n)
{
    if (ptr)
        md::check_magic_number(ptr, old_size);

    void* new_ptr;

    if (n > 0)
    {
        new_ptr = md::allocator_impl::realloc(ptr, n + sizeof(Integer));

        if (!new_ptr)
        {
            matcl::free_caches();
            new_ptr     = md::allocator_impl::realloc(ptr,n + sizeof(Integer));
        }

        if (new_ptr)
            md::set_magic_number(new_ptr, n);
    }
    else
    {
        new_ptr = md::allocator_impl::realloc(ptr,n);        
    };

    #if MATCL_DEBUG_MEMORY
        if (new_ptr != ptr)
        {
            if (ptr)
                details::leak_detector::report_free(ptr);

            if (new_ptr)
                details::leak_detector::report_malloc(new_ptr);
        }
    #endif

    md::check_ptr<Throw_bad_alloc>::eval(new_ptr, n);
    return new_ptr;
}

template<bool Throw_bad_alloc>
void* default_allocator<Throw_bad_alloc>
                ::simple_realloc(void* ptr, size_t n)
{
    void* new_ptr = md::allocator_impl::realloc(ptr,n);

    #if MATCL_DEBUG_MEMORY
        if (new_ptr != ptr)
        {
            if (ptr)
                details::leak_detector::report_free(ptr);

            if (new_ptr)
                details::leak_detector::report_malloc(new_ptr);
        }
    #endif

    md::check_ptr<Throw_bad_alloc>::eval(new_ptr, n);
    return new_ptr;
}

template<bool Throw_bad_alloc>
void* default_allocator<Throw_bad_alloc>
                ::aligned_realloc(void* ptr, size_t old_size, size_t n)
{
    if (ptr)
        md::check_magic_number(ptr, old_size);

    old_size    += sizeof(Integer);

    void* new_ptr;

    if (n > 0)
    {
        new_ptr = md::allocator_impl::aligned_realloc(ptr, old_size,
                                 n + sizeof(Integer), MATCL_CACHE_LINE_SIZE);

        if (!new_ptr)
        {
            matcl::free_caches();
            new_ptr     = md::allocator_impl::aligned_realloc(ptr, old_size, 
                                n + sizeof(Integer), MATCL_CACHE_LINE_SIZE);
        }

        if (new_ptr)
            md::set_magic_number(new_ptr, n);
    }
    else
    {
        new_ptr = md::allocator_impl::aligned_realloc(ptr, old_size, n, 
                                            MATCL_CACHE_LINE_SIZE);
    };

    #if MATCL_DEBUG_MEMORY
        if (new_ptr != ptr)
        {
            if (ptr)
                details::leak_detector::report_free(ptr);

            if (new_ptr)
                details::leak_detector::report_malloc(new_ptr);
        }
    #endif

    md::check_ptr<Throw_bad_alloc>::eval(new_ptr, n);
    return new_ptr;
}

template<bool Throw_bad_alloc>
void 
default_allocator<Throw_bad_alloc>::free(void* ptr, size_t bytes)
{
    if (ptr)
    {
        md::check_magic_number(ptr, bytes);

        #if MATCL_DEBUG_MEMORY
            details::leak_detector::report_free(ptr);
        #endif
    }

    return md::allocator_impl::free(ptr);
}

template<bool Throw_bad_alloc>
void default_allocator<Throw_bad_alloc>::simple_free(void* ptr)
{
    #if MATCL_DEBUG_MEMORY
        if (ptr)
            details::leak_detector::report_free(ptr);
    #endif

    return md::allocator_impl::free(ptr);
}

template<bool Throw_bad_alloc>
void default_allocator<Throw_bad_alloc>::aligned_free(void* ptr, size_t bytes)
{
    if (ptr)
    {
        md::check_magic_number(ptr, bytes);

        #if MATCL_DEBUG_MEMORY
            details::leak_detector::report_free(ptr);
        #endif
    };

    return md::allocator_impl::aligned_free(ptr);
}

template MATCL_CORE_EXPORT struct default_allocator<false>;
template MATCL_CORE_EXPORT struct default_allocator<true>;

}

#pragma warning(pop)