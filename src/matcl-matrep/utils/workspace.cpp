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

#include "matcl-matrep/utils/workspace.h"
#include "matcl-core/general/thread.h"
#include "matcl-core/error/exception_classes.h"
#include "../base/alloc.h"
#include "matcl-core/memory/memory.h"

namespace matcl {namespace details
{

//helper class to register workspace_provider as cache
struct workspace_cache : cache_registerer
{
    virtual void        release_cache() override;
};

//simple implementation, all stored workspaces have the same size WORK_SIZE
class workspace_provider
{
    private:
        using reg_ptr   = std::shared_ptr<workspace_cache>;
        using allocator = default_allocator<true>;

    private:
        static const Integer MAX_ARRAYS = 20;
        static const size_t  WORK_SIZE  = 1024;

    private:
        struct item
        {
            void*   m_ptr;
            size_t  m_sizeof;
        };

    private:
        item        m_arrays[MAX_ARRAYS];
        Integer     m_pos;
        reg_ptr     m_registerer;

    public:
        workspace_provider();
        ~workspace_provider();

        static workspace_provider&  get();

        void*       get_workspace_ptr(size_t& size_of);
        void        release_workspace_ptr(void* ptr, size_t size_of);
        void*       resize_workspace_ptr(void* old, size_t& size_of, size_t new_sizeof);

        //release unused memory
        void        release_cache();

    private:
        workspace_provider(const workspace_provider&) = delete;
        workspace_provider& operator=(const workspace_provider&) = delete;
};

workspace_provider::workspace_provider()
    :m_pos(0)
{
    using atomic_int   = matcl::atomic<int>;

    static atomic_int registered = 0;

    if (registered.fetch_add(1) == 0)
    {
        //avoid registering this cache many times
        m_registerer = reg_ptr(new workspace_cache());
        register_cache(m_registerer);
    };
};

workspace_provider::~workspace_provider()
{
    for (Integer i = 0; i < m_pos; ++i)
        allocator::aligned_free(m_arrays[i].m_ptr, m_arrays[i].m_sizeof);
    m_pos = 0;
};

void workspace_provider::release_cache()
{
    for (Integer i = 0; i < m_pos; ++i)
        allocator::aligned_free(m_arrays[i].m_ptr, m_arrays[i].m_sizeof);
};


MATCL_THREAD_LOCAL
workspace_provider* g_workspace_prov = nullptr;

void workspace_cache::release_cache()
{
    //release only unused memory in this thread
    if (g_workspace_prov)
        g_workspace_prov->release_cache();
};

workspace_provider& workspace_provider::get()
{
    if (!g_workspace_prov)
        g_workspace_prov = new workspace_provider();

    return *g_workspace_prov;
};

void* workspace_provider::resize_workspace_ptr(void* old, size_t& size_of, size_t new_sizeof)
{
    if (new_sizeof <= size_of)
        return old;

    size_of         = new_sizeof;
    void* new_ptr   = allocator::aligned_realloc(old, size_of, new_sizeof);
    return new_ptr;    
}

void  workspace_provider::release_workspace_ptr(void* ptr, size_t size_of)
{
    if (m_pos == MAX_ARRAYS)
        allocator::aligned_free(ptr, size_of);
    else if (size_of != WORK_SIZE)
        allocator::aligned_free(ptr, size_of);
    else
        m_arrays[m_pos++] = item{ptr, size_of};
};

void* workspace_provider::get_workspace_ptr(size_t& size_of)
{
    if (size_of > WORK_SIZE)
    {
        void* ptr   = allocator::aligned_malloc(size_of);
        return ptr;
    }
    if (m_pos == 0)
    {
        size_of     = WORK_SIZE;
        void* ptr   = allocator::aligned_malloc(size_of);
        return ptr;
    }

    size_of     = WORK_SIZE;
    void* ptr   = m_arrays[--m_pos].m_ptr;
    return ptr;
};

workspace_base::workspace_base()
    :m_size_of(0), m_total_size_of(0)
{
    set_magic_number(m_data.m_data, 0);
};

workspace_base::workspace_base(size_t size_of)
    :m_size_of(size_of), m_total_size_of(0)
{
    if (size_of > MAX_ON_STACK)
    {
        m_total_size_of = size_of + sizeof(Integer);
        m_data.m_ptr    = workspace_provider::get().get_workspace_ptr(m_total_size_of);
        set_magic_number(m_data.m_ptr, size_of);
    }
    else
    {
        set_magic_number(m_data.m_data, size_of);
    }
}

workspace_base::workspace_base(workspace_base&& other)
    : m_size_of(other.m_size_of), m_total_size_of(other.m_total_size_of)
{
    if (m_size_of > MAX_ON_STACK)
        m_data.m_ptr = other.m_data.m_ptr;
    else
        memcpy(m_data.m_data, other.m_data.m_data, MAX_ON_STACK + sizeof(Integer));
};

workspace_base::~workspace_base()
{
    if (m_size_of > MAX_ON_STACK)
    {
        check_magic_number(m_data.m_ptr, m_size_of);
        workspace_provider::get().release_workspace_ptr(m_data.m_ptr, m_total_size_of);
    }
    else
    {
        check_magic_number(m_data.m_data, m_size_of);
    }
}

size_t workspace_base::size_of() const
{
    return m_size_of;
}

void* workspace_base::get_ptr()
{
    if (m_size_of > MAX_ON_STACK)
        return m_data.m_ptr;
    else
        return m_data.m_data;

}

static const Integer MAGIC_NUMBER = 918273645;

void workspace_base::set_magic_number(void* data, size_t size_of)
{
    Integer* ptr    = (Integer*)((char*)data + size_of);
    *ptr            = MAGIC_NUMBER;
}

void workspace_base::check_magic_number(void* data, size_t size_of)
{
    Integer* ptr    = (Integer*)((char*)data + size_of);

    if (*ptr != MAGIC_NUMBER)
        error::memory_corrupted();
};

void workspace_base::resize(size_t new_size_of)
{
    if (this->m_total_size_of != 0)
        check_magic_number(m_data.m_ptr, m_size_of);
    else
        check_magic_number(m_data.m_data, m_size_of);

    m_size_of   = new_size_of;

    if (new_size_of > MAX_ON_STACK)
    {
        if (this->m_total_size_of > 0)
        {
            m_data.m_ptr    = workspace_provider::get().resize_workspace_ptr(m_data.m_ptr, m_total_size_of, 
                                                                             new_size_of + sizeof(Integer));            
        }
        else
        {
            m_total_size_of = new_size_of + sizeof(Integer);
            m_data.m_ptr    = workspace_provider::get().get_workspace_ptr(m_total_size_of);
        };        

        set_magic_number(m_data.m_ptr, m_size_of);
    }
    else
    {
        if (this->m_total_size_of > 0)
        {
            workspace_provider::get().release_workspace_ptr(m_data.m_ptr, m_total_size_of);
            m_total_size_of = 0;
        }

        set_magic_number(m_data.m_data, m_size_of);
    }
}

};};
