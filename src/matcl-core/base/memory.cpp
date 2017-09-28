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

#include "matcl-core/general/memory.h"
#include "matcl-core/general/thread.h"
#include "matcl-dynamic/utils.h"

#include <vector>

namespace matcl { namespace details
{

//cache management
class cache_handler
{
    private:
        using mutex_type        = matcl::default_mutex;
        using registerer_ptr    = std::weak_ptr<cache_registerer>;
        using vec_registerers   = std::vector<registerer_ptr>;

    private:
        mutex_type              m_mutex;
        vec_registerers         m_caches;

    public:
        static void             register_cache(const registerer_ptr& reg);
        static void             free_caches();

    private:
        cache_handler();
        static cache_handler*   get();

        void                    implement_free_caches();
};

cache_handler::cache_handler()
{
    //explicitly register caches from matcl::dynamic
    //registerer_ptr reg_dyn = matcl::dynamic::get_cache_registerer();
    //m_caches.push_back(reg_dyn);
}

cache_handler* cache_handler::get()
{
    static cache_handler* inst = new cache_handler();
    return inst;
};

void cache_handler::register_cache(const std::weak_ptr<cache_registerer>& reg)
{
    cache_handler* ch = get();

    using lock_type = std::unique_lock<mutex_type>;
    lock_type lock(ch->m_mutex);

    ch->m_caches.push_back(reg);
};

void cache_handler::free_caches()
{
    cache_handler* ch = get();

    using lock_type = std::unique_lock<mutex_type>;
    lock_type lock(ch->m_mutex);

    ch->implement_free_caches();
};

void cache_handler::implement_free_caches()
{
    for (auto& pos : m_caches)
    {
        auto reg    = pos.lock();

        if (!reg)
            continue;

        reg->release_cache();
    };
};

}};

void matcl::register_cache(const std::weak_ptr<cache_registerer>& reg)
{
    details::cache_handler::register_cache(reg);
};

// explicitly release all memory stored in caches
void matcl::free_caches()
{
    details::cache_handler::free_caches();
};
