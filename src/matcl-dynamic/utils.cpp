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

#include "matcl-dynamic/utils.h"
#include "matcl-core/memory/memory.h"
#include "type_table.h"

namespace matcl { namespace dynamic
{

// registerer of caches in matcl::dynamic
class dyn_chache_registerer : public matcl::cache_registerer
{
    virtual void release_cache() override
    {
        free_cache();
    }
};

// explicitly release all memory stored in caches
void dynamic::free_cache()
{
    details::type_table::get()->free_cache();
}

struct register_cache_dynamic
{
    std::shared_ptr<cache_registerer> reg;

    register_cache_dynamic()
    {
        reg = std::shared_ptr<cache_registerer>(new dyn_chache_registerer());

        matcl::register_cache(reg);
    };
};

static register_cache_dynamic g_reg_cache;

};};

