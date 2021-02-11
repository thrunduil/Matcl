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

#pragma once

#include "matcl-core/config.h"
#include "matcl-core/memory/alloc.h"
#include <memory>

namespace matcl
{

// an intermediary class implementing cache clearing
class cache_registerer : public matcl_new_delete
{
    public:
        virtual ~cache_registerer(){};

        // release all unused memory
        virtual void    release_cache() = 0;
};

// register on object, that stores memory, which can be released at
// any time; this memory will be released when free_caches function
// is called or malloc fails; only thread local caches can be registered
// in this way; one can also register global caches protected by a mutex
// if they can be cleared at any time from any thread; this function is 
// mostly for internal use; reg is a pointer to class that implements cache
// clearing; if given cache is destroyed; then pointer reg should be explired
MATCL_CORE_EXPORT void   register_cache(const std::weak_ptr<cache_registerer>& reg);

// explicitly release all memory stored in caches
MATCL_CORE_EXPORT void   free_caches();

};
