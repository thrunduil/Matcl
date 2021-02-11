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

#include "matcl-matrep/general/config.h"
#include "matcl-core/general/thread.h"

#pragma warning(push)
#pragma warning(disable:4251)    //needs to have dll-interface 

namespace matcl { namespace details
{

using atomic_long = matcl::atomic<long>;

template<class refcount_pool>
struct refcount_str
{
    public:
        refcount_str(long val = 0);

        long                    get_count() const;
        refcount_str*           increase();
        bool                    decrease();
        bool                    is_unique() const;
        void                    destroy();

        static refcount_str*    create(long val);

    protected:
        atomic_long             m_refcount;

    private:
        refcount_str(const refcount_str&)            = delete;
        refcount_str& operator=(const refcount_str&) = delete;
};

struct default_item_counter
{
    long    m_count;

    void    reset()     { m_count = 0; };
    void    increase()  { ++m_count; };
    void    decrease()  {--m_count; };
    long    count()     { return m_count; };
};

struct empty_item_counter
{
    void    reset()     {};
    void    increase()  {};
    void    decrease()  {};
    long    count()     { return 0; };
};

struct default_refcount_pool
{
    static const int cache_capacity = 50;

  #if (DEBUG_MEMORY == 1)
    using item_counter = default_item_counter;
  #else
    using item_counter = empty_item_counter;
  #endif
};

struct scatter_refcount_pool
{
    static const int cache_capacity = 10;
    using item_counter = empty_item_counter;
};

using refcount_str_default = refcount_str<default_refcount_pool>;
using refcount_str_scatter = refcount_str<scatter_refcount_pool>;

//not thread safe
long MATCL_MATREP_EXPORT no_existing_objects();

//not thread safe
void MATCL_MATREP_EXPORT close_refcount();

};};

#include "matcl-matrep/details/refcount.inl"

#pragma warning(pop)