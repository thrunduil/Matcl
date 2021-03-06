/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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

#include "matcl-mp/cache.h"
#include "utils/impl_types.h"
#include "matcl-core/error/exception_classes.h"

namespace matcl
{

unique_id::unique_id(const std::string& optional_name)
    :m_name(optional_name)
{
    using atomic_int    = matcl::atomic<int>;

    static atomic_int global_code = 0;

    this->m_code        = ++global_code;
};

bool unique_id::operator<(const unique_id& other) const
{ 
    return this->m_code < other.m_code; 
}

bool unique_id::operator>(const unique_id& other) const
{ 
    return this->m_code > other.m_code; 
}

bool unique_id::operator<=(const unique_id& other) const
{ 
    return this->m_code <= other.m_code; 
}

bool unique_id::operator>=(const unique_id& other) const
{ 
    return this->m_code >= other.m_code; 
}

bool unique_id::operator==(const unique_id& other) const
{
    return this->m_code == other.m_code; 
}

bool unique_id::operator!=(const unique_id& other) const
{ 
    return this->m_code != other.m_code; 
}

const std::string& unique_id::to_string() const
{ 
    return m_name; 
};

cache* cache::get()
{
    MATCL_THREAD_LOCAL
    static cache* c = new cache();
    return c;
};

bool cache::has(unique_id id, precision prec)
{
    const cache* c  = cache::get();

    auto pos = c->m_map_float.find(id);
    if (pos == c->m_map_float.end())
        return false;

    precision save_prec = pos->second.first;
    if (prec <= save_prec)
        return true;
    else
        return false;
}

#pragma warning(push)
#pragma warning(disable : 4702) //unreachable code

mp_float cache::get(unique_id id, precision prec)
{
    const cache* c  = cache::get();

    auto pos = c->m_map_float.find(id);

    if (pos == c->m_map_float.end())
        throw error::value_not_in_cache(id.to_string(), (Integer)prec, -1);

    precision save_prec = pos->second.first;

    if (prec > save_prec)
    {
        throw error::value_not_in_cache(id.to_string(), (Integer)prec, 
                                (Integer)save_prec);
    }

    return pos->second.second;
};

#pragma warning(pop)

void cache::set(unique_id id, precision prec, const mp_float& val)
{
    cache* c    = cache::get();
    auto pos    = c->m_map_float.find(id);

    using map_value_type = map_float::value_type;

    if (pos == c->m_map_float.end())
    {
        c->m_map_float.insert(pos, map_value_type(id, value_float(prec, val)));
        return;
    }

    precision save_prec = pos->second.first;

    //do not update values with higher precision
    if (prec < save_prec)
        return;
    
    c->m_map_float.insert(pos, map_value_type(id, value_float(prec, val)));
    return;
};

void cache::remove(unique_id id)
{
    cache* c    = cache::get();

    c->m_map_float.erase(id);
};

void cache::clear()
{
    cache* c    = cache::get();

    c->m_map_float.clear();
    mpfr_free_cache();
};

void cache::clear_global()
{
    clear();
}

void cache::close_global()
{
    delete this;
}

};