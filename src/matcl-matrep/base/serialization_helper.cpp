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

#include "matcl-matrep/IO/serialization_helper.h"
#include "matcl-core/general/thread.h"

#include <map>

namespace matcl { namespace details
{

class helper_map
{
    private:
        using ptr_type  = serialization_helper_base;
        using key_type  = std::pair<std::string, std::string>;
        using map_type  = std::map<key_type, ptr_type*>;
        
    private:        
        map_type    m_map;

    public:
        ptr_type*   get_registered(const std::string& unique_id, 
                                const std::string& base_class_name);
        void        make_register(const std::string& unique_id, 
                                const std::string& base_class_name, ptr_type* ptr);
};

helper_map::ptr_type* helper_map::get_registered(const std::string& unique_id, 
                        const std::string& base_class_name)
{
    auto pos = m_map.find(key_type(unique_id,base_class_name));

    if (pos == m_map.end())
        return nullptr;
    else
        return pos->second;
};

void helper_map::make_register(const std::string& unique_id, 
                        const std::string& base_class_name, ptr_type* ptr)
{
    using value_type = map_type::value_type;

    m_map.insert(value_type(key_type(unique_id,base_class_name), ptr));
};

using spin_lock = default_spinlock_mutex;
using lock_type = std::unique_lock<spin_lock>;

static helper_map   g_helper_map;
static spin_lock    g_mutex;

serialization_helper_base* 
register_serialization_impl::get_registered(const std::string& unique_id, 
                        const char* base_class_name0)
{
    lock_type lock(g_mutex);

    std::string base_class_name = base_class_name0? std::string(base_class_name0) : "";
    return g_helper_map.get_registered(unique_id,base_class_name);
};

void register_serialization_impl::make_register(const std::string& unique_id, 
                        const char* base_class_name0, ptr_type* ptr)
{
    lock_type lock(g_mutex);

    std::string base_class_name = base_class_name0? std::string(base_class_name0) : "";
    return g_helper_map.make_register(unique_id,base_class_name,ptr);
};

};};
