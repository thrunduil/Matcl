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

#include "matcl-dynamic/details/register_object.inl"
#include "matcl-core/config.h"
#include "type_table.h"

namespace matcl { namespace dynamic { namespace details
{

void register_object_helper::get_type_data(const char* obj_name, Type& ty,
                                           pool_type*& pool)
{
    ty = type_table::get()->get_type_data(obj_name, pool);
};

std::string	register_object_helper::get_name(const type_info& ti)
{
    return type_table::get_class_name(ti.name());
};

register_object_helper::register_object_helper(const type_info& ti, constr_type constr)
{
    const char* obj_name = ti.name();

    type_table::open();
    type_table::get()->register_object(obj_name, constr);
    get_type_data(obj_name, m_type, m_pool);
};

};};};
