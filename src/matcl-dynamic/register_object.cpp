/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017
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

#include "matcl-dynamic/details/register_object.h"
#include "matcl-core/config.h"
#include "type_table.h"

namespace matcl { namespace dynamic { namespace details
{

Type register_object_helper::get_object(const char* obj_name)
{
    return type_table::get()->get_type(obj_name);
};

std::string	register_object_helper::get_name(const char* name)
{
    return type_table::get_class_name(name);
};

register_object_helper::register_object_helper(const char* obj_name, Type (*creator)())
{
    type_table::get()->register_object(obj_name,creator);
};

};};};
