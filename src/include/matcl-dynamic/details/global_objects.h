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

#include "matcl-dynamic/config.h"
#include "matcl-dynamic/details/type_object.h"
#include "matcl-dynamic/details/object_data_pool.h"
#include "matcl-core/memory/alloc.h"

namespace matcl { namespace dynamic { namespace details
{

class global_objects
{
    private:
        using pool_type     = object_data_pool_impl;

    public:
        template<class T>
        static type_impl*   initialize_type(pool_type*& pool);
};

template<class T>
type_impl* global_objects::initialize_type(pool_type*& pool)
{
    type_impl* type = new details::type_object<T>();
    pool            = new object_data_pool_impl(sizeof(object_data<T>));

    return type;
};

};};};