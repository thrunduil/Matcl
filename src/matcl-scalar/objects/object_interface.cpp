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

#include "matcl-scalar/details/object_interface.h"

namespace matcl { namespace details
{

struct matrix_object_interface_default : matrix_object_interface_impl
{
    bool displayed_object_matrix(const disp_stream_ptr& os, 
                        const dynamic::object&) const override
    {
        (void)os;
        return false;
    };
};

matrix_object_interface_impl* g_object_impl = nullptr;

void details::set_matrix_object_intrface(matrix_object_interface_impl* oi)
{
    g_object_impl = oi;
};

static matrix_object_interface_impl* get_object_impl()
{
    if (!g_object_impl)
    {
        static matrix_object_interface_default g_default_matrix_interface;
        return &g_default_matrix_interface;
    }
    else
        return g_object_impl;
};

matrix_object_interface::matrix_object_interface(const dynamic::object* obj)
    :m_object(obj)
{};

bool matrix_object_interface::displayed_object_matrix(const disp_stream_ptr& os) const
{
    return get_object_impl()->displayed_object_matrix(os, *m_object);
};

}}
