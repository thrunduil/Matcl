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

#include "matcl-scalar/config.h"
#include "matcl-dynamic/details/object.inl"

namespace matcl { namespace details
{

class MATCL_SCALAR_EXPORT matrix_object_interface
{
    private:
        const dynamic::object*  m_object;

    public:
        matrix_object_interface(const dynamic::object* obj);

        bool        displayed_object_matrix(const disp_stream_ptr& os) const;
};

class matrix_object_interface_impl
{
    public:
        virtual bool    displayed_object_matrix(const disp_stream_ptr& os, 
                            const dynamic::object&) const = 0;
};

MATCL_SCALAR_EXPORT 
void set_matrix_object_intrface(matrix_object_interface_impl* oi);

}}
