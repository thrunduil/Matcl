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

#include "matcl-core/details/object_interface.h"

namespace matcl { namespace details
{

object_interface_impl* g_object_impl = nullptr;

void details::set_object_intrface(object_interface_impl* oi)
{
    g_object_impl = oi;
};

const_object_interface::const_object_interface(const dynamic::object* obj)
    :m_object(obj)
{};

bool const_object_interface::is_zero() const
{
    return g_object_impl->is_zero(*m_object);
};

void const_object_interface::disp_elem(printer& p, Integer w, align_type at, 
                Integer value_pos) const
{
    g_object_impl->disp_elem(*m_object, p, w, at, value_pos);
};

void const_object_interface::write(std::ostream& os)
{
    g_object_impl->write(*m_object, os);
};

nonconst_object_interface::nonconst_object_interface(dynamic::object* obj)
    :m_object(obj)
{};

void nonconst_object_interface::read(std::istream& is)
{
    g_object_impl->read(*m_object, is);
};

const_type_interface::const_type_interface(const dynamic::Type* ty)
    :m_type(ty)
{};

void const_type_interface::write(std::ostream& os)
{
    g_object_impl->write_type(*m_type, os);
};

nonconst_type_interface::nonconst_type_interface(dynamic::Type* ty)
    :m_type(ty)
{};

void nonconst_type_interface::read(std::istream& is)
{
    g_object_impl->read_type(*m_type, is);
}

}}
