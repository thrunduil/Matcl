/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include <cassert>
#include "matcl-sqlite-cpp/builder/column.h"
#include "sql/matcl-sqlite-cpp/builder/details/column_impl.h"

namespace matcl { namespace sql
{

column::column(const std::string& name)
    : m_impl(new details::column_impl(name))
{}

std::string column::to_str() const
{
    return m_impl->to_str();
}

void column::set_table(const table& t)
{
    m_impl->set_table(t);
}

std::string column::get_name() const
{
    return m_impl->get_name();
}

namespace details
{

column_impl::column_impl(const std::string& name)
    : m_name(name), m_table("")
{
    assert(name.find(' ') == std::string::npos);
}

std::string column_impl::to_str() const
{
    std::string ref = m_table.get_ref();

    if( "" != ref )
        return ref + '.' + m_name;

    return m_name;
}

void column_impl::set_table(const table& t)
{
    m_table = t;
}

std::string column_impl::get_name() const
{
    return m_name;
}

}

}}
