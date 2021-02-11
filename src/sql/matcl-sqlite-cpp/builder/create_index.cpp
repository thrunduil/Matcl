/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-sqlite-cpp/builder/create_index.h"
#include "sql/matcl-sqlite-cpp/builder/details/create_index_impl.h"

namespace matcl { namespace sql
{

namespace details
{

create_index_template_base::create_index_template_base()
    : m_impl(new details::create_index_impl)
{}

create_index_template_base::create_index_template_base(const impl_ptr& impl)
    : m_impl(impl)
{}

create_index_impl::create_index_impl()
    : m_unique(false), m_if_not_exists(false), m_table("")
{}

std::string create_index_impl::to_str() const
{
    std::string sql = "CREATE ";
    
    if( m_unique )
        sql += "UNIQUE ";
 
    sql += "INDEX ";

    if( m_if_not_exists )
        sql += "IF NOT EXISTS ";

    sql += m_name;
    sql += " ON ";
    sql += m_table.get_name();
    sql += " (";
    sql += m_cols.to_str();
    sql += ")";

    return sql;
}

void create_index_impl::unique(bool b)
{
    m_unique = b;
}

void create_index_impl::if_not_exists(bool b)
{
    m_if_not_exists = b;
}

void create_index_impl::name(const std::string& name)
{
    m_name = name;
}

void create_index_impl::on(const table& table)
{
    m_table = table;
}

void create_index_impl::cols(const column_list& cols)
{
    assert(!cols.empty());
    m_cols = cols;
}

}
}}
