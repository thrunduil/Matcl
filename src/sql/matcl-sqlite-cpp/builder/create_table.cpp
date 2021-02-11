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

#include <cassert>
#include "matcl-sqlite-cpp/builder/create_table.h"
#include "sql/matcl-sqlite-cpp/builder/details/create_table_impl.h"
#include "matcl-sqlite-cpp/builder/Select.h"

namespace matcl { namespace sql
{

namespace details
{

create_table_template_base::create_table_template_base()
    : m_impl(new details::create_table_impl())
{}

create_table_template_base::create_table_template_base(const impl_ptr& impl)
    : m_impl(impl)
{}

create_table_impl::create_table_impl()
    : m_temp(false), m_if_not_exists(false), m_table("")
{}

void create_table_impl::temp(bool b)
{
    m_temp = b;
}

void create_table_impl::if_not_exists(bool b)
{
    m_if_not_exists = b;
}

void create_table_impl::tab(const table& t)
{
    m_table = t;
}

void create_table_impl::cols(const column_decl_list& cdl)
{
    m_cdl = cdl;
}

void create_table_impl::data(const details::select_template_base& sel)
{
    m_sel.reset( new details::select_template_base(sel) );
}

void create_table_impl::constraints(const table_constraint_list& tcl)
{
    m_tcl = tcl;
}

std::string create_table_impl::to_str() const
{
    assert( m_cdl.empty() == static_cast<bool>(m_sel) ); // no cols iff select exists
 
    std::string sql = "CREATE";

    if( m_temp )
        sql += " TEMPORARY";

    sql += " TABLE";

    if( m_if_not_exists )
        sql += " IF NOT EXISTS";

    sql += ' ';
    sql += m_table.get_name();
    sql += ' ';

    if( !m_cdl.empty() )
    {
        sql += "( ";
        sql += m_cdl.to_str();
    
        if( !m_tcl.empty() )
        {
            sql += ", ";
            sql += m_tcl.to_str();
        }
        
        sql += ")";   
    }
    else if( m_sel )
    {
        sql += "AS ";
        sql += m_sel->to_str();
    }

    return sql + "\n";
}

}
}}
