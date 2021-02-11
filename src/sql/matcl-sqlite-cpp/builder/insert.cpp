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

#include "matcl-sqlite-cpp/builder/insert.h"
#include "matcl-sqlite-cpp/builder/Select.h"
#include "sql/matcl-sqlite-cpp/builder/details/insert_impl.h"

namespace matcl { namespace sql
{

namespace details
{

insert_template_base::insert_template_base()
    : m_impl(new details::insert_impl())
{}

insert_template_base::insert_template_base(const impl_ptr& impl)
    : m_impl(impl)
{}

insert_impl::insert_impl()
    : m_cb(conflict_behaviour::default_cb), m_tab("")
{}

void insert_impl::on_conflict(conflict_behaviour cb)
{
    m_cb = cb;
}

void insert_impl::columns(const column_list& cols)
{
    m_cl = cols;
}

void insert_impl::into(const table& t)
{
    m_tab = t;
}

void insert_impl::values(const details::select_template_base& sel)
{
    m_sel.reset(new details::select_template_base(sel));
}

void insert_impl::values(const expression_list& el)
{
    m_vals = el;
}

std::string insert_impl::to_str() const
{
    std::string sql = "INSERT";
    if (conflict_behaviour::default_cb != m_cb )
        sql += " OR ";

    switch( m_cb )
    {
        case conflict_behaviour::rollback:
            sql += "ROLLBACK";
            break;
        case conflict_behaviour::abort:
            sql += "ABORT";
            break;
        case conflict_behaviour::replace:
            sql += "REPLACE";
            break;
        case conflict_behaviour::fail:
            sql += "FAIL";
            break;
        case conflict_behaviour::ignore:
            sql += "IGNORE";
            break;
    }

    sql += " INTO ";
    sql += m_tab.get_name();

    if( !m_cl.empty() )
    {
        sql += " (";
        sql += m_cl.to_str(&column::get_name);
        sql += ")";
    }

    if( m_sel )
    {
        sql += ' ';
        sql += m_sel->to_str();
    }
    else if( !m_vals.empty() )
    {
        sql += " VALUES (";
        sql += m_vals.to_str();
        sql += ")";
    }
    else
    {
        assert( m_cl.empty() ); // assured by external interface
        sql += " DEFAULT VALUES";
    }
    return sql + "\n";
}

}
}}
