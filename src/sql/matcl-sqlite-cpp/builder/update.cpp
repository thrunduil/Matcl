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

#include "matcl-sqlite-cpp/builder/update.h"
#include "sql/matcl-sqlite-cpp/builder/details/update_impl.h"

namespace matcl { namespace sql
{

namespace details
{

update_template_base::update_template_base()
    : m_impl(new details::update_impl)
{}

update_template_base::update_template_base(const impl_ptr& impl)
    : m_impl(impl)
{}

update_impl::update_impl()
    : m_cb(conflict_behaviour::default_cb), m_tab("")
{
}

void update_impl::on_conflict(conflict_behaviour cb)
{
    m_cb = cb;
}

void update_impl::tab(const table& t)
{
    m_tab = t;
}

void update_impl::set(const assignment_list& al)
{
    m_al = al;
}

void update_impl::where(const expression& e)
{
    m_cond = e;
}

std::string update_impl::to_str() const
{
    std::string sql = "UPDATE";

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

    sql += " ";
    sql += m_tab.get_name();
    sql += " SET ";
    sql += m_al.to_str();

    if( m_cond.exists() )
    {
        sql += " WHERE ";
        sql += m_cond.to_str();
    }

    return sql + "\n";
}

}
}}
