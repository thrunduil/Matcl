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

#include "matcl-sqlite-cpp/builder/delete.h"
#include "sql/matcl-sqlite-cpp/builder/details/delete_impl.h"

namespace matcl { namespace sql
{

namespace details
{

delete_template_base::delete_template_base()
    : m_impl(new details::delete_impl())
{}

delete_template_base::delete_template_base(const impl_ptr& impl)
    : m_impl(impl)
{}

delete_impl::delete_impl()
    : m_tab("")
{}

void delete_impl::tab(const table& t)
{
    m_tab = t;
}

void delete_impl::where(const expression& e)
{
    m_cond = e;
}

std::string delete_impl::to_str() const
{
    std::string sql = "DELETE FROM ";
    sql += m_tab.get_name();
 
    if( m_cond.exists() )
    {
        sql += " WHERE ";
        sql += m_cond.to_str();
    }
    
    return sql + "\n";
}

}
}}
