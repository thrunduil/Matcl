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

#include "matcl-sqlite-cpp/builder/column_decl.h"
#include "sql/matcl-sqlite-cpp/builder/details/column_decl_impl.h"
#include "sql/matcl-sqlite-cpp/builder/details/conflict_impl.h"
#include "sql/matcl-sqlite-cpp/builder/details/expression_impl.h"

namespace matcl { namespace sql
{

column_decl::column_decl(const std::string& name)
    : m_impl(new details::column_decl_impl(name))
{}

std::string column_decl::to_str() const
{
    return m_impl->to_str();
}

column_decl& column_decl::type(column_type ct)
{
    m_impl->type(ct);
    return *this;
}

column_decl& column_decl::primary_key(bool ascending, conflict_behaviour cb,
                                    bool autoincrement)
{
    m_impl->primary_key(ascending, cb, autoincrement);
    return *this;
}

column_decl& column_decl::not_null(conflict_behaviour cb)
{
    m_impl->not_null(cb);
    return *this;
}

column_decl& column_decl::unique(conflict_behaviour cb)
{
    m_impl->unique(cb);
    return *this;
}

column_decl& column_decl::check(const expression& e)
{
    m_impl->check(e);
    return *this;
}

column_decl& column_decl::default_val(const expression& e)
{
    m_impl->default_val(e);
    return *this;
}

column_decl& column_decl::collate(collation c)
{
    m_impl->collate(c);
    return *this;
}

column_decl& column_decl::foreign_key(const table& t, const column_list& cl)
{
    m_impl->foreign_key(t, cl);
    return *this;
}

namespace details
{

column_decl_impl::column_decl_impl(const std::string& name)
    : m_name(name), m_type(column_type::no_type)
{}

std::string column_decl_impl::to_str() const
{
    std::string sql = m_name;
 
    switch(m_type)
    {
        case column_type::integer:
            sql += " INTEGER";
            break;
        case column_type::text:
            sql += " TEXT";
            break;
        case column_type::blob:
            sql += " BLOB";
            break;
        case column_type::numeric:
            sql += " NUMERIC";
            break;
        case column_type::real:
            sql += " REAL";
            break;
    }

    if( m_pkey )
    {
        sql += " PRIMARY KEY";

        if(m_pkey->m_ascending)
            sql += " ASC";
        else
            sql += " DESC";

        sql += cb_to_str(m_pkey->m_cb);

        if(m_pkey->m_autoincrement)
            sql += " AUTOINCREMENT";
    }

    if( m_not_null )
    {
        sql += " NOT NULL";
        sql += cb_to_str(*m_not_null);
    }

    if( m_unique )
    {
        sql += " UNIQUE";
        sql += cb_to_str(*m_unique);
    }

    if( m_check.exists() )
    {
        sql += " CHECK (";
        sql += m_check.to_str();
        sql += ")";
    }

    if( m_default_val.exists() )
    {
        sql += "DEFAULT (";
        sql += m_default_val.to_str();
        sql += ")";
    }

    if( m_collate )
    {
        sql += collation_to_str(*m_collate);
    }

    if( m_fkey )
    {
        sql += " REFERENCES ";
        sql += m_fkey->m_tab.get_name();
        if( !m_fkey->m_cols.empty() )
        {
            sql += " (";
            sql += m_fkey->m_cols.to_str();
            sql += ')';
        }
    }
    return sql;
}

void column_decl_impl::type(column_type ct)
{
    m_type = ct;
}

void column_decl_impl::primary_key(bool ascending, conflict_behaviour cb,
                                bool autoincrement)
{
    m_pkey.reset( new primary_key_info(ascending, cb, autoincrement) );
}

void column_decl_impl::not_null(conflict_behaviour cb)
{
    m_not_null.reset( new conflict_behaviour(cb) );
}

void column_decl_impl::unique(conflict_behaviour cb)
{
    m_unique.reset( new conflict_behaviour (cb) );
}

void column_decl_impl::check(const expression& e)
{
    m_check = e;
}

void column_decl_impl::default_val(const expression& e)
{
    m_default_val = e;
}

void column_decl_impl::collate(collation c)
{
    m_collate.reset( new collation (c) );
}

void column_decl_impl::foreign_key(const table& t, const column_list& cl)
{
    m_fkey.reset( new foreign_key_info(t,cl) );
}

}
}}

