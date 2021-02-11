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

#include "matcl-sqlite-cpp/builder/table_constraint.h"
#include "sql/matcl-sqlite-cpp/builder/details/table_constraint_impl.h"
#include "sql/matcl-sqlite-cpp/builder/details/conflict_impl.h"

namespace matcl { namespace sql
{

table_constraint::table_constraint()
    : m_impl(new details::table_constraint_impl)
{}

std::string table_constraint::to_str() const
{
    return m_impl->to_str(); 
}

table_constraint& table_constraint::primary_key(const column_list& cols, conflict_behaviour cb)
{
    m_impl->primary_key(cols, cb);
    return *this;
}

table_constraint& table_constraint::unique(const column_list& cols, conflict_behaviour cb)
{
    m_impl->unique(cols, cb);
    return *this;
}

table_constraint& table_constraint::check(const expression& e)
{
    m_impl->check(e);
    return *this;
}

table_constraint& table_constraint::foreign_key(const column_list& local_cols, const table& foreign_table,
                                            const column_list& foreign_cols)
{
    m_impl->foreign_key(local_cols, foreign_table, foreign_cols);
    return *this;
}

namespace details
{

table_constraint_impl::table_constraint_impl()
{}

std::string table_constraint_impl::to_str() const
{
    std::string sql;
 
    if( m_pkey )
    {
        sql += "PRIMARY KEY (";
        sql += m_pkey->m_cols.to_str();
        sql += ")";
        sql += cb_to_str(m_pkey->m_cb);
    }
    
    for( unsigned i = 0; i < m_uniques.size(); ++i )
    {
        if( "" != sql ) sql += ", ";
        sql += "UNIQUE (";
        sql += m_uniques[i].m_cols.to_str();
        sql += ")";
        sql += cb_to_str(m_uniques[i].m_cb);
    }
    
    for( unsigned i = 0; i < m_fkeys.size(); ++i )
    {
        if( "" != sql ) sql += ", ";
        sql += "FOREIGN KEY (";
        sql += m_fkeys[i].m_local_cols.to_str();
        sql += ") REFERENCES ";
        sql += m_fkeys[i].m_tab.get_name();
        if( !m_fkeys[i].m_foreign_cols.empty() )
        {
            sql += " (";
            sql += m_fkeys[i].m_foreign_cols.to_str();
            sql += ')';
        }
    }
    
    for( unsigned i = 0; i < m_checks.size(); ++i )
    {
        if( "" != sql ) sql += ", ";
        sql += " CHECK (";
        sql += m_checks[i].to_str();
        sql += ")";
    }
    
    return sql;
}

void table_constraint_impl::primary_key(const column_list& cols, conflict_behaviour cb)
{
    m_pkey.reset( new primary_key_info(cols,cb) );
}

void table_constraint_impl::unique(const column_list& cols, conflict_behaviour cb)
{
    m_uniques.push_back( primary_key_info(cols, cb) );
}

void table_constraint_impl::check(const expression& e)
{
    m_checks.push_back( e );
}

void table_constraint_impl::foreign_key(const column_list& local_cols, const table& foreign_table,
                                        const column_list& foreign_cols)
{
    m_fkeys.push_back( foreign_key_info(local_cols, foreign_table, foreign_cols) );
}

}
}}
