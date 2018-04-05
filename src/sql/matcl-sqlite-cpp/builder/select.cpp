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

#include "matcl-sqlite-cpp/builder/Select.h"
#include "sql/matcl-sqlite-cpp/builder/details/select_impl.h"

namespace matcl { namespace sql
{

namespace details
{

select_template_base::select_template_base()
    : m_impl(new details::select_impl)
{}

select_template_base::select_template_base(const impl_ptr& impl)
    : m_impl(impl)
{}

std::string select_template_base::to_str() const
{
    return m_impl->to_str();
}

select_impl::select_impl()
    : m_source(table("")), m_source_exists(false), m_distinct(false)
{}

std::string select_impl::to_str() const
{
    std::string query = to_str_internal();

    for(unsigned i = 0; i < m_compounds.size(); ++i)
    {
        switch(m_compounds[i].second)
        {
            case compose_type::union_type:
                query += " UNION ";
            case compose_type::union_all:
                query += "ALL ";
                break;
            case compose_type::except:
                query += " EXCEPT ";
                break;
            case compose_type::intersect:
                query += " INTERSECT ";
                break;
        }

        query += m_compounds[i].first->to_str_internal();
    }

    if(!m_ordering.empty())
    {
        query += " ORDER BY ";
        query += m_ordering.to_str();
    }

    if(m_limit.exists())
    {
        query += " LIMIT ";
        query += m_limit.to_str();

        if(m_offset.exists())
        {
            query += " OFFSET ";
            query += m_offset.to_str();
        }
    }

    return query + "\n";
}

void select_impl::print_col(std::string& query, const result_column& rc) const
{
    query += " ";
    query += rc.first.to_str();

    if( "" != rc.second )
    {
        query += " AS ";
        query += rc.second;
    }
}

std::string select_impl::to_str_internal() const
{
    std::string query = "SELECT";

    if(m_distinct)
        query += " DISTINCT";

    if(!m_cols.empty())
    {
        result_column_list::const_iterator it = m_cols.begin();
        print_col(query, *it);

        for(++it; it != m_cols.end(); ++it)
        {
            query += ",";
            print_col(query, *it);
        }
    }
    else
    {
        query += " *";
    }

    if(m_source_exists)
    {
        query += " FROM ";
        query += m_source.to_str();
    }

    if(m_cond.exists())
    {
        query += " WHERE ";
        query += m_cond.to_str();
    }

    if(!m_grouping.empty())
    {
        query += " GROUP BY ";
        query += m_grouping.to_str();
        if(m_having.exists())
        {
            query += " HAVING ";
            query += m_having.to_str();
        }
    }

    return query;
}

void select_impl::from(const join_source& js)
{
    m_source = js;
    m_source_exists = true;
}

void select_impl::where(const expression& e)
{
    m_cond = e;
}

void select_impl::add_column(const expression& e, const std::string& alias)
{
    m_cols.push_back(make_pair(e, alias));
}

void select_impl::add_table_columns(const table& t)
{
    column col("*");
    col.set_table(t);
    m_cols.push_back(std::make_pair(col, ""));
}

void select_impl::distinct(bool b)
{
    m_distinct = b;
}

void select_impl::group_by(const expression_list& e, const expression &constraint)
{
    m_grouping = e;
    m_having = constraint;
}

void select_impl::order_by(const ord_expression_list& o)
{
    m_ordering = o;
}

void select_impl::limit(const expression& e, const expression& offset)
{
    m_limit = e;
    m_offset = offset;
}

void select_impl::compose(const std::shared_ptr<select_base_impl>& impl_base, compose_type c)
{
    auto impl = std::static_pointer_cast<select_impl>(impl_base);
    m_compounds.push_back(std::make_pair(impl, c));
    m_compounds.insert(m_compounds.end(), impl->m_compounds.begin(), impl->m_compounds.end());
}

}
}}
