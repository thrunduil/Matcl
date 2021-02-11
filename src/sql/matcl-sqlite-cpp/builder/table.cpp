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

#include "matcl-sqlite-cpp/builder/table.h"
#include "sql/matcl-sqlite-cpp/builder/details/table_impl.h"

namespace matcl { namespace sql
{

table::table(const std::string& name, const std::string& alias)
    : join_source(impl_ptr(new details::table_impl(name,alias)))
{}

std::string table::get_ref() const
{
    return m_impl->get_ref();
}

std::string table::get_name() const
{
    return m_impl->get_name();
}

join_source::join_source(const impl_ptr& impl)
    : m_impl(impl)
{}

std::string join_source::to_str() const
{
    return m_impl->to_str();
}

join_source join_source::join(const join_source& t, bool natural, join_type jt)
{
    impl_ptr impl_copy(new details::table_impl(*m_impl));
    impl_copy->join(*t.m_impl,natural,jt);
    return join_source(impl_copy);
}

join_source join_source::join(const join_source& t, const column_list& cl, join_type jt)
{
    impl_ptr impl_copy(new details::table_impl(*m_impl));
    impl_copy->join(*t.m_impl,cl,jt);
    return join_source(impl_copy);
}

join_source join_source::join(const join_source& t, const expression& e, join_type jt)
{
    impl_ptr impl_copy(new details::table_impl(*m_impl));
    impl_copy->join(*t.m_impl,e,jt);
    return join_source(impl_copy);
}

namespace details
{

table_impl::table_impl(const std::string& name, const std::string& alias)
    : m_name(name), m_alias(alias)
{}

std::string table_impl::to_str() const
{
    std::string txt = m_name;
 
    if( "" != m_alias )
    {
        txt += " AS ";
        txt += m_alias;
    }

    for(join_info j : m_joins)
    {
        txt += ' ';
        txt += j.m_join->op();
        txt += ' ';
        txt += j.m_joinee;
        std::string constraint = j.m_join->constraint();
    
        if("" != constraint)
        {
            txt += ' ';
            txt += constraint;
        }
    }

    return txt;
}

std::string table_impl::get_ref() const
{
    return "" == m_alias ? m_name : m_alias;
}

std::string table_impl::get_name() const
{
    return m_name;
}

void table_impl::join(const table_impl& t, join_def_base* jdef)
{
    join_info ji;
    ji.m_joinee = t.m_name;
    
    if( "" != t.m_alias )
    {
        ji.m_joinee += " AS ";
        ji.m_joinee += t.m_alias;
    }
    
    ji.m_join.reset(jdef);
    m_joins.push_back(ji);
    m_joins.insert( m_joins.end(), t.m_joins.begin(), t.m_joins.end() );
}

void table_impl::join(const table_impl& t, bool natural, join_type jt)
{
    join(t, natural ? new details::natural_join(jt) : new details::simple_join(jt));
}

void table_impl::join(const table_impl& t, const column_list& c, join_type jt)
{
    join(t, new details::column_join(c, jt));
}

void table_impl::join(const table_impl& t, const expression& e, join_type jt)
{
    join(t, new details::expr_join(e, jt));
}

simple_join::simple_join(join_type jt)
    : m_jt(jt)
{}

std::string simple_join::op() const
{
    switch(m_jt)
    {
        case join_type::left:
            return "LEFT JOIN";
        case join_type::full:
            return "JOIN";
    }

    assert(!"unexpected join_type");
    return "JOIN"; // fallback to something that makes sense
}

std::string simple_join::constraint() const
{
    return "";
}

natural_join::natural_join(join_type jt)
    : simple_join(jt)
{}

std::string natural_join::op() const
{
    return "NATURAL " + simple_join::op();
}

column_join::column_join(const column_list& cl, join_type jt)
    : simple_join(jt), m_cols(cl)
{
    assert( !m_cols.empty() );
}

std::string column_join::constraint() const
{
    return "USING ( " + m_cols.to_str() + " )";
}

expr_join::expr_join(const expression& e, join_type jt)
    : simple_join(jt), m_expr(e)
{}

std::string expr_join::constraint() const
{
    return "ON " + m_expr.to_str();
}

}
}}
