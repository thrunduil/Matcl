/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2018
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

#include "matcl-sqlite-cpp/builder/ord_expression.h"
#include "sql/matcl-sqlite-cpp/builder/details/ord_expression_impl.h"

namespace matcl { namespace sql
{

ord_expression::ord_expression(const expression& e, ordering o)
    : m_impl(new details::ord_expression_impl(e,o))
{}

std::string ord_expression::to_str() const
{
    return m_impl->to_str();
}

namespace details
{

ord_expression_impl::ord_expression_impl(const expression& e, ordering o)
    : m_expr(e), m_ord(o)
{}

std::string ord_expression_impl::to_str() const
{
    std::string txt = m_expr.to_str();
    switch(m_ord)
    {
        case ordering::desc:
            txt += " DESC";
            break;
    }
    return txt;
}

}
}}
