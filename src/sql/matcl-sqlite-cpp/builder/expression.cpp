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

#include "matcl-sqlite-cpp/builder/expression.h"
#include "sql/matcl-sqlite-cpp/builder/details/expression_impl.h"
#include <sstream>

namespace matcl { namespace sql
{

null_t null;

expression::expression()
{}

expression::expression(const column& c)
    : m_impl(new details::column_expression(c))
{}

expression::expression(const bind_param& b)
    : m_impl(new details::bind_param_expression(b))
{}

expression::expression(double d)
    : m_impl(new details::numeric_expression(d))
{}

expression::expression(const std::string& s)
    : m_impl(new details::string_expression(s))
{}

expression::expression(const char* c)
    : m_impl(new details::string_expression(c))
{}

expression::expression(const null_t)
    : m_impl(new details::null_expression())
{}

expression::expression(details::expression_base_impl* impl)
    : m_impl(impl)
{}

std::string expression::to_str() const
{
    return m_impl->to_str();
}

bool expression::exists() const
{
    return m_impl ? true : false ;
}

#define CPP2SQL1(FUNC, OP)                                      \
expression FUNC(const expression& e)                            \
{                                                               \
    return expression(new details::unary_op_expression(e,OP));  \
}

CPP2SQL1(operator-, "-");
CPP2SQL1(operator+, "+");
CPP2SQL1(operator~, "~");
CPP2SQL1(operator!, "NOT");

#define CPP2SQL2(FUNC, OP)                                      \
expression FUNC(const expression& e1, const expression& e2)     \
{                                                               \
    return expression(new details::binary_op_expression(e1,e2,OP));\
}

CPP2SQL2(concat, "||");
CPP2SQL2(operator*, "*");
CPP2SQL2(operator/, "/");
CPP2SQL2(operator%, "%");
CPP2SQL2(operator+, "+");
CPP2SQL2(operator-, "-");
CPP2SQL2(operator<<, "<<");
CPP2SQL2(operator>>, ">>");
CPP2SQL2(operator&, "&");
CPP2SQL2(operator|, "|");
CPP2SQL2(operator<, "<");
CPP2SQL2(operator<=, "<=");
CPP2SQL2(operator>, ">");
CPP2SQL2(operator>=, ">=");
CPP2SQL2(operator==, "==");
CPP2SQL2(operator!=, "!=");
CPP2SQL2(is, "IS");
CPP2SQL2(is_not, "IS NOT");
CPP2SQL2(in, "IN");
CPP2SQL2(like, "LIKE");
CPP2SQL2(glob, "GLOB");
CPP2SQL2(match, "MATCH");
CPP2SQL2(regexp, "REGEXP");
CPP2SQL2(operator&&, "AND");
CPP2SQL2(operator||, "OR");

expression sql::call(const std::string& func, wildcard)
{
    expression_list el;
    el.push_back( expression(new details::star_expression()) );
    return expression(new details::call_expression(func, el));
}

expression sql::call(const std::string& func)
{
    expression_list el;
    return expression(new details::call_expression(func, el));
}

expression sql::call(const std::string& func, const expression& e1)
{
    expression_list el;
    el.push_back(e1);
    return expression(new details::call_expression(func, el));
}

expression sql::call(const std::string& func, const expression& e1, const expression& e2)
{
    expression_list el;
    el.push_back(e1);
    el.push_back(e2);
    return expression(new details::call_expression(func, el));
}

expression sql::call(const std::string& func, const expression& e1, const expression& e2, 
                const expression& e3)
{
    expression_list el;
    el.push_back(e1);
    el.push_back(e2);
    el.push_back(e3);
    return expression(new details::call_expression(func, el));
}

expression sql::call(const std::string& func, const expression_list& el, bool distinct)
{
    return expression(new details::call_expression(func, el, distinct));
}

expression sql::collate(const expression& e, collation c)
{
    return expression(new details::collate_expression(e, c));
}

namespace details
{

column_expression::column_expression(const column& c)
    : m_c(c)
{}

std::string column_expression::to_str() const
{
    return m_c.to_str();
}

numeric_expression::numeric_expression(double d)
    : m_num(d)
{}

std::string numeric_expression::to_str() const
{
    std::stringstream ss;
    ss << m_num;
    return ss.str();
}

string_expression::string_expression(const std::string& s)
    : m_str(s)
{}

std::string string_expression::to_str() const
{
    return '\'' + m_str + '\'';
}

std::string null_expression::to_str() const
{
    return "NULL";
}

bind_param_expression::bind_param_expression(const bind_param& b)
    : m_bp(b)
{}

std::string bind_param_expression::to_str() const
{
    return m_bp.to_str();
}

binary_op_expression::binary_op_expression(const expression& e1, const expression& e2, 
                                            const std::string& op)
    : m_e1(e1), m_e2(e2), m_op(op)
{}

std::string binary_op_expression::to_str() const
{
    return '(' + m_e1.to_str() + ' ' + m_op + ' ' + m_e2.to_str() + ')';
}

unary_op_expression::unary_op_expression(const expression& e, const std::string& op)
    :m_e(e), m_op(op)
{}

std::string unary_op_expression::to_str() const
{
    return '(' + m_op + m_e.to_str() + ')';
}

call_expression::call_expression(const std::string& name, const expression_list& el, bool distinct)
    : m_name(name), m_args(el), m_distinct(distinct)
{
}

std::string call_expression::to_str() const
{
    std::string call_text = m_name;
    call_text += '(';

    if( !m_args.empty() )
    {
        if( m_distinct )
            call_text += "DISTINCT ";
    
        expression_list::const_iterator it = m_args.begin();
        call_text += it->to_str();
        
        for(++it; it != m_args.end(); ++it)
        {
            call_text += ", ";
            call_text += it->to_str();
        }
    }

    call_text += ')';
    return call_text;
}

std::string star_expression::to_str() const
{
    return "*";
};

collate_expression::collate_expression(const expression& e, collation c)
    : m_expr(e), m_collation(c)
{}

std::string collate_expression::to_str() const
{
    std::string txt = m_expr.to_str();
    return txt += collation_to_str(m_collation);
}

std::string collation_to_str(collation c)
{
    std::string txt;
    
    if(collation::no_collation != c)
        txt += " COLLATE";
 
    switch(c)
    {
        case collation::binary:
            txt += " BINARY";
            break;
        case collation::nocase:
            txt += " NOCASE";
            break;
        case collation::rtrim:
            txt += " RTRIM";
            break;
    }
    return txt;
}

}
}}
