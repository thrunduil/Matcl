/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-dynamic/type.h"
#include "matcl-dynamic/details/type_impl.h"
#include "matcl-dynamic/function.h"
#include "matcl-core/IO/archive.h"
#include "type_table.h"
#include "type_reference.h"

namespace matcl { namespace dynamic
{

Type::Type()
    :m_impl(nullptr)
{};

Type::Type(const type_impl* impl)
    :m_impl(impl)
{};

std::string Type::to_string() const
{ 
    return m_impl? m_impl->to_string() : "null"; 
};

void Type::serialize(oarchive_impl & ar, const unsigned int) const
{
    if (m_impl == nullptr)
    {
        ar << "null";
        return;
    };

    std::string name(m_impl->get_class_name());
    ar << name;
};

void Type::serialize(iarchive_impl & ar, const unsigned int)
{    
    std::string name;
    ar >> name;

    *this = details::type_table::get()->get_type(name.c_str());
};

bool Type::is_reference(Type t)
{
    if (t == Type())
        return false;
    
    return t.m_impl->is_reference();
};

Type Type::decay(Type t)
{
    if (t == Type())
        return false;

    return t.m_impl->decay();
};

std::ostream& dynamic::operator<<(std::ostream& os, Type t)
{
    if (t == Type())
        os << '\'' << '\'';
    else
        os << '\'' << t.get_impl()->get_class_name() << '\'';
    return os;
};

std::istream& dynamic::operator>>(std::istream& is, Type& t)
{
    std::stringbuf buf;
    std::string s;
    char c;
    is >> c;
    
    if(c == '\'')
    {
        is.get(buf, '\'');
        s = buf.str();
        is >> c;
    }
    
    if (s.size() == 0)
    {
        t = Type();
        return is;
    };

    t = details::type_table::get()->get_type(s.c_str());

    if (t == Type())
    {
        std::string msg = "unregistered class " + s;
        throw error::matcl_dynamic_exception(msg);
    };

    return is;
};

bool operations::has_trivial_assignment(Type lhs, Type rhs)
{
    if (lhs != rhs)
        return false;

    if (lhs == Type())
        return rhs == Type();

    return lhs.get_impl()->has_trivial_assign();
};

bool operations::has_one(Type t)
{
    if (t == Type())
        return false;

    return t.get_impl()->has_one();
};

Type operations::unify_types(Type t1, Type t2)
{
    if (t1 == t2)
        return t1;

    Type t = details::type_table::get()->unify_types(t1, t2);
    return t;
};

Type operations::make_reference_type(Type t)
{
    return details::type_table::get()->make_reference_type(t);
};

function operations::get_overload(const function_name& func, const Type t[], int n_args)
{
    switch(n_args)
    {
        case 1:
            return details::type_table::get()->get_overload_1(func, t);
        case 2:
            return details::type_table::get()->get_overload_2(func, t);
        default:
            return details::type_table::get()->get_overload_n(func, t, n_args);
    };
};

Type operations::return_type(const function_name& func, int n_args, const Type* t)
{
    function f   = get_overload(func, t ,n_args);
    return f.return_type();
}

function
operations::get_template_overload(const function_name& func, int n_templ, 
                   const Type templates[], int n_args, const Type arg_types[])
{
    return details::type_table::get()->get_template_overload
                        (func,n_templ, templates, n_args, arg_types);
}

Type predefined::type_int()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_int);
    return t; 
};

Type predefined::type_bool()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_bool);
    return t; 
};

Type predefined::type_null()
{
    return Type(); 
};

Type predefined::type_unit()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_unit);
    return t; 
};

Type predefined::type_any()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_any);
    return t; 
};

Type predefined::type_real()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_real);
    return t; 
};

Type predefined::type_float()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_float);
    return t; 
};

Type predefined::type_complex()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_compl);
    return t; 
};

Type predefined::type_float_complex()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_fcompl);
    return t; 
};

Type predefined::type_string()
{
    static Type t = details::type_table::get()->get_predefined(details::type_table::ti_string);
    return t; 
};

};};
