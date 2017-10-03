/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "type_table_cache.h"

namespace matcl { namespace dynamic { namespace details
{

type_table_cache::type_table_cache()
    :m_unifiers(false), m_overloads(false), m_template_overloads(false)
    ,m_convert(false), m_assign(false)
{};

Type type_table_cache::get_unifier(Type t1, Type t2) const
{
    unifier_ptr ret = m_unifiers.get_existing(unifier_info{t1,t2}); 
    if (ret.m_ptr)
        return ret.m_ptr->get();
    else
        return Type();
};

Type type_table_cache::set_unifier(Type t1, Type t2, Type t)
{
    m_unifiers.get(unifier_info{t1,t2,t}); 
    return t;
};

function type_table_cache::get_converter(Type to, Type from, converter_type ctype) const
{
    convert_ptr ret = m_convert.get_existing(convert_info{to,from,ctype}); 
    if (ret.m_ptr)
        return ret.m_ptr->get();
    else
        return function();
};

function type_table_cache::set_converter(Type to, Type from, converter_type ctype, 
                                         const function& f)
{
    m_convert.get(convert_info{to,from,ctype, f}); 
    return f;
};

function type_table_cache::get_assigner(Type to, Type from) const
{
    assign_ptr ret = m_assign.get_existing(assign_info{to,from}); 
    if (ret.m_ptr)
        return ret.m_ptr->get();
    else
        return function();
};

function type_table_cache::set_assigner(Type to, Type from, const function& f)
{
    m_assign.get(assign_info{to,from,f}); 
    return f;
};

function type_table_cache::get_overload(const function_name& func, 
                                int n_args, const Type t[]) const
{
    overload_ptr ret = m_overloads.get_existing(overload_info{func,n_args,t}); 
    if (ret.m_ptr)
        return ret.m_ptr->get();
    else
        return function();
};

function type_table_cache::set_overload(const function_name& func, 
                            int n_args, const Type t[], const function& f)
{
    m_overloads.get(overload_info{func,n_args,t, f}); 
    return f;
};

function type_table_cache::get_template_overload(const function_name& func, int n_templ, 
                            const Type templates[], int n_args, const Type args[]) const
{
    toverload_ptr ret = m_template_overloads.get_existing
                            (toverload_info{func,n_templ,templates,n_args,args}); 

    if (ret.m_ptr)
        return ret.m_ptr->get();
    else
        return function();
};

function type_table_cache::set_template_overload(const function_name& func, int n_templ,
                    const Type templates[], int n_args, const Type t[], const function& f)
{
    m_template_overloads.get(toverload_info{func,n_templ, templates, n_args,t, f}); 
    return f;
};

void type_table_cache::clear()
{
    m_unifiers.clear();
    m_overloads.clear();
    m_convert.clear();
    m_assign.clear();
    m_template_overloads.clear();
};

};};};
