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

//--------------------------------------------------------------------------
//                         last_call_info
//--------------------------------------------------------------------------

inline last_call_info::last_call_info()
    :n_args(0)
{};

inline last_call_info::last_call_info(int args, Type t1, Type t2, function f)
    :n_args(args), m_ty_1(t1), m_ty_2(t2), m_func(f)
{};

//--------------------------------------------------------------------------
//                         last_call_cache
//--------------------------------------------------------------------------
inline void last_call_cache::clear()
{
    m_cache = vec_type();
}

force_inline
function last_call_cache::get_last_function(size_t code, int n_args, Type t1, Type t2) const
{
    if (code >= m_cache.size())
        return function();

    const last_call_info& info  = m_cache[code];

    bool eq1    = info.n_args   == n_args;
    bool eq2    = info.m_ty_1   == t1;
    bool eq3    = info.m_ty_2   == t2;
    
    if (eq1 && eq2 && eq3)
        return info.m_func;
    else
        return function();
};

force_inline
void last_call_cache::set_last_function(size_t code, int n_args, Type t1, Type t2, 
                    function f)
{
    if (code >= m_cache.size())
        m_cache.resize(code + 1);

    m_cache[code] = last_call_info(n_args, t1, t2, f);
};

//--------------------------------------------------------------------------
//                         type_table_cache
//--------------------------------------------------------------------------
inline
Type type_table_cache::get_unifier(Type t1, Type t2) const
{
    unifier_ptr ret = m_unifiers.get_existing(unifier_info{t1,t2}); 
    if (ret.m_ptr)
        return ret.m_ptr->get();
    else
        return Type();
};

inline
void type_table_cache::set_unifier(Type t1, Type t2, Type t)
{
    m_unifiers.get(unifier_info{t1,t2,t}); 
};

inline
function type_table_cache::get_converter(Type to, Type from, converter_type ctype) const
{
    convert_ptr ret = m_convert.get_existing(convert_info{to,from,ctype}); 

    if (ret.m_ptr)
        return ret.m_ptr->get();
    else
        return function();
};

inline
function type_table_cache::set_converter(Type to, Type from, converter_type ctype, 
                                         function f)
{
    convert_ptr ret = m_convert.get(convert_info{to,from,ctype, f}); 
    return ret.m_ptr->get();
};

inline
function type_table_cache::get_assigner(Type to, Type from) const
{
    assign_ptr ret = m_assign.get_existing(assign_info{to,from}); 

    if (ret.m_ptr)
        return ret.m_ptr->get();
    else
        return function();
};

inline
function type_table_cache::set_assigner(Type to, Type from, function f)
{
    assign_ptr ret = m_assign.get(assign_info{to,from,f}); 
    return ret.m_ptr->get();
};

force_inline
function type_table_cache::get_overload(const function_name& func, 
                                int n_args, const Type t[])
{
    Type t1, t2;

    switch(n_args)
    {
        case 2:
        {
            t1  = t[0];
            t2  = t[1];
            break;
        }
        case 1:
        {
            t1  = t[0];
            break;
        }
    }

    function f;

    size_t code = func.get_unique_code();    
    f           = m_last_call.get_last_function(code, n_args, t1, t2);

    if (f.is_null() == false)
        return f;

    overload_ptr ret = m_overloads.get_existing(overload_info{func,n_args,t}); 

    if (ret.m_ptr)
    {
        f = ret.m_ptr->get();
        m_last_call.set_last_function(code, n_args, t1, t2, f);
    }

    return f;
};

inline
function type_table_cache::set_overload(const function_name& func, 
                            int n_args, const Type t[], function f)
{
    overload_ptr ret = m_overloads.get(overload_info{func,n_args,t, f}); 
    return ret.m_ptr->get();
};

inline
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

inline
function type_table_cache::set_template_overload(const function_name& func, int n_templ,
                    const Type templates[], int n_args, const Type t[], function f)
{
    toverload_ptr ret = m_template_overloads.get(toverload_info{func,n_templ, templates, n_args,t, f}); 
    return ret.m_ptr->get();
};

};};};
