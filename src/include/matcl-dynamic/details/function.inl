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

#pragma once

#include "matcl-dynamic/function.h"

namespace matcl { namespace dynamic
{

inline function::function()
    :m_evaler()
{};        

inline bool function::is_null() const
{ 
    return m_evaler == nullptr;
};

inline int function::number_arguments() const
{ 
    return m_evaler? m_evaler->m_args_size : 0; 
};

inline Type function::argument_type(int narg) const
{ 
    return m_evaler? m_evaler->m_arg_ti[narg] : Type();
};

inline Type function::return_type() const
{ 
    return m_evaler? m_evaler->m_ret_ti : Type(); 
};

inline const function::ptr_type& function::get_evaler() const
{ 
    return m_evaler; 
};

inline function::function(evaler* ev)
    :m_evaler(ev)
{};

template<class ... Object>
force_inline 
void eval_function::eval(const function_name& func, object& ret, Object&& ... in_args)
{
    static const int n      = sizeof...(Object);

    const object* args[n]   = {(const object*)(&in_args) ... };
    Type types[n]           = {details::get_arg_type_eval<Object>::eval(in_args) ... };    

    using evaler            = typename details::select_evaler<n>::type;
    evaler::eval(func, ret, types, args, n);
};

force_inline
void eval_function::eval(const function_name& func, object& ret)
{
    Type* types         = nullptr;
    const object** args = nullptr;
    
    details::eval_function_n::eval(func, ret, types, args, 0);
};

inline eval_function_template::eval_function_template(std::initializer_list<Type> types)
    :m_types{types}
{};

template<class ... Object>
force_inline
void eval_function_template::eval(const function_name& func, object& ret, Object&& ... in_args)
{
    Type types[]            = {details::get_arg_type_eval<Object>::eval(in_args) ... };
    const object* args[]    = {(const object*)&in_args...};
    int n_args              = sizeof...(Object);

    eval_impl(func, ret, types, args, n_args);
};

force_inline 
void eval_function_template::eval(const function_name& func, object& ret)
{
    Type* types         = nullptr;
    const object** args = nullptr;
    int n_args          = 0;

    eval_impl(func, ret, types, args, n_args);
};

};};
