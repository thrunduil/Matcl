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
inline object eval_function::eval(const function_name& func, Object&& ... args)
{
    Type types[]            = {details::get_arg_type_eval<Object>::eval(args) ... };
    int n_args              = sizeof...(Object);


    function f              = operations::get_overload(func, n_args, types);
    const object* args[]    = {(const object*)&args...};

    object ret;
    f.make(n_args, args, ret);

    return ret;
};

inline object eval_function::eval(const function_name& func)
{
    Type* types         = nullptr;
    int n_args          = 0;
    
    function f          = operations::get_overload(func, n_args, types);

    const object** args = nullptr;

    object ret;
    f.make(n_args, args, ret);
    
    return ret;
};

inline eval_function_template::eval_function_template(std::initializer_list<Type> types)
    :m_types{types}
{};

template<class ... Object>
inline object eval_function_template::eval(const function_name& func, Object&& ... args)
{
    Type types[]            = {details::get_arg_type_eval<Object>::eval(args) ... };
    int n_args              = sizeof...(Object);
    int n_types             = (int)m_types.size();

    function f              = operations::get_template_overload(func, n_types, 
                                m_types.data(), n_args, types);

    const object* args[]    = {(const object*)&args...};

    object ret;    
    f.make(n_args, args, ret);

    return ret;
};

inline object eval_function_template::eval(const function_name& func)
{
    Type* types         = nullptr;
    int n_args          = 0;
    int n_types         = (int)m_types.size();

    function f          = operations::get_template_overload(func, n_types, 
                                m_types.data(), n_args, types);

    const object** args = nullptr;
    
    object ret;
    f.make(n_args, args, ret);

    return ret;
};

};};
