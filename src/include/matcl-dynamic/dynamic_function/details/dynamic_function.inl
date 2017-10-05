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

#include "function_evaler.h"

namespace matcl { namespace dynamic { namespace details
{

template<typename Fun, class base_type>
template<class data_constructor>
force_inline
void dynamic_function<Fun,base_type>::eval(data_constructor arg) const
{
    static const 
    int n_arg       = function_traits::n_inputs;

    static const 
    bool void_return= std::is_same<void,typename function_traits::return_type>::value;

    using ds_type   = typename details::get_reference_type<data_constructor>::type;
	using ret_t     = typename function_traits::return_type;

    ds_type m_ds    = details::get_reference_type<data_constructor>::get_value(arg);

    function_evaler<n_arg,void_return,ds_type,function_traits>::eval(m_ds, m_fun);
};

template<typename Fun, class base_type>
template<class class_type,class data_constructor>
force_inline
void dynamic_function<Fun,base_type>::eval(class_type object, data_constructor arg) const
{
    static const 
    int n_arg       = function_traits::n_inputs;

    static const 
    bool void_return= std::is_same<void,typename function_traits::return_type>::value;

    using ds_type   = typename details::get_reference_type<data_constructor>::type;
    using cl_type   = typename details::get_reference_type<class_type>::type;

    cl_type m_cl    = details::get_reference_type<class_type>::get_value(object);
    ds_type m_ds    = details::get_reference_type<data_constructor>::get_value(arg);

    details::member_evaler<n_arg,void_return,ds_type,function_traits,cl_type>::eval(m_ds,m_fun,m_cl);
};

template<typename Fun, class base_type>
dynamic_function<Fun,base_type>::dynamic_function(Fun f)
	: m_fun(f)
{};

};};};