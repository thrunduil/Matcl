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

#include "null_conversions_table.h"
#include "matcl-core/error/exception_classes.h"

namespace matcl { namespace dynamic { namespace details
{

class fun_conv_null : public evaler
{
    private:
        fun_conv_null(const fun_conv_null&) = delete;
        fun_conv_null& operator=(const fun_conv_null&) = delete;

    public:
        fun_conv_null(Type to);

        ~fun_conv_null() override;

        bool        make_eval(const object** _args, object& ret) const override;
        void        make_eval(const object** _args) const override;
        function    make_converter(int n_deduced, const Type deduced[], Type deduced_ret, 
                        const std::vector<function>& arg_converters) const override;
};

fun_conv_null::fun_conv_null(Type to)
{
    m_args_size = 1;
    m_ret_ti    = to;
    m_arg_ti    = new Type();
    m_arg_ti[0] = Type();
}

fun_conv_null::~fun_conv_null()
{
    delete m_arg_ti;
}

bool fun_conv_null::make_eval(const object** _args, object& ret) const
{
    object& obj_null = *const_cast<object*>(_args[0]);
    obj_null.reset(object(m_ret_ti));
    ret.reset(obj_null);  
    return false;
};

void fun_conv_null::make_eval(const object** _args) const
{
    object& obj_null = *const_cast<object*>(_args[0]);
    obj_null.reset(object(m_ret_ti));
    return;
};

function fun_conv_null::make_converter(int, const Type[], Type, 
                             const std::vector<function>&) const
{
    matcl_assert(false, "Should never be here!");
    return function();
};

null_conversion::null_conversion(Type to)
    :m_to(to), m_func(new fun_conv_null(to))
{};

};};};

