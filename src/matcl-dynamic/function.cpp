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

#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"

#include "matcl-dynamic/function.h"
#include "matcl-dynamic/type.h"
#include "type_table.h"

namespace matcl { namespace dynamic
{

void function::make(int n_args, const object* args[], object& ret) const
{
    (void)n_args;
    m_evaler->make_eval(args, ret);
};

void function::make(int n_args, const object* args[]) const
{
    (void)n_args;
    m_evaler->make_eval(args);
};

void details::eval_function_1::eval(const function_name& func, object& ret, const Type* types,
                const object ** args, int)
{
    function f  = details::type_table::get()->get_overload_1(func, types);
    f.get_evaler()->make_eval(args, ret);
};

void details::eval_function_2::eval(const function_name& func, object& ret, const Type* types,
                const object ** args, int)
{
    function f  = details::type_table::get()->get_overload_2(func, types);
    f.get_evaler()->make_eval(args, ret);
};

void details::eval_function_n::eval(const function_name& func, object& ret, const Type* types,
                const object ** args, int n_args)
{
    function f  = details::type_table::get()->get_overload_n(func, types, n_args);
    f.get_evaler()->make_eval(args, ret);
};

void eval_function_template::eval_impl(const function_name& func, object& ret, const Type* types,
            const object ** args, int n_args)
{
    int n_types = (int)m_types.size();
    function f  = details::type_table::get()->get_template_overload
                        (func, n_types, m_types.data(), n_args, types);

    f.get_evaler()->make_eval(args, ret);
};

};};

