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

#include "func_conv_link.h"
#include "matcl-dynamic/details/func_evaler_conv.h"

#include <cassert>

namespace matcl { namespace dynamic { namespace details
{

fun_conv_link::fun_conv_link(const std::vector<function>& convs)
    :m_converters(convs)
{
    m_args_size = 1;
    m_ret_ti    = convs.back().return_type();
    m_arg_ti    = new Type();
    m_arg_ti[0] = convs.front().argument_type(0);
};

bool fun_conv_link::make_eval(const object** _args, object& ret) const
{
    size_t n = m_converters.size();

    const object* args = _args[0];

    for (size_t i = 0; i < n; ++i)
    {
        const evaler* ev = m_converters[i].get_evaler();
        ev->make_eval(&args, ret);
        args = &ret;
    };

    return true;
}

function fun_conv_link::make_converter(int n_deduced, const Type deduced[], 
                    Type ded_ret, const std::vector<function>& conv_vec) const
{
    return new details::fun_evaler_conv(this,n_deduced, deduced, ded_ret, conv_vec);
};

};};};