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

#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"
#include "matcl-core/general/fwd_decls.h"

#include "matcl-dynamic/function.h"
#include "matcl-dynamic/type.h"

namespace matcl { namespace dynamic
{

function function::make_converter(int n_deduced, const Type deduced[], Type deduced_ret, 
                                  const std::vector<function>& arg_converters) const
{
    if (!m_evaler)
        return function();

    return m_evaler->make_converter(n_deduced, deduced, deduced_ret, arg_converters);
};

object function::make(int n_args, const object* args[]) const
{
    (void)n_args;

    object ret;
    m_evaler->make_eval(args, ret);

    return ret;
};

};};

