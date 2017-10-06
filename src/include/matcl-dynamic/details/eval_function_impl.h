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

#include "matcl-core/config.h"
#include "matcl-dynamic/details/fwd_decls.h"

#include <vector>

#pragma warning(push)
#pragma warning(disable : 4251) // needs to have dll-interface to be used by clients

namespace matcl { namespace dynamic { namespace details
{

struct eval_function_1
{
    // object 'ret' cannot be initialized (otherwise memory leaks are possible)
    // types and args are arrays of size 1; types[i] == args[i]->get_type() 
    MATCL_DYN_EXPORT
    static void eval(const function_name& func, object& ret, const Type* types,
                    const object ** args, int);
};

struct eval_function_2
{
    // object 'ret' cannot be initialized (otherwise memory leaks are possible)
    // types and args are arrays of size 2; types[i] == args[i]->get_type() 
    MATCL_DYN_EXPORT
    static void eval(const function_name& func, object& ret, const Type* types,
                    const object ** args, int);
};

struct eval_function_n
{
    // object 'ret' cannot be initialized (otherwise memory leaks are possible)
    // types and args are arrays of size n_args; types[i] == args[i]->get_type() 
    MATCL_DYN_EXPORT
    static void eval(const function_name& func, object& ret, const Type* types,
                    const object ** args, int n_args);
};

template<int N_args>
struct      select_evaler           { using type = eval_function_n; };
template<>  struct select_evaler<1> { using type = eval_function_1; };
template<>  struct select_evaler<2> { using type = eval_function_2; };


};};};

#pragma warning(pop)

