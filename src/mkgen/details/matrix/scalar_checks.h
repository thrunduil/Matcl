/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019
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

#include "mkgen/details/mkgen_fwd.h"
#include "mkgen/details/utils/has_function.h"
#include "matcl-core/details/mpl.h"

namespace matcl { namespace mkgen { namespace details
{

has_static_member_template_function_x(print)
has_static_member_template_function_x(eval)

struct subs_context_dummy
{
};

struct local_storage_dummy
{
};

//-----------------------------------------------------------------------
//                      check_scalar_data
//-----------------------------------------------------------------------
template<class Data>
struct check_valid_scalar_data
{
    static_assert(matcl::details::dependent_false<Data>::value, "type is not scalar_data<>");
};

template<class Data>
struct check_valid_scalar_data<scalar_data<Data>>
{
    // Data is already checked
    using type = void;
};

// check if Data parameter supplied to ct_scalar is valid
template<class Data>
struct check_scalar_data
{
    // Data must be specialization of scalar_data<T>
    using type = typename check_valid_scalar_data<Data>::type;
};

// check if Data parameter supplied to ct_scalar is valid
template<class Data>
struct check_scalar_data_impl
{
    // Data must implement:

    // template<class Subs_Context>
    // static void print(std::ostream& os, int prior)
    using func_print_type       = void (std::ostream& os, int prior);
    static const bool has_print = has_static_member_template_function_print
                                    <Data, subs_context_dummy, func_print_type>::value;

    static_assert(has_print == true, "Data must implement function print");

    // template<class Val, class Local_Storage>
    // static Val eval(const Local_Storage& ls);

    using func_eval_type       = double (const local_storage_dummy&);
    static const bool has_eval = has_static_member_template_function_eval
                                    <Data, double, func_eval_type>::value;

    static_assert(has_eval == true, "Data must implement function eval");


    //TODO
    using type = void;
};

//-----------------------------------------------------------------------
//                      check_scalar_deps
//-----------------------------------------------------------------------
// check if Deps parameter supplied to ct_scalar is valid
template<class Deps>
struct check_scalar_deps
{
    using type = typename check_valid_dps<Deps>::type;
};

//-----------------------------------------------------------------------
//                      check_computation_tag
//-----------------------------------------------------------------------
// check tag supplied to compute function
template<class Tag>
struct check_computation_tag
{
    // TODO: collect requirements
    using type = void;
};

}}}
