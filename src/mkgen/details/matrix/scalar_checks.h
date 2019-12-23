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

namespace mk = matcl :: mkgen;

has_static_member_template_function_x(print)
has_static_member_template_function_x(eval)
has_static_member_template_function_x(value)
has_static_member_function_x(accept)
has_template_alias_x(simplify)
has_constexpr_static_member_function_x(is_simplified)

has_static_member_function_x(print)

//TODO
struct subs_context_dummy
{
};

struct local_storage_dummy
{
};

struct visitor_dummy
{
};

// return true if T can be passed to ct_scalar<T, ...>;
template<class T>
struct is_valid_scalar_data
{
    static const bool value = std::is_base_of<mkd::scalar_data<T>, T>::value;
};

//-----------------------------------------------------------------------
//                      check_valid_scalar_data
//-----------------------------------------------------------------------
// check if Data parameter supplied to ct_scalar is valid
template<class Data>
struct check_valid_scalar_data
{
    static const bool is_sd = std::is_base_of<mkd::scalar_data<Data>, Data>::value;

    static_assert(is_sd == true, "type is not scalar_data<>");

    using type  = typename Data::template check_scalar_data<void>;
};

// check if Data parameter supplied to ct_scalar is valid
template<class Data, class Ret>
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

    // template<class Visitor>
    // static void accept(Visitor& vis);
    using func_accept_type      = void (visitor_dummy&);
    static const bool has_accept= has_static_member_function_accept
                                    <Data, func_accept_type>::value;

    static_assert(has_accept == true, "Data must implement function accept");

    // template<class Void>
    // using simplify;

    static const bool has_simpl = has_template_alias_simplify<Data, void>::value;
    static_assert(has_simpl == true, "Data must implement template alias: simplify");

    // static constexpr bool is_simplified();

    using func_accept_is_simpl  = bool ();
    static const bool has_issimpl = has_constexpr_static_member_function_is_simplified
                                    <Data, func_accept_is_simpl>::value;
    static_assert(has_issimpl == true, "Data must implement consexpr function: is_simplified");

    using type = Ret;
};

//-----------------------------------------------------------------------
//                      check_valid_const_data_tag
//-----------------------------------------------------------------------
// check if Tag parameter supplied to scal_data_const_value is valid
template<class Tag>
struct check_valid_const_data_tag
{
    using tag_type  = Tag;
    using base_type = mk::scal_data_const_value_tag<tag_type>;
    static const bool is_valid_tag  = std::is_base_of<base_type, tag_type>::value;

    static_assert(is_valid_tag == true, "Tag is not scal_data_const_value_tag<>");

    using type  = typename Tag::template check_scal_data_const_value_tag<void>;
};

// check if Tag parameter supplied to scal_data_const_value_tag is valid
template<class Tag, class Ret>
struct check_const_data_tag_impl
{
    // Data must implement:

    // template<class Val>
    // static constexpr Val value();

    using func_eval_type       = double ();
    static const bool has_eval = has_static_member_template_function_value
                                    <Tag, double, func_eval_type>::value;

    static_assert(has_eval == true, "Tag must implement function value()");

    //check if Tag::value() is constexpr
    static_assert(Tag::template value<double>() != std::numeric_limits<double>::quiet_NaN(), 
                  "Tag::value must be constexpr");

    using type = Ret;
};

//-----------------------------------------------------------------------
//                      check_valid_data_tag
//-----------------------------------------------------------------------
// check if Tag parameter supplied to scal_data_value is valid
template<class Tag>
struct check_valid_data_tag
{
    static const bool is_valid_tag  = std::is_base_of<mk::scal_data_value_tag<Tag>, Tag>
                                                ::value;

    static_assert(is_valid_tag == true, "Tag is not scal_data_value_tag<>");

    using type  = typename Tag::template check_scal_data_value_tag<void>;
};

// check if Tag parameter supplied to scal_data_const_value_tag is valid
template<class Tag, class Ret>
struct check_data_tag_impl
{
    // Data must implement:

    // template<class Val>
    // static Val value();

    using func_eval_type       = double ();
    static const bool has_eval = has_static_member_template_function_value
                                    <Tag, double, func_eval_type>::value;

    static_assert(has_eval == true, "Tag must implement function value()");

    using type = Ret;
};

//-----------------------------------------------------------------------
//                      check_valid_gen_data_tag
//-----------------------------------------------------------------------
// check if Tag parameter supplied to scal_data_gen_value is valid
template<class Tag>
struct check_valid_gen_data_tag
{
    using tag_type  = Tag;
    using base_type = mk::scal_data_gen_value_tag<tag_type>;

    static const bool is_valid_tag  = std::is_base_of<base_type, tag_type>
                                                ::value;

    static_assert(is_valid_tag == true, "Tag is not scal_data_gen_value_tag<>");

    using type  = typename tag_type::template check_scal_data_gen_value_tag<void>;
};

// check if Tag parameter supplied to scal_data_gen_value_tag is valid
template<class Tag, class Ret>
struct check_gen_data_tag_impl
{
    // Data must implement:

    // template<class Val, class Local_Storage>
    // static Val eval(const Local_Storage& ls);

    using func_eval_type       = double (const local_storage_dummy&);
    static const bool has_eval = has_static_member_template_function_eval
                                    <Tag, double, func_eval_type>::value;

    static_assert(has_eval == true, "Tag must implement function eval");

    using type = Ret;
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
