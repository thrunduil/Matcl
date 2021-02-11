/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2019 - 2021
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

#include "matcl-core/matrix/scalar_types.h"

#include <iosfwd>

namespace matcl { namespace mkgen 
{

namespace mk    = matcl :: mkgen;

//-----------------------------------------------------------------------
//                      forward declarations
//-----------------------------------------------------------------------
template<class... T>
struct dps;

}};

namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                      forward declarations
//-----------------------------------------------------------------------
template<class Data, Integer Row, Integer Col>
struct matrix_array_get_elem;

//-----------------------------------------------------------------------
//                      check helpers
//-----------------------------------------------------------------------

//TODO
struct subs_context_dummy
{};

struct local_storage_dummy
{};

struct visitor_dummy
{};

// check if v is compile time constant
template<class T, T v>
struct check_constexpr_helper
{};

template<class T, class V, Integer Row, Integer Col>
concept has_static_member_value = 
    requires { {T::template value<V, Row, Col>()} -> std::convertible_to<V>; };

// TODO: check constexpr
template<class T, class V, Integer Row, Integer Col>
concept has_static_member_value_constexpr = 
    requires { {T::template value<V, Row, Col>()} -> std::convertible_to<V>; };

template<class T, class V>
concept has_static_member_value2 = 
    requires { {T::template value<V>()} -> std::convertible_to<V>; };

// TODO: check constexpr
template<class T, class V>
concept has_static_member_value2_constexpr = 
    requires { {T::template value<V>()} -> std::convertible_to<V>; };

template<class T, Integer Row, Integer Col>
concept has_static_member_get_element_impl = 
    requires { typename T :: template get_element_impl<Row, Col>; };

template<class T, class Subs_ctx>
concept has_static_member_template_function_print = 
    requires { T :: template print<Subs_ctx>(std::declval<std::ostream&>(), 
                                                std::declval<int>()); };

template<class T>
concept has_static_member_function_print = 
    requires { T :: template print(std::declval<std::ostream&>(), 
                                                std::declval<int>()); };

template<class T, class Val, class Local_Storage>
concept has_static_member_template_function_eval = 
    requires { {T :: template eval<Val, Local_Storage>
                (std::declval<const Local_Storage&>())} -> std::convertible_to<Val>; };

template<class T, class Visitor>
concept has_static_member_function_accept = 
    requires { T :: template accept<Visitor>(std::declval<Visitor&>()); };

template<class T, class Void>
concept has_template_alias_simplify = 
    requires { T :: template simplify<Void>; };

template<class T>
concept has_constexpr_static_member_function_is_simplified = 
    requires { T :: is_simplified(); 
               typename check_constexpr_helper<bool, T :: is_simplified()>; 
             };

//-----------------------------------------------------------------------
//                      matrix_array
//-----------------------------------------------------------------------
// base class for ct_matrix arrays
template<class Data>
struct matrix_array
{
    // ct_matrix arrays must implement:
    //
    //  template<Integer Row, Integer Col>
    //  using get_element_impl   = [impl]
    // 
    // where T := get_element_impl<Row, Col>::type return type of element at (Row, Col)
    // T must be derived from scalar_data<T> (and satisfy scalar_data requirements)

    template<Integer Row, Integer Col>
    using get_element  = matrix_array_get_elem<Data, Row, Col>;
};

template<class Arr>
struct matrix_array_check_impl
{
    static const bool value = has_static_member_get_element_impl<Arr, 1, 1>;

    // type T = get_element_impl<Row,Col>::type will be checked, when fully constructed
};

struct matrix_array_check
{
    // check if Arr has interface required for matrix_array

    template<class Arr>
    static const bool is_valid = matrix_array_check_impl<Arr> :: value;
};

//-----------------------------------------------------------------------
//                      scalar_data
//-----------------------------------------------------------------------
// base class for ct_scalar arrays
template<class Data>
struct scalar_data
{
    // ct_scalar arrays must implement:
    // template<class Subs_Context>
    // static void print(std::ostream& os, int prior);
    // 
    // template<class Val, class Local_Storage>
    // static Val eval(const Local_Storage& ls);
    //
    // template<class Visitor>
    // static void accept(Visitor& vis);
    //
    // template<class Void>
    // using simplify;
    //
    // template<class Void>
    // static constexpr bool is_simplified;

    //TODO: remove or implement
    // append to Arr_List all arrays required by this scalar
    template<Integer Step, class Arr_List>
    using get_arrays    = Arr_List;
};

template<class Arr>
struct scalar_data_check_impl
{
    static const bool value 
        = mkd::has_static_member_template_function_print<Arr, subs_context_dummy>
        && mkd::has_static_member_template_function_eval<Arr, double, local_storage_dummy>
        && mkd::has_static_member_function_accept<Arr, visitor_dummy>
        && mkd::has_template_alias_simplify<Arr, void>
        && mkd::has_constexpr_static_member_function_is_simplified<Arr>;
};

struct scalar_data_check
{
    // check if Arr has interface required for scalar_data

    template<class Arr>
    static const bool is_valid = scalar_data_check_impl<Arr> :: value;
};

//-----------------------------------------------------------------------
//                      dependencies
//-----------------------------------------------------------------------
template<class Deps>
struct dps_check
{
    static const bool value = false;
};

template<class ...T>
struct dps_check<mk::dps<T...>>
{
    // arguments T are already checked
    // TODO
    static const bool value = true;
};

}}}

namespace matcl { namespace mkgen 
{

//-----------------------------------------------------------------------
//                      matrix_data_const_value_tag
//-----------------------------------------------------------------------
// base class for Tags used in creating const_value_mat
template<class Tag>
struct matrix_data_const_value_tag
{
    // Tag must implement:
    // 
    // template<class Val, Integer Row, Integer Col>
    // static constexpr Val value();
};

template<class Tag>
struct matrix_data_const_value_tag_check_impl
{
    static const bool value = mkd::has_static_member_value_constexpr
                                    <Tag, double, 1, 1>;
};

struct matrix_data_const_value_tag_check
{
    // check if Arr has interface required for matrix_data_const_value_tag

    template<class Tag>
    static const bool is_valid = matrix_data_const_value_tag_check_impl<Tag> :: value;
};

//-----------------------------------------------------------------------
//                      matrix_data_value_tag
//-----------------------------------------------------------------------
// base class for Tags used in creating value_mat
template<class Tag>
struct matrix_data_value_tag
{
    // Tag must implement:
    // 
    // template<class Val, Integer Row, Integer Col>
    // static Val value();
};

template<class Tag>
struct matrix_data_value_tag_check_impl
{
    static const bool value = mkd::has_static_member_value<Tag, double, 1, 1>;
};

struct matrix_data_value_tag_check
{
    // check if Arr has interface required for matrix_data_value_tag

    template<class Tag>
    static const bool is_valid = matrix_data_value_tag_check_impl<Tag> :: value;
};
                        
//-----------------------------------------------------------------------
//                      scal_data_const_value_tag
//-----------------------------------------------------------------------
// base class for Tags used in creating const_value_scalar
template<class Tag>
struct scal_data_const_value_tag
{
    // Tag must implement:
    // template<class Val>
    // static constexpr Val value();
};

template<class Tag>
struct scal_data_const_value_tag_check_impl
{
    static const bool value = mkd::has_static_member_value2_constexpr<Tag, double>;
};

struct scal_data_const_value_tag_check
{
    // check if Tag has interface required for scal_data_const_value_tag
    template<class Tag>
    static const bool is_valid = scal_data_const_value_tag_check_impl<Tag> :: value;
};

//-----------------------------------------------------------------------
//                      scal_data_value_tag
//-----------------------------------------------------------------------
// base class for Tags used in creating value_scalar
template<class Tag>
struct scal_data_value_tag
{
    // Tag must implement:
    // template<class Val>
    // static Val value();
};

template<class Tag>
struct scal_data_value_tag_check_impl
{
    static const bool value = mkd::has_static_member_value2<Tag, double>;
};

struct scal_data_value_tag_check
{
    // check if Tag has interface required for scal_data_value_tag
    template<class Tag>
    static const bool is_valid = scal_data_value_tag_check_impl<Tag> :: value;
};

//-----------------------------------------------------------------------
//                      scal_data_gen_value_tag
//-----------------------------------------------------------------------
// base class for Tags used in creating gen_scalar
template<class Tag>
struct scal_data_gen_value_tag
{
    // Tag must implement:
    // static void print(std::ostream& os, int prior);
    // 
    // template<class Val, class Local_Storage>
    // static Val eval(const Local_Storage& ls);
};

template<class Tag>
struct scal_data_gen_value_tag_check_impl
{
    static const bool value 
        = mkd::has_static_member_function_print<Tag>
        && mkd::has_static_member_template_function_eval<Tag, double, mkd::local_storage_dummy>;
};

struct scal_data_gen_value_tag_check
{
    // check if Tag has interface required for scal_data_gen_value_tag
    template<class Tag>
    static const bool is_valid = scal_data_gen_value_tag_check_impl<Tag> :: value;
};

}}