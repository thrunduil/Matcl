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

#include "matcl-core/matrix/scalar_types.h"

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

template<class T, class V, Integer Row, Integer Col>
concept has_static_member_value = 
    requires { {T::template value<V, Row, Col>()} -> std::convertible_to<V>; };

template<class T, Integer Row, Integer Col>
concept has_static_member_get_element_impl = 
    requires { typename T :: template get_element_impl<Row, Col>; };

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
    static const bool value = mkd::has_static_member_value<Tag, double, 1, 1>;
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
                        
}}