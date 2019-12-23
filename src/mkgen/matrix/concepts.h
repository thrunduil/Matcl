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

#include <type_traits>
#include <concepts>

#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace mkgen
{

namespace mk = matcl :: mkgen;

//------------------------------------------------------------------------------
//                      fordward declarations
//------------------------------------------------------------------------------
template<class Tag>
struct matrix_data_const_value_tag;

template<class Tag>
struct matrix_data_value_tag;

template<class T, class V, Integer Row, Integer Col>
concept has_static_member_value = 
    requires { {T::template value<V, Row, Col>()} -> std::convertible_to<V>; };

//------------------------------------------------------------------------------
//                      matrix tags
//------------------------------------------------------------------------------

// concept of Tag used to create const_value_mat
// Tag must be derived from matrix_data_const_value_tag<Tag> and implement:
//    return value at position (Row, Col)
//    template<class V, Integer Row, Integer Col>
//    static constexpr V value();
template<class Tag>
concept Tag_matrix_const_data = std::is_base_of<mk::matrix_data_const_value_tag<Tag>,
                                            Tag>::value
                              && has_static_member_value<Tag, double, 1, 1>;

// concept of Tag used to create value_mat
// Tag must be derived from matrix_data_value_tag<Tag> and implement:
//    return value at position (Row, Col)
//    template<class V, Integer Row, Integer Col>
//    static V value();
template<class Tag>
concept Tag_matrix_data = std::is_base_of<mk::matrix_data_value_tag<Tag>, Tag>::value
                        && has_static_member_value<Tag, double, 1, 1>;

}}

