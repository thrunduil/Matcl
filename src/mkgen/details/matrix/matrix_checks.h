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

template<class T>
struct make_true_type
{
    using type = std::true_type;
};

template<class C>
struct has_template_alias
{                                                                          
    /* attempt to call it */         
    template<typename T>                                                   
    static constexpr auto check(T*)                                        
        -> typename make_true_type<typename T::template get_element_impl<1,1>>::type;
                                                                           
    template<typename>                                                     
    static constexpr std::false_type check(...);                           
                                                                           
    using type  = decltype(check<C>(0));                                   
                                                                           
    static const bool value     = type::value;                             
};

//-----------------------------------------------------------------------
//                      check_valid_matrix_array
//-----------------------------------------------------------------------
// check if Array_t parameter supplied to ct_matrix is valid
template<class Array_t>
struct check_valid_matrix_array
{
    static const bool is_array  = std::is_base_of<mkd::matrix_array<Array_t>, Array_t>::value;

    static_assert(is_array == true, "Array_t is not matrix_array<>");

    using type  = typename Array_t::template check_matrix_array<void>;
};

// check if Array_t parameter supplied to matrix_array is valid
template<class Array, class Ret>
struct check_matrix_array_impl
{
    // Array must implement:

    // template<Integer Row, Integer Col>
    // using get_element_impl   = [impl]

    static const bool has_get   = has_template_alias<Array>::value;

    static_assert(has_get == true, "Array must implement template alias get_element_impl");

    using type = Ret;

    // type T = get_element_impl<Row,Col>::type will be checked, when fully constructed
};


}}}
