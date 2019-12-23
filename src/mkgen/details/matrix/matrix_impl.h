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

#include "mkgen/matrix/matrix.h"
#include "mkgen/details/matrix/matrix_arrays.h"
#include "mkgen/details/matrix/colon_func.h"

namespace matcl { namespace mkgen { namespace details
{

//------------------------------------------------------------------------------
//                      forward declarations
//------------------------------------------------------------------------------
template<class Array>
struct is_virtual_matrix_array;

//------------------------------------------------------------------------------
//                      submatrix_maker_1
//------------------------------------------------------------------------------
// implements ct_matrix<>::sub(Colon_1)
template<class A, class Colon_1>
struct submatrix_maker_1
{
    static const bool is_mat    = is_matrix<A>::value;
    static_assert(is_mat == true, "A must be ct_matrix");

    static const Integer M      = A::rows;
    static const Integer N      = A::cols;
    using array_t               = typename A::array_type;
    using deps                  = typename A::dps_type;

    using dum                   = typename colon_func::check_colon<Colon_1, M * N>::type;

    static const Integer size   = colon_func::size<Colon_1, M>::value;
    static const Integer offset = colon_func::offset<Colon_1>::value;
    static const Integer step   = colon_func::step<Colon_1>::value;

    using new_array = sub_array_1<array_t, offset, step>;
    using type      = ct_matrix<size,1, new_array, deps>;
};

//------------------------------------------------------------------------------
//                      submatrix_maker_2
//------------------------------------------------------------------------------
template<class A, class Colon_1, class Colon_2>
struct submatrix_maker_2
{
    static const bool is_mat    = is_matrix<A>::value;
    static_assert(is_mat == true, "A must be ct_matrix");

    static const Integer M          = A::rows;
    static const Integer N          = A::cols;
    using array_t                   = typename A::array_type;
    using deps                      = typename A::dps_type;

    using dum1                      = typename colon_func::check_colon<Colon_1, M>::type;
    using dum2                      = typename colon_func::check_colon<Colon_2, N>::type;

    static const Integer size1      = colon_func::size<Colon_1, M>::value;
    static const Integer size2      = colon_func::size<Colon_2, N>::value;
    static const Integer offset1    = colon_func::offset<Colon_1>::value;
    static const Integer offset2    = colon_func::offset<Colon_2>::value;
    static const Integer step1      = colon_func::step<Colon_1>::value;
    static const Integer step2      = colon_func::step<Colon_2>::value;

    using new_array = mkd::sub_array_2<array_t, offset1, offset2, step1, step2>;
    using type      = ct_matrix<size1, size2, new_array, deps>;
};

//------------------------------------------------------------------------------
//                      submatrix_elem_1
//------------------------------------------------------------------------------
// implements ct_matrix<>:: elem(colon<Pos>)
template<class A, Integer Pos>
struct submatrix_elem_1
{
    static const bool is_mat    = is_matrix<A>::value;
    static_assert(is_mat == true, "A must be ct_matrix");

    static const Integer rows = A :: rows;
    static const Integer col = (Pos-1)/rows + 1;
    static const Integer row = Pos - (col-1) * rows;

    using array_t       = typename A :: array_type;
    using deps_t        = typename A :: dps_type;

    using elem_type     = typename array_t::template get_element<row, col>::type;
    using type          = ct_scalar<elem_type, deps_t>;
};

//------------------------------------------------------------------------------
//                      submatrix_elem_2
//------------------------------------------------------------------------------
// implements ct_matrix<>:: elem(colon<Row>, colon<Col>) 
template<class A, Integer Row, Integer Col>
struct submatrix_elem_2
{
    static const bool is_mat    = is_matrix<A>::value;
    static_assert(is_mat == true, "A must be ct_matrix");

    using array_t       = typename A :: array_type;
    using deps_t        = typename A :: dps_type;

    using elem_type     = typename array_t::template get_element<Row, Col>::type;
    using type          = ct_scalar<elem_type, deps_t>;
};

//------------------------------------------------------------------------------
//                      is_virtual_matrix_array
//------------------------------------------------------------------------------

// return true if Array is virtual_array<Tag, Assigns...>
template<class Array>
struct is_virtual_matrix_array
{
    static const bool value = false;
};

template<class Tag, class... Assign_List>
struct is_virtual_matrix_array<mkd::virtual_array<Tag, Assign_List...>>
{
    static const bool value = true;
};

}}}

namespace matcl { namespace mkgen
{

//------------------------------------------------------------------------------
//                      isa functions
//------------------------------------------------------------------------------

// return true if T is ct_matrix
template<class T>
struct mkgen::is_matrix
{
    static const bool value = false;
};

template<Integer M, Integer N, Mat_array Array_t, class Deps>
struct mkgen::is_matrix<ct_matrix<M, N, Array_t, Deps>>
{
    static const bool value = true;
};

// return true if T is virtual ct_matrix, i.e. has type virtual_mat<M, N, Tag>
template<class T>
struct mkgen::is_virtual_matrix
{
    static const bool value = false;
};

template<Integer M, Integer N, Mat_array Array_t, class Deps>
struct mkgen::is_virtual_matrix<ct_matrix<M, N, Array_t, Deps>>
{
    static const bool value = mkd::is_virtual_matrix_array<Array_t>::value;
};

}}