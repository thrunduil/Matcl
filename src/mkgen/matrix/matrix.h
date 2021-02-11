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

#include <iosfwd>

#include "mkgen/matrix/scalar.h"
#include "mkgen/matrix/concepts.h"

namespace matcl { namespace mkgen
{

// access matrix element at position Pos (1-based)
template<Integer Pos>
struct colon{};

// access all elements
struct colon_all{};

// access matrix elements from index Start to index End (1-based),
// i.e. in Matlab's notation Start : End
template<Integer Start, Integer End>
struct colon2{};

// access matrix elements from index Start to index End with step Step
// (1-based), i.e. in Matlab's notation  Start : Step : End
template<Integer Start, Integer Step, Integer End>
struct colon3{};

// compile time matrix with size M x N storing symbolic elements in array Array_t.
//
// Matrix 1x1 is not a scalar, for example multiplying 1x1 matrix and 2x2 matrix will
// produce an error. Array_t must be derived from matrix_array
// Deps is a type representing dependencies from runtime values; must be specialization
// of dps type
template<Integer M, Integer N, Mat_array Array_t, DPS Deps>
struct ct_matrix
{
    public:
        // number of rows
        static const Integer    rows        = M;

        // number of columns
        static const Integer    cols        = N;

        // number of elements
        static const Integer    size        = M * N;

        // type storing elements of this matrix
        using array_type                    = Array_t;

        // type representing dependencies from runtime values
        using dps_type                      = Deps;

    public:
        // get element at position Pos as a scalar
        template<Integer Pos>
        static auto elem(colon<Pos>)        -> typename mkd::submatrix_elem_1<ct_matrix, Pos>::type;

        // get element at row Row and column Col
        template<Integer Row, Integer Col>
        static auto elem(colon<Row>, colon<Col>) 
                                            -> typename mkd::submatrix_elem_2<ct_matrix, Row, Col>::type;

        // get submatrix Mat(Col_t)
        template<Colon Col_t>
        static auto sub(Col_t)              -> typename mkd::submatrix_maker_1<ct_matrix, Col_t>::type;

        // get submatrix Mat(Col_1, Col_2)
        template<Colon Col_1, Colon Col_2>
        static auto sub(Col_1, Col_2)       -> typename mkd::submatrix_maker_2<ct_matrix, Col_1, Col_2>::type;

        // make assignment This(Col_t) = Mat; 
        // this type must be a virtual_matrix
        template<Colon Col_t, Mat_or_scalar Mat>
        static auto assign_1(Col_t, Mat)    -> typename mkd::mat_virtual_assign_1<ct_matrix, Mat, Col_t>::type;

        // store current results in temporary matrix (placed on the stack) if cond == true,
        // otherwise return this matrix.
        template<class Tag, bool Force = false>
        static auto make_temp()             -> typename mat_temporary<ct_matrix, Tag, Force>::type;

        // print matrix; add nspaces white spaces at the beginning
        template<class Subs_Context>
        static void print(std::ostream& os, int nspaces = 0);

        //TODO: add compute function
};

//------------------------------------------------------------------------------
//                      predefined matrices
//------------------------------------------------------------------------------

// matrix storing a value of type Value_type defined by the tag Tag; value
// cannot depend on external data; Tag must be derived from matrix_data_const_value_tag
// Tag::get_elem<Row,Col> must evaluate at compile time
template<Integer M, Integer N, Tag_matrix_cdata Tag, Value Val_t>
using const_value_mat = ct_matrix<M, N, mkd::matrix_array_const_value<Tag, Val_t>, empty_deps>;

// matrix storing a values of type Value_type defined by the tag Tag; values
// cannot depend on external data; Tag must be derived from matrix_data_value_tag
template<Integer M, Integer N, Tag_matrix_data Tag, Value Val_t>
using value_mat = ct_matrix<M, N, mkd::matrix_array_value<Tag, Val_t>, empty_deps>;

//TODO:

//stores generic data unknown statically, usually supplied by some data_provider
template<Integer M, Integer N, class Tag>
using gen_mat = ct_matrix<M, N, mkd::gen_array<Tag>, extern_deps<Tag>>;

//stores results
template<Integer M, Integer N, class Tag>
using output_mat = ct_matrix<M,N, mkd::output_array<Tag>, extern_deps<Tag>>;

template<Integer M, Integer N, class Tag>
using temp_output_mat = ct_matrix<M,N, mkd::temp_output_array<Tag,M,N>, return_deps<Tag,M,N>>;

//matrix for storing local computation results, without creating
//runtime buffers.
template<Integer M, Integer N, class Tag>
using virtual_mat = ct_matrix<M,N, mkd::virtual_array<Tag>,empty_deps>;

//------------------------------------------------------------------------------
//                      isa functions
//------------------------------------------------------------------------------
// return true if T is ct_matrix
template<class T>
struct is_matrix;

// return true if T is virtual ct_matrix, i.e. has type virtual_mat<M, N, Tag>
template<class T>
struct is_virtual_matrix;

}}

#include "mkgen/details/matrix/matrix_impl.h"
