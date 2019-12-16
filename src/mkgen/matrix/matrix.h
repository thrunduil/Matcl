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

#include <iosfwd>

#include "mkgen/matrix/scalar.h"
#include "mkgen/details/matrix/matrix_checks.h"

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
// produce an error. For a matrix with generic elements Array_t type is basically array<Tag>
// with unique Tag for each matrix. 
// Deps is a type representing dependencies from runtime values; must be specialization
// of dps type
template<Integer M, Integer N, class Array_t, class Deps>
struct ct_matrix
{
    public:
        // check arguments
        using check1    = typename mkd::check_valid_matrix_array<Array_t>::type;
        using check2    = typename mkd::check_deps<Deps>::type;

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
        static auto elem(colon<Pos>)        -> ct_scalar<mkd::scalar_mat_elem_1<ct_matrix, Pos>, Deps>;

        // get element at row Row and column Col
        template<Integer Row, Integer Col>
        static auto elem(colon<Row>, colon<Col>) 
                                            -> ct_scalar<mkd::scalar_mat_elem_2<ct_matrix, Row, Col>, Deps>;

        // get submatrix Mat(Colon_1)
        template<class Colon_1>
        static auto sub(Colon_1)            -> typename mkd::submatrix_maker_1<ct_matrix, Colon_1>::type;

        // get submatrix Mat(Colon_1, Colon_2)
        template<class Colon_1, class Colon_2>
        static auto sub(Colon_1, Colon_2)   -> typename mkd::submatrix_maker_2<ct_matrix, 
                                                    Colon_1, Colon_2>::type;

        // make assignment This(Colon_1) = Mat; 
        // this type must be a virtual_matrix
        template<class Colon_1, class Mat>
        static auto assign_1(Colon_1, Mat)  -> typename mkd::mat_virtual_assign_1<ct_matrix, Mat, Colon_1>::type;

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
//TODO:

//stores statically known data, accessible through tag argument
template<Integer M, Integer N, class Tag>
using const_mat = ct_matrix<M,N, mkd::const_array<Tag>,empty_deps>;

//stores generic data unknown statically, usually supplied by some data_provider
template<Integer M, Integer N, class Tag>
using gen_mat = ct_matrix<M,N, mkd::gen_array<Tag>, extern_deps<Tag>>;

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
