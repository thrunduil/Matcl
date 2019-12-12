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

namespace matcl { namespace mkgen { namespace details
{

//-----------------------------------------------------------------------
//                      matrix_array
//-----------------------------------------------------------------------
// base class for ct_matrix arrays
template<class Data>
struct matrix_array
{
    //TODO
    //check arguments
    template<class Dummy>
    using check     = Dummy;

    //TODO
    /*
    //check arguments
    template<class Dummy>
    using check     = typename details::check_scalar_data_impl<Data, Dummy>::type;
    */
};

//-----------------------------------------------------------------------
//                      arrays
//-----------------------------------------------------------------------

template<class Array_t, Integer Offset, Integer Step>
struct sub_array_1 : public matrix_array<sub_array_1<Array_t, Offset, Step>>
{};

template<class Array_t, Integer Offset1, Integer Offset2, Integer Step1, Integer Step2>
struct sub_array_2 : public matrix_array<sub_array_2<Array_t, Offset1, Offset2, Step1, Step2>>
{};

// const_mat array
template<class Tag>                                     
struct const_array : public matrix_array<const_array<Tag>>
{};

template<class Tag>                                     
struct gen_array : public matrix_array<gen_array<Tag>>
{};

template<class Tag> 
struct output_array : public matrix_array<output_array<Tag>>
{};

template<class Tag, Integer Rows, Integer Cols>
struct temp_output_array : public matrix_array<temp_output_array<Tag, Rows, Cols>>
{};

template<class Tag, class... Assign_List>
struct virtual_array : public matrix_array<virtual_array<Tag, Assign_List ...>>
{};

template<class Tag, Integer M, Integer N, class Array>
struct mat_ufunc_array : public matrix_array<mat_ufunc_array<Tag, M, N, Array>>
{};

template<Integer K, class Array1, class Array2>
struct mat_mult_array : public matrix_array<mat_mult_array<K, Array1, Array2>>
{};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_scal_plus_array : public matrix_array<mat_scal_plus_array<M, N, Array1, Array2>>
{};

template<class Tag, Integer Mat_Rows, Integer Mat_Cols, bool Force>
struct mat_temp_array : public matrix_array<mat_temp_array<Tag, Mat_Rows, Mat_Cols, Force>>
{};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_scal_minus_array : public matrix_array<mat_scal_minus_array<M, N, Array1, Array2>>
{};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_plus_array : public matrix_array<mat_plus_array<M, N, Array1, Array2>>
{};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_minus_array : public matrix_array<mat_minus_array<M, N, Array1, Array2>>
{};

template<class Array1, class Array2>
struct mult_rows_array : public matrix_array<mult_rows_array<Array1, Array2>>
{};

template<class Ret_Tag>
struct empty_array : public matrix_array<empty_array<Ret_Tag>>
{};

template<Integer M, Integer N, class Array1, class Array2>
struct mat_assign_array : public matrix_array<mat_assign_array<M, N, Array1, Array2>>
{};

template<Integer M,Integer N,class Array1,class Colon, Integer M2, Integer N2, class Array2> 
struct mat_assign_array_colon : public matrix_array<
                                mat_assign_array_colon<M, N, Array1, Colon, M2, N2, Array2>>
{};

//TODO
/*
// Symbolic Array of generic elements of type element<Tag, row, col, ...> for all 
// pairs of (row, col) available in given matrix.
template<class Tag>
struct array {};

template<class Tag>
struct output_array {};
*/

}}}
