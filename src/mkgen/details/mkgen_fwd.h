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

#include "mkgen/mkgen_config.h"
#include "mkgen/matrix/concepts.h"
#include "mkgen/matrix/scalar.h"

namespace matcl { namespace mkgen
{

// compile time scalar value
template<Scal_data Data, DPS Deps> 
class ct_scalar;

// compile time matrix
template<Integer M, Integer N, Mat_array Array_t, DPS Deps>
struct ct_matrix;

// compute result of ct_scalar::compute()
template<Scal_data Data, DPS Deps, Tag_comp Tag>
struct make_evaled_scalar;

//TODO
template<class Mat_L, class Tag, bool Force>
struct mat_temporary;

// TODO:
template<class Array, Integer Row, Integer Col>
struct sub_array_2_get_elem;

//TODO
template<class Array, Integer Row, Integer Col>
struct mat_ufunc_array_get_elem;

//TODO
template<class Array, Integer Row, Integer Col>
struct sub_array_1_get_elem;

//TODO:
template <class Tag, Integer Mat_Rows, Integer Mat_Cols, Integer Row, Integer Col>
struct get_temporary;

//TODO
template<class Array, Integer Row, Integer Col>
struct mat_assign_array_colon_get_elem;

//TODO
template<class Array, Integer Row, Integer Col>
struct mat_assign_array_get_elem;

//TODO
template<class Array, Integer Row, Integer Col>
struct mat_scal_assign_array_get_elem;

//TODO
template<class Array, Integer Row, Integer Col>
struct mat_bfunc_array_get_elem;

//TODO
template<class Array, Integer Row, Integer Col>
struct mat_scal_bfunc_array_get_elem;

//TODO
template<class Array, Integer Row, Integer Col>
struct scal_mat_bfunc_array_get_elem;

}}

namespace matcl { namespace mkgen { namespace details
{
    
// represents a scalar storing rational values N / D
template<Integer N, Integer D>
struct scal_data_rational;

// represents a scalar storing values of type Value_type
// known at compile time
template<Tag_scalar_const_value Tag, class Value_type>
struct scal_data_const_value;

// represents a scalar storing values of type Value_type
template<Tag_scalar_value Tag, class Value_type>
struct scal_data_value;

// represents a scalar storing external values
template<Tag_scalar_gen_value Tag>
struct scal_data_gen_value;

// append to Arr_List all arrays required by this scalar; 
// implements ct_scalar::get_arrays
template<Scal_data Data, DPS Deps, Integer Step, class Arr_List>
struct get_arrays_scalar;

template<class Elem, Integer Step, class Type>
struct array_item;

// implements ct_scalar::eval_loop
template<class Loop_Storage, class Data>
struct eval_loop_scalar;

// make submatrix
template<class Mat, class Colon_1, class Colon_2>
struct submatrix_maker_2;

// make submatrix
template<class Mat, class Colon_1>
struct submatrix_maker_1;

// get element from a matrix
template<class Mat, Integer Pos>
struct submatrix_elem_1;

// get element from a matrix
template<class Mat, Integer Row, Integer Col>
struct submatrix_elem_2;

// value_mat array
template<Tag_matrix_data Tag, class Value_type>                                     
struct matrix_array_value;

// const_value_mat array
template<Tag_matrix_const_data Tag, class Value_type>                                     
struct matrix_array_const_value;

// gen_mat array
template<class Tag>                                     
struct gen_array;

// output_mat array
template<class Tag> 
struct output_array;

// temp_output_mat array
template<class Tag, Integer Rows, Integer Cols>
struct temp_output_array;

// virtual_mat array
template<class Tag, class... Assign_List>
struct virtual_array;

// make virtual assignment Mat_L(Colon_1) = Mat_R
template<class Mat_L, class Mat_R, class Colon_1>
struct mat_virtual_assign_1;

// implements virtual_array::get_element
template<Integer Row, Integer Col, class... Items>
struct get_virtual_array_assignment;

// implements matrix or scalar addition
template<class Mat1, class Mat2>
struct mat_plus_impl;

// implements matrix or scalar substraction
template<class Mat1, class Mat2>
struct mat_minus_impl;

// implements matrix or scalar multiplication
template<class Mat, class Mat2>
struct mat_mult_impl;

// implements matrix or scalar element by element division
template<class Mat, class Mat2>
struct mat_div_impl;

// implements mult_rows function
template<class Mat, class D>
struct mult_rows_impl;

// implements mult_cols function
template<class Mat, class D>
struct mult_cols_impl;

// implements mul function
template<class Mat, class D>
struct mul_impl;

// representation of mult expression
template<bool Flag, class S, class ... T>
struct expr_mult_sd;

// representation of plus expression
template<bool Flag, class S, class ... T>
struct expr_plus_sd;

// checks
template<class Temp_Tag, Integer Size, class Type>
struct check_dep_argument;

template<class... T>
struct check_dps_argument;

template<class Tag>
struct check_valid_dep_tag;

template<Integer Size>
struct check_valid_dep_size;

template<class Type>
struct check_valid_dep_type;

}}}