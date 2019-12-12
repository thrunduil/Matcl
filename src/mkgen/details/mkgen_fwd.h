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
#include "matcl-core/matrix/scalar_types.h"
#include "mkgen/matrix/scalar.h"

namespace matcl { namespace mkgen
{

// compile time scalar value
template<class Data, class Deps> 
class ct_scalar;

// compile time matrix
template<Integer M, Integer N, class Array_t, class Deps>
struct ct_matrix;

// compute result of ct_scalar::compute()
template<class Data, class Deps, class Tag>
struct make_evaled_scalar;

// get element from an array
template<class Array, Integer Row, Integer Col>
struct get_array_elem;

//TODO:
template<class Mat_L, class Mat_R, class Colon_1>
struct mat_virtual_assign_1;

//TODO
template<class Mat_L, class Tag, bool Force>
struct mat_temporary;

}}

namespace matcl { namespace mkgen { namespace details
{
    
// represents a scalar storing rational values
template<Integer N, Integer D>
struct scal_data_rational;

// represents a scalar storing values of type Value_type
template<class Tag, class Value_type>
struct scal_data_value;

// append to Arr_List all arrays required by this scalar; 
// implements ct_scalar::get_arrays
template<class Data, class Deps, Integer Step, class Arr_List>
struct get_arrays_scalar;

template<class Elem, Integer Step, class Type>
struct array_item;

// implements ct_scalar::eval_loop
template<class Loop_Storage, class Data, class Deps>
struct eval_loop_scalar;

// store data in ct_scalar
template<class Data>
struct scalar_data;

// store data in ct_matrix
template<class Array>
struct matrix_array;

// make submatrix
template<class Mat, class Colon_1, class Colon_2>
struct submatrix_maker_2;

// make submatrix
template<class Mat, class Colon_1>
struct submatrix_maker_1;

// const_mat array
template<class Tag>                                     
struct const_array;

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

// dep check
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