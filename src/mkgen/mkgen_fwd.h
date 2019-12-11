/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2019
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

namespace matcl { namespace mkgen 
{

using Integer   = matcl::Integer;

//------------------------------------------------------------------------
//                              arrays
//------------------------------------------------------------------------
//data arrays
template<class Tag>                                     struct array;
template<class Tag>                                     struct output_array;
template<class Tag>                                     struct const_array;
template<class Tag, Integer Rows, Integer Cols>         struct temp_output_array;
template<class Tag, Integer Mat_Rows, Integer Mat_Cols,bool Force> 
                                                        struct mat_temp_array;

//computation arrays
template<Integer K, class Array1, class Array2>         struct mat_mult_array;
template<class Array1, class Scalar2>                   struct mat_scal_mult_array;
template<class Array1, class Array2>                    struct mult_rows_array;
template<class Array1, class Array2>                    struct mult_cols_array;
template<class Array1, class Array2>                    struct mult_array;
template<class Array1, class Array2>                    struct div_array;
template<class Array1, class Scal2>                     struct div_array_mat_scal;
template<class Array2, class Scal1>                     struct div_array_scal_mat;
template<Integer M, Integer N, class Array>             struct mat_trans_array;
template<Integer M, Integer N, class Array>             struct mat_ctrans_array;
template<class Tag, Integer M, Integer N, class Array>  struct mat_ufunc_array;
template<class Tag, Integer M, Integer N, class Array1, class Array2>
                                                        struct mat_bfunc_array;
template<class Tag, Integer M, Integer N, class Array1, class Array2>
                                                        struct mat_scal_bfunc_array;
template<class Tag, Integer M, Integer N, class Array1, class Array2>
                                                        struct scal_mat_bfunc_array;
template<Integer M,Integer N,class Array1,class Array2> struct mat_plus_array;
template<Integer M,Integer N,class Array1,class Array2> struct mat_scal_plus_array;
template<Integer M,Integer N,class Array1,class Array2> struct mat_minus_array;
template<Integer M,Integer N,class Array1,class Array2> struct mat_scal_minus_array;
template<Integer M,Integer N,class Array1,class Array2> struct scal_mat_minus_array;
template<Integer M, Integer N, class Array>             struct mat_uminus_array;

//assignments
template<class Tag, class... Assign_List>               struct virtual_array;
template<Integer M,Integer N,class Array1,class Array2> struct mat_assign_array;
template<Integer M,Integer N,class Array1,class Array2> struct mat_scal_assign_array;
template<Integer M,Integer N,class Array1,class Colon, Integer M2, Integer N2, class Array2> 
                                                        struct mat_assign_array_colon;
template<class Array_t, Integer Offset, Integer Step>   struct sub_array_1;
template<class Array_t, Integer Offset1, Integer Offset2, Integer Step1, Integer step2>
                                                        struct sub_array_2;

//technical arrays
template<class Ret_Tag>                                 struct empty_array;
                                                        struct call_array_type;

//------------------------------------------------------------------------
//                              other
//------------------------------------------------------------------------
template<Integer M, Integer N, class Array_t, class Deps>
struct ct_matrix;

template<class ...Items>
struct list;

template<class Tag, class Mat, class Assignments = list<>>
struct computation;

template<Integer Pos> struct colon;
struct colon_all;

template<Integer Start, Integer End> 
struct colon2;

template<Integer Start, Integer Step, Integer End> 
struct colon3;

template<class Mat1, class Mat2>
struct mat_assign;

template<class Subject, class Scalar, class Colon>
struct comp_assign_1;

template<class Comp>
struct make_comp_result;

template<class Mat_L, class Mat_R, class Colon_1>
struct mat_virtual_assign_1;

template<class Mat_L, class Tag, bool Force>
struct mat_temporary;

template<class T1, class T2>
struct expr_plus;

template<class T1, class T2>
struct expr_minus;

template<class ... T>
struct expr_mult;

template<class List_1, class List_2>
struct expr_dot;

template<class T, bool With_Forced>
struct is_temporary_mat;

template<class Temp_Tag, class Ret_Tag, class Colon, Integer Rows, Integer Cols, bool Init>
struct dps_modif;

template<class Tag, class Colon, bool Init>
struct modif;

template<class Tag, class Colon, Integer Mat_Rows, Integer Mat_Cols>
struct modif2;

template <class Tag, Integer Mat_Rows, Integer Mat_Cols, Integer Row, Integer Col>
struct get_temporary;

template<class Elem, Integer Step>
struct element_step;

template<class Array, Integer Row, Integer Col>
struct get_array_elem;

template<class Val, class Elem>
struct loop_context_data_scalar;

}}