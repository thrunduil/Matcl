/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/details/mpl.h"
#include "matcl-core/matrix/enums.h"
#include "matcl-core/details/isa.h"

namespace matcl { namespace details
{


template<class T>	struct is_matrix    {	static const bool value = false; };
template<class T>	struct is_sparse    {	static const bool value = false; };
template<class T>	struct is_band      {	static const bool value = false; };

template<class value_type,class struct_type>
struct is_matrix<raw::Matrix<value_type,struct_type>>
{	
    static const bool value = true;	 
};

template<class value_type>
struct is_matrix<raw::dense_matrix_base<value_type>>
{	
    static const bool value = true;	 
};

template<class value_type>
struct is_matrix<raw::sparse_matrix_base<value_type>>
{	
    static const bool value = true;	 
};

template<class value_type>
struct is_sparse<raw::Matrix<value_type,struct_sparse>>
{	
    static const bool value = true;	 
};

template<class value_type>
struct is_sparse<raw::sparse_matrix_base<value_type>>
{	
    static const bool value = true;	 
};

template<class value_type>
struct is_band<raw::Matrix<value_type,struct_banded>>
{	
    static const bool value = true;	 
};

template<class T>
struct is_typed_matrix
{ 
    static const bool value = false; 
};

template<class T, bool Is_safe>
struct is_typed_matrix<dense_matrix<T,Is_safe>>
{ 
    static const bool value = true; 
};

template<class T, bool Is_safe>
struct is_typed_matrix<sparse_matrix<T,Is_safe>>
{ 
    static const bool value = true; 
};

template<class T, bool Is_safe>
struct is_typed_matrix<band_matrix<T,Is_safe>>
{ 
    static const bool value = true; 
};

template<class T> struct is_object_matrix                       { static const bool value = false; };
template<class T> struct is_object_matrix<object_matrix<T>>     { static const bool value = true; };

template<class T> struct is_submatrix                           { static const bool value = false; };
template<> struct is_submatrix<sub_matrix>                      { static const bool value = true; };
template<> struct is_submatrix<sub_matrix_1>                    { static const bool value = true; };
template<> struct is_submatrix<sub_matrix_2>                    { static const bool value = true; };

template<class T>
struct is_typed_submatrix                                       { static const bool value = false; };
template<class T>
struct is_typed_submatrix<sub_dense_matrix<T>>                  { static const bool value = true; };
template<class T>
struct is_typed_submatrix<sub_sparse_matrix<T>>                 { static const bool value = true; };
template<class T>
struct is_typed_submatrix<sub_sparse_matrix_1<T>>               { static const bool value = true; };
template<class T>
struct is_typed_submatrix<sub_sparse_matrix_2<T>>               { static const bool value = true; };
template<class T>
struct is_typed_submatrix<sub_band_matrix<T>>                   { static const bool value = true; };
template<class T>
struct is_typed_submatrix<sub_band_matrix_1<T>>                 { static const bool value = true; };
template<class T>
struct is_typed_submatrix<sub_band_matrix_2<T>>                 { static const bool value = true; };

template<class T>
struct is_object_submatrix                                      { static const bool value = false; };
template<class T>
struct is_object_submatrix<sub_object_matrix<T>>                { static const bool value = true; };
template<class T>
struct is_object_submatrix<sub_object_matrix_1<T>>              { static const bool value = true; };
template<class T>
struct is_object_submatrix<sub_object_matrix_2<T>>              { static const bool value = true; };

template<class T> struct is_concat                              { static const bool value = false; };
template<> struct is_concat<mat_row>                            { static const bool value = true; };
template<> struct is_concat<mat_col>                            { static const bool value = true; };

template<class T> struct is_typed_concat                        { static const bool value = false; };
template<class T> struct is_typed_concat<dense_row<T>>          { static const bool value = true; };
template<class T> struct is_typed_concat<dense_col<T>>          { static const bool value = true; };
template<class T> struct is_typed_concat<sparse_row<T>>         { static const bool value = true; };
template<class T> struct is_typed_concat<sparse_col<T>>         { static const bool value = true; };

template<class T> struct is_object_concat                       { static const bool value = false; };
template<class T> struct is_object_concat<object_row<T>>        { static const bool value = true; };
template<class T> struct is_object_concat<object_col<T>>        { static const bool value = true; };

template<class T>
struct is_convertible_to_matrix
{
    static const bool value  
                        = is_scalar<T>::value || is_typed_matrix<T>::value 
                        || is_object_matrix<T>::value || is_submatrix<T>::value 
                        || is_typed_submatrix<T>::value || is_object_submatrix<T>::value
                        || is_concat<T>::value || is_typed_concat<T>::value 
                        || is_object_concat<T>::value || std::is_same<T, unique_matrix>::value;
};

};};