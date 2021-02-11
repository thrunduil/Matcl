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

#include "matcl-blas-lapack/level1/config.h"
#include "matcl-blas-lapack/level1/utils.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace level1
{

//-----------------------------------------------------------------------
//                  		GENERAL ASSUMPTIONS
//-----------------------------------------------------------------------
// -   different arrays must be disjoint (no aliassing)
// -   negative steps are allowed, i-th element of array Y with stepping 'step'
//         is assumed to be Y[i*step] (not that this is different than in BLAS)
// -   no optimization against scalar arguments of functions unless otherwise
//         states
// -   functions are optimized for step = 1 
// -   if the resultig array stores floating point numbers of type T, then all
//         scalars are converted to floating point scalars of the same precision
//         if this is possible (i.e. the scalar is convertible to the floating
//         point type T)

//-----------------------------------------------------------------------
//                  eval_mat_func_Y
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = func(Y) for every elements in matrix Y
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y matrix
*  Rows         - number of rows of the matrix, 0 if unknown statically
*  Cols         - number of columns of the matrix, 0 if unknown statically
*  Continuous   = 1 - leading dimensions of source and destination matrices
*                   are equal to number of rows
*               = 2 - leading dimension of one of matrices is greater than
*                   number of rows            
*               = 0 - determine dynamically
*  Func         - class with member function implementing:
*                   template<Integer Rows>
*                   void eval(T1* y, Integer rows) const
*                 that evaluates the function func for every element in array y
*                 of static length Rows and dynamic length rows.
*
*  Inputs
*  =======
*  Y            - pointer to Y matrix
*  Y_ld         - leading dimension of Y matrix
*  rows         - number of rows of Y
*  cols         - number of columns of Y
*  func         - instance of the class Func
*/
template<bool Select, class TY, Integer Rows, Integer Cols, Integer Continuous, class Func, 
            class Main = details::true_t>
struct eval_mat_func_Y
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, Integer rows, Integer cols, const Func& fun);
};

//-----------------------------------------------------------------------
//                  eval_mat_func_XY
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = func(X, Y) for every elements in matrices X and Y
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in destination matrix
*  TX           - type of elements in source matrix
*  Rows         - number of rows of the matrix, 0 if unknown statically
*  Cols         - number of columns of the matrix, 0 if unknown statically
*  Continuous   = 1 - leading dimensions of source and destination matrices
*                   are equal to number of rows
*               = 2 - leading dimension of one of matrices is greater than
*                   number of rows            
*               = 0 - determine dynamically
* Func          - class with member function implementing:
*                   template<Integer Rows>
*                   void eval(T1* Y, const T2* X, Integer rows) const
*                 that evaluates the function func for every element in arrays X and Y
*                 of static length Rows and dynamic length rows.
*
*  Inputs
*  =======
*  Y            - pointer to Y matrix
*  Y_ld         - leading dimension of Y matrix
*  X            - pointer to X matrix
*  X_ld         - leading dimension of X matrix
*  rows         - number of rows of X and Y
*  cols         - number of columns of X and Y
*  func         - instance of the class Func
*/
template<bool Select, class TY, class TX, Integer Rows, Integer Cols, Integer Continuous,
        class Func, class Main = details::true_t>
struct eval_mat_func_XY
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, Integer rows, Integer cols,
                     const Func& fun);
};

//-----------------------------------------------------------------------
//                  			set_val
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  Rows         - number of elements of the array, 0 if unknown statically
*
*  Inputs
*  =======
*  Y         	- pointer to Y array
*  Y_step       - optional step between two elements
*  rows         - number of elements
*  a            - scalar
*/
template<class TY, Integer Rows, 
        bool Use_simd = details::check_simd<TY>::value, class Main = details::true_t>
struct set_val
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer rows, const TY& a);
    static void eval(TY* Y, Integer Y_step, Integer rows, const TY& a);
};

//-----------------------------------------------------------------------
//                          set_val_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y array
*  Rows         - number of rows of the matrix, 0 if unknown statically
*  Cols         - number of columns of the matrix, 0 if unknown statically
*  Continuous   = 1 - leading dimensions of source and destination matrices
*                   are equal to number of rows
*               = 2 - leading dimension of one of matrices is greater than
*                   number of rows            
*               = 0 - determine dynamically
*
*  Inputs
*  =======
*  Y         	- pointer to Y array
*  Y_ld         - leading dimension of Y matrix
*  rows         - number of rows of X and Y
*  cols         - number of columns of X and Y
*  a 			- scalar
*/
template<bool Select, class TY, Integer Rows, Integer Cols, Integer Continuous, 
         class Main = details::true_t>
struct set_val_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, Integer rows, Integer cols, const TY& a);
};

//-----------------------------------------------------------------------
//                  copy
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  copy array to destrination from source
*
*  Template arguments
*  =======
*  TY           - type of elements in destination array
*  TX           - type of elements in source array
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  dest         - pointer to destination array
*  dest_step    - optional step between two elements in dest array
*  source       - pointer to source array
*  source_step  - optional step between two elements in source array
*  rows         - number of rows to copy
*/
template<class TY, class TX, Integer Rows, 
        bool Use_simd = details::check_simd<TY,TX>::value, class Main = details::true_t>
struct copy
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* dest, const TX* source, Integer rows);
    static void eval(TY* dest, Integer dest_step, const TX* source, Integer source_step, 
                     Integer rows);
};
    
//-----------------------------------------------------------------------
//                  copy_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  copy matrix to destrination from source
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in destination matrix
*  TX           - type of elements in source matrix
*  Rows         - number of rows of the matrix, 0 if unknown statically
*  Cols         - number of columns of the matrix, 0 if unknown statically
*  Continuous   = 1 - leading dimensions of source and destination matrices
*                   are equal to number of rows
*               = 2 - leading dimension of one of matrices is greater than
*                   number of rows            
*               = 0 - determine dynamically
*
*  Inputs
*  =======
*  dest         - pointer to destination array
*  dest_ld      - leading dimension of destination matrix
*  source       - pointer to source array
*  source_ld    - leading dimension of source matrix
*  rows         - number of rows to copy
*  cols         - number of columns to copy
*/
template<bool Select, class TY, class TX, Integer Rows, Integer Cols, Integer Continuous, 
    bool Use_simd = details::check_simd<TY,TX>::value, class Main = details::true_t>
struct copy_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* dest, Integer dest_ld, const TX* source, Integer source_ld, 
                     Integer rows, Integer cols);
};

//-----------------------------------------------------------------------
//                  				swap
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  swap content of two vectors
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y_1, Y_2 arrays
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y1         	- pointer to Y_1 array
*  Y1_step      - optional step between two elements in Y1 array
*  Y2         	- pointer to Y_2 array
*  Y2_step      - optional step between two elements in Y2 array
*  rows         - number of elements
*/
template<class TY, Integer Rows, 
        bool Use_simd = details::check_simd<TY>::value, class Main = details::true_t>
struct swap
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y1, TY* Y2, Integer rows);
    static void eval(TY* Y1, Integer Y1_step, TY* Y2, Integer Y2_step, Integer rows);
};

};};


#include "matcl-blas-lapack/level1/details/level1_basic.inl"
