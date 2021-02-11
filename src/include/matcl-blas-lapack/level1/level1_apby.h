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
//                  				apby
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a + b * Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array and the scalar a
*  TB           - type of scalar b
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y         	- pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  rows         - number of elements
*  a 			- scalar, tests for special cases are not performed
*  b 			- scalar, tests for special cases are not performed
*/
template<class TY, class TB, Integer Rows, 
            bool Use_simd = details::check_simd_scal<TB,TY>::value, class Main = details::true_t>
struct apby
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer rows, const TY& a, const TB& b);
    static void eval(TY* Y, Integer Y_step, Integer rows, const TY& a, const TB& b);
};

//-----------------------------------------------------------------------
//                      apby_test_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a + b * Y
*  test all special cases for a, b (test values 0, +1, -1)
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y array and the scalar a
*  TB           - type of scalar b
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
*  Y            - pointer to Y array
*  Y_ld         - leading dimension of Y matrix
*  rows         - number of rows of X and Y
*  cols         - number of columns of X and Y
*  a            - scalar
*  b            - scalar
*/
template<bool Select, class TY, class TB, Integer Rows, Integer Cols, Integer Continuous, 
         class Main = details::true_t>
struct apby_test_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, Integer rows, Integer cols, const TY& a, const TB& b);
};

//-----------------------------------------------------------------------
//                          apy
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a + Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array and the scalar a
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  rows         - number of elements
*  a            - scalar, tests for special cases are not performed
*/
template<class TY, Integer Rows, 
        bool Use_simd = details::check_simd<TY>::value, class Main = details::true_t>
struct apy
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer rows, const TY& a);
    static void eval(TY* Y, Integer Y_step, Integer rows, const TY& a);
};

//-----------------------------------------------------------------------
//                  				amy
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a - Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array and the scalar a
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y         	- pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  rows         - number of elements
*  a 			- scalar, tests for special cases are not performed
*/
template<class TY, Integer Rows, 
            bool Use_simd = details::check_simd<TY>::value, class Main = details::true_t>
struct amy
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer rows, const TY& a);
    static void eval(TY* Y, Integer Y_step, Integer rows, const TY& a);
};

};};

#include "matcl-blas-lapack/level1/details/level1_apby.inl"
