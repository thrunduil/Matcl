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
//                  				xy
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Z = X*Y
*
*  Template arguments
*  =======
*  TZ           - type of elements in Z array
*  TX           - type of elements in X array
*  TY           - type of elements in Y array
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Z            - pointer to Z array
*  Z_step       - optional step between two elements in Z array
*  X         	- pointer to X array
*  X_step       - optional step between two elements in X array
*  Y         	- pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  rows         - number of elements
*/
template<class TZ, class TX, class TY, Integer Rows, 
        bool Use_simd = details::check_simd<TZ,TX,TY>::value, class Main = details::true_t>
struct xy
{
    static_assert(md::dependent_false<TZ>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TZ* Z, const TX* X, const TY* Y, Integer rows);
    static void eval(TZ* Z, Integer Z_step, const TX* X, Integer X_step, const TY* Y, 
                     Integer Y_step, Integer rows);
};

//-----------------------------------------------------------------------
//                  				axypz
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Z = a * X * Y + Z
*
*  Template arguments
*  =======
*  TZ           - type of elements in Z array
*  TX           - type of elements in X array
*  TY           - type of elements in Y array
*  TA           - type of scalar a
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Z         	- pointer to Z array
*  Z_step       - optional step between two elements in Z array
*  X         	- pointer to X array
*  X_step       - optional step between two elements in X array
*  Y       		- pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  rows         - number of elements
*  a 			- scalar, tests for special cases are not performed
*/
template<class TZ, class TX, class TY, class TA, Integer Rows, 
        bool Use_simd = details::check_simd_scal<TA, TZ,TX,TY>::value, class Main = details::true_t>
struct axypz
{
    static_assert(md::dependent_false<TZ>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TZ* Z, const TX* X, const TY* Y, Integer rows, const TA& a);
    static void eval(TZ* Z, Integer Z_step, const TX* X, Integer X_step, const TY* Y, 
                     Integer Y_step, Integer rows, const TA& a);
};

//-----------------------------------------------------------------------
//                          dot
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form val = sum X[i] * Y[i]
*
*  Template arguments
*  =======
*  T1           - type of elements in X and Y arrays
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*/
template<class T1, Integer Rows, 
            bool Use_simd = details::check_simd<T1>::value, class Main = details::true_t>
struct dot
{
    static_assert(md::dependent_false<T1>::value,
                    "this version not implemented, probably one of argument is invalid");

    static T1 eval(const T1* Y, const T1* X, Integer rows);
    static T1 eval(const T1* Y, Integer Y_step, const T1* X, Integer X_step, Integer rows);
};

//-----------------------------------------------------------------------
//                  				x_abs
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = abs(x)
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X         	- pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*/
template<class TY, class TX, Integer Rows, 
        bool Use_simd = details::check_simd<TY,TX>::value && details::is_complex<TY>::value == false, 
        class Main = details::true_t>
struct apx_abs
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, Integer rows);
};

//-----------------------------------------------------------------------
//                  				apx_abs
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a + abs(x)
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array and the scalar a
*  TX           - type of elements in X array
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X         	- pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  a         	- scalar
*/
template<class TY, class TX, Integer Rows, 
        bool Use_simd = details::check_simd<TY,TX>::value && details::is_complex<TY>::value == false, 
        class Main = details::true_t>
struct x_abs
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TY& a);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TY& a);
};

//-----------------------------------------------------------------------
//                      apx_abs_test_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a + abs(x)
*  test special cases (i.e. a = 0)
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y array and the scalar a
*  TX           - type of elements in X array
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
*  X            - pointer to X array
*  X_ld         - leading dimension of X matrix
*  rows         - number of rows of X and Y
*  cols         - number of columns of X and Y
*  a            - the scalar a
*/
template<bool Select, class TY, class TX, Integer Rows, Integer Cols, Integer Continuous, 
         class Main = details::true_t>
struct apx_abs_test_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, Integer rows, Integer cols,
                     const TY& a);
};


};};

#include "matcl-blas-lapack/level1/details/level1_other.inl"
