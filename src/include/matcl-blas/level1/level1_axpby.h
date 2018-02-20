/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-blas/level1/config.h"
#include "matcl-blas/level1/utils.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace level1
{

//=======================================================================
// Form:
//
//          Y = a * X + b * Y
//
//  and all special cases
//=======================================================================

//-----------------------------------------------------------------------
//                          GENERAL ASSUMPTIONS
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
//                      axpby_test
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * X + b * Y
*  test all special cases for steps and a, b (test values 0, +1, -1)
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
*  TB           - type of b scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y             - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X               - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  a            - scalar
*  b            - scalar
*/
template<class TY, class TX, class TA, class TB, Integer Rows, 
        class Main = details::true_t>
struct axpby_test
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TA& a, const TB& b);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TA& a, const TB& b);
};

//-----------------------------------------------------------------------
//                      axpby_test_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * X + b * Y
*  test all special cases for a, b (test values 0, +1, -1)
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
*  TB           - type of b scalar
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
*  a            - scalar
*  b            - scalar
*/
template<bool Select, class TY, class TX, class TA, class TB, Integer Rows, Integer Cols, Integer Continuous, 
         class Main = details::true_t>
struct axpby_test_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, Integer rows, Integer cols, 
                     const TA& a, const TB& b);
};

//-----------------------------------------------------------------------
//                      axpy_test_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * X + Y
*  test all special cases for a (test values 0, +1, -1)
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
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
*  a            - scalar
*/
template<bool Select, class TY, class TX, class TA, Integer Rows, Integer Cols, Integer Continuous, 
         class Main = details::true_t>
struct axpy_test_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, Integer rows, Integer cols, 
                     const TA& a);
};

//-----------------------------------------------------------------------
//                      ax_test_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * X
*  test all special cases for a (test values 0, +1, -1)
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of the scalar a
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
*  X            - pointer to Y array
*  X_ld         - leading dimension of Y matrix
*  rows         - number of rows of X and Y
*  cols         - number of columns of X and Y
*  a            - scalar
*/
template<bool Select, class TY, class TX, class TA, Integer Rows, Integer Cols, Integer Continuous, 
         class Main = details::true_t>
struct ax_test_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, Integer rows, Integer cols,
                     const TA& a);
};

//-----------------------------------------------------------------------
//                      ay_test_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * Y
*  test all special cases for a (test values 0, +1, -1)
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y array
*  TA           - type of the scalar a
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
*/
template<bool Select, class TY, class TA, Integer Rows, Integer Cols, Integer Continuous, 
         class Main = details::true_t>
struct ay_test_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, Integer rows, Integer cols, const TA& a);
};

//-----------------------------------------------------------------------
//                      axpby
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * X + b * Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
*  TB           - type of b scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  a            - scalar, tests for special cases are not performed
*  b            - scalar, tests for special cases are not performed
*/
template<class TY, class TX, class TA, class TB, Integer Rows, 
        bool Use_simd = details::check_simd_scal2<TA,TB,TY,TX>::value,
        class Main = details::true_t>
struct axpby
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TA& a, const TB& b);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TA& a, const TB& b);
};

//-----------------------------------------------------------------------
//                          ax
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * X
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  a            - the scalar, tests for special cases are not performed
*/
template<class TY, class TX, class TA, Integer Rows, 
        bool Use_simd = details::check_simd_scal<TA,TY,TX>::value, class Main = details::true_t>
struct ax
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TA& a);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, Integer rows, const TA& a);
};

//-----------------------------------------------------------------------
//                          mx
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = - X
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
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*/
template<class TY, class TX, Integer Rows, 
        bool Use_simd = details::check_simd<TY,TX>::value, class Main = details::true_t>
struct mx
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, Integer rows);
};

//-----------------------------------------------------------------------
//                          ay
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array and the scalar a
*  TA           - type of a scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  rows         - number of array
*  a            - the scalar, tests for special cases are not performed
*/
template<class TY, class TA, Integer Rows, 
        bool Use_simd = details::check_simd_scal<TA, TY>::value, class Main = details::true_t>
struct ay
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer rows, const TA& a);
    static void eval(TY* Y, Integer Y_step, Integer rows, const TA& a);
};

//-----------------------------------------------------------------------
//                          my
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = - Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  rows         - number of array
*/
template<class TY, Integer Rows, 
        bool Use_simd = details::check_simd<TY>::value, class Main = details::true_t>
struct my
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer rows);
    static void eval(TY* Y, Integer Y_step, Integer rows);
};

//-----------------------------------------------------------------------
//                          xmy
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = X - Y
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
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*/
template<class TY, class TX, Integer Rows, 
        bool Use_simd = details::check_simd<TY,TX>::value, class Main = details::true_t>
struct xmy
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, Integer rows);
};

//-----------------------------------------------------------------------
//                          ymx
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = Y - X
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
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*/
template<class TY, class TX, Integer Rows, 
        bool Use_simd = details::check_simd<TY,TX>::value, class Main = details::true_t>
struct ymx
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, Integer rows);
};

//-----------------------------------------------------------------------
//                      ymx_mat
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = Y - X
*
*  Template arguments
*  =======
*  Select       - choose version dynamically based on values of arguments,
*                 true or false.
*  TY           - type of elements in Y array
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
*/
template<bool Select, class TY, class TX, Integer Rows, Integer Cols, Integer Continuous, 
         class Main = details::true_t>
struct ymx_mat
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, Integer Y_ld, const TX* X, Integer X_ld, Integer rows, Integer cols);
};

//-----------------------------------------------------------------------
//                          ypx
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = Y + X
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
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*/
template<class TY, class TX, Integer Rows, 
        bool Use_simd = details::check_simd<TY,TX>::value, class Main = details::true_t>
struct ypx
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, Integer rows);
};

//-----------------------------------------------------------------------
//                          ypxm
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = -(Y + X)
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
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*/
template<class TY, class TX, Integer Rows, 
        bool Use_simd = details::check_simd<TY,TX>::value, class Main = details::true_t>
struct ypxm
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, Integer rows);
};

//-----------------------------------------------------------------------
//                          axpy
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * X + Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  a            - scalar, tests for special cases are not performed
*/
template<class TY, class TX, class TA, Integer Rows, 
            bool Use_simd = details::check_simd_scal<TA,TY,TX>::value, class Main = details::true_t>
struct axpy
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TA& a);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TA& a);
};

//-----------------------------------------------------------------------
//                          axmy
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * X - Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  a            - scalar, tests for special cases are not performed
*/
template<class TY, class TX, class TA, Integer Rows, 
            bool Use_simd = details::check_simd_scal<TA,TY,TX>::value, class Main = details::true_t>
struct axmy
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TA& a);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TA& a);
};

//-----------------------------------------------------------------------
//                          xpya
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * (X + Y)
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  a            - scalar, tests for special cases are not performed
*/
template<class TY, class TX, class TA, Integer Rows, 
            bool Use_simd = details::check_simd_scal<TA,TY,TX>::value, class Main = details::true_t>
struct xpya
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TA& a);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TA& a);
};

//-----------------------------------------------------------------------
//                          xmya
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = a * (X - Y)
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TA           - type of a scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  a            - scalar, tests for special cases are not performed
*/
template<class TY, class TX, class TA, Integer Rows, 
            bool Use_simd = details::check_simd_scal<TA,TY,TX>::value, class Main = details::true_t>
struct xmya
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TA& a);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TA& a);
};

//-----------------------------------------------------------------------
//                          xpby
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = X + b * Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TB           - type of b scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  b            - scalar, tests for special cases are not performed
*/
template<class TY, class TX, class TB, Integer Rows, 
            bool Use_simd = details::check_simd_scal<TB,TY,TX>::value, class Main = details::true_t>
struct xpby
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TB& b);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TB& b);
};

//-----------------------------------------------------------------------
//                          mxpby
//-----------------------------------------------------------------------
/**
*  Purpose
*  =======
*  form Y = -X + b * Y
*
*  Template arguments
*  =======
*  TY           - type of elements in Y array
*  TX           - type of elements in X array
*  TB           - type of b scalar
*  Rows         - number of elements, 0 if unknown statically
*
*  Inputs
*  =======
*  Y            - pointer to Y array
*  Y_step       - optional step between two elements in Y array
*  X            - pointer to X array
*  X_step       - optional step between two elements in X array
*  rows         - number of elements
*  b            - scalar, tests for special cases are not performed
*/
template<class TY, class TX, class TB, Integer Rows, 
            bool Use_simd = details::check_simd_scal<TB,TY,TX>::value, class Main = details::true_t>
struct mxpby
{
    static_assert(md::dependent_false<TY>::value,
                    "this version not implemented, probably one of argument is invalid");

    static void eval(TY* Y, const TX* X, Integer rows, const TB& b);
    static void eval(TY* Y, Integer Y_step, const TX* X, Integer X_step, 
                     Integer rows, const TB& b);
};

}}

#include "matcl-blas/level1/details/level1_axpby.inl"
