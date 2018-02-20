/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/matrix/matrix.h"

namespace matcl
{

//--------------------------------------------------------------
//      functions operating on rows or columns
//--------------------------------------------------------------

// number of nonzero elements along dimension dim
// return Integer matrix of size 1 x A.cols() if dim == 1
// or Integer matrix of size A.rows() x 1 if dim == 2
MATCL_MATREP_EXPORT Matrix  nnz(const Matrix& A, Integer dim);

// check if all elements in a matrix are nonzero; if dim == 1 then
// function operates on columns and return Integer matrix of 
// size 1 x A.cols(); if dim == 2 then function operators on rows
// and return Integer matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  all(const Matrix& v,int dim = 1);

// check if all elements in a matrix pass the test function t; 
// if dim == 1 then function operates on columns and return
// Integer matrix of size 1 x A.cols(); if dim == 2 then function
// operators on rows and return Integer matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  all(const Matrix& v,const test_function& t,int dim = 1);

// check if any element in a matrix is nonzero; if dim == 1 then
// function operates on columns and return Integer matrix of 
// size 1 x A.cols(); if dim == 2 then function operators on rows
// and return Integer matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  any(const Matrix& v,int dim = 1);

// check if any element in a matrix pass the test function t; 
// if dim == 1 then function operates on columns and return
// Integer matrix of size 1 x A.cols(); if dim == 2 then function
// operators on rows and return Integer matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  any(const Matrix& v,const test_function& t,int dim = 1);

// sum of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  sum(const Matrix& v, int dim = 1);

// product of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  prod(const Matrix& v, int dim = 1);

// cumulative sum of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  cumsum(const Matrix& v, int dim = 1);

// cumulative product of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  cumprod(const Matrix& v, int dim = 1);

// mean value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
MATCL_MATREP_EXPORT Matrix  mean(const Matrix& v, int dim = 1);

// standard deviation of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// if unbiaded == true, then unbiased estimator of standard deviation
// is returned (i.e. normalization by N-1 if N > 1, where N is the sample
// size), otherwise biased estimator is returned (normalization by N)
MATCL_MATREP_EXPORT Matrix  std(const Matrix& v, int dim = 1, bool unbiased = true);

// minimum value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// complex numbers are compared by real part and then by imaginary part
MATCL_MATREP_EXPORT Matrix  min_d(const Matrix& v, int dim = 1);

// maximum value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// complex numbers are compared by real part and then by imaginary part
MATCL_MATREP_EXPORT Matrix  max_d(const Matrix& v, int dim = 1);

// minimum absolute value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
MATCL_MATREP_EXPORT Matrix  min_abs_d(const Matrix& v, int dim = 1);

// maximum absolute value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
MATCL_MATREP_EXPORT Matrix  max_abs_d(const Matrix& v, int dim = 1);

// minimum value of elements in a matrix and positions of minimum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// complex numbers are compared by real part and then by imaginary part;
// return a 2-tupple of minimum values (via .get<1>()) and a 
// indices (via .get<2>())
MATCL_MATREP_EXPORT mat_tup_2 min2(const Matrix& v, int dim = 1);

// maximum value of elements in a matrix and positions of maximum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// complex numbers are compared by real part and then by imaginary part;
// return a 2-tupple of maximum values (via .get<1>()) and a 
// indices (via .get<2>())
MATCL_MATREP_EXPORT mat_tup_2 max2(const Matrix& v, int dim = 1);

// minimum absolute value of elements in a matrix and positions of minimum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// complex numbers are compared by real part and then by imaginary part;
// return a 2-tupple of minimum values (via .get<1>()) and a 
// indices (via .get<2>())
MATCL_MATREP_EXPORT mat_tup_2 min_abs2(const Matrix& v, int dim = 1);

// maximum absolute value of elements in a matrix and positions of maximum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// complex numbers are compared by real part and then by imaginary part;
// return a 2-tupple of maximum values (via .get<1>()) and a 
// indices (via .get<2>())
MATCL_MATREP_EXPORT mat_tup_2 max_abs2(const Matrix& v, int dim = 1);

// function equivalent to min_d
inline Matrix               min(const Matrix& v)    { return min_d(v,1); };

// function equivalent to max_d
inline Matrix               max(const Matrix& v)    { return max_d(v,1); };

// function equivalent to min_abs_d
inline Matrix               min_abs(const Matrix& v){ return min_abs_d(v,1); };

// function equivalent to max_abs_d
inline Matrix               max_abs(const Matrix& v){ return max_abs_d(v,1); };

//--------------------------------------------------------------
//      functions operating on all elements
//--------------------------------------------------------------

// number of nonzero elements;
MATCL_MATREP_EXPORT Integer nnz_vec(const Matrix& A);

// check if all elements in a matrix are nonzero; 
// return true for empty matrix
MATCL_MATREP_EXPORT bool    all_vec(const Matrix& v);

// check if all elements in a matrix pass the test function t; 
// return true for empty matrix
MATCL_MATREP_EXPORT bool    all_vec(const Matrix& v,const test_function& t);

// check if any element in a matrix is nonzero;
// return false for empty matrix
MATCL_MATREP_EXPORT bool    any_vec(const Matrix& v);

// check if any element in a matrix pass the test function t; 
// return false for empty matrix
MATCL_MATREP_EXPORT bool    any_vec(const Matrix& v,const test_function& t);

// sum of elements in a matrix; return scalar
// return zero for empty matrix
MATCL_MATREP_EXPORT Matrix  sum_vec(const Matrix& v);

// product of elements in a matrix; return scalar
// return zero for empty matrix
MATCL_MATREP_EXPORT Matrix  prod_vec(const Matrix& v);

// mean value of elements in a matrix; return scalar
// return zero for empty matrix
MATCL_MATREP_EXPORT Matrix  mean_vec(const Matrix& v);

// standard deviation of elements in a matrix; return scalar
// if unbiaded == true, then unbiased estimator of standard deviation
// is returned (i.e. normalization by N-1 if N > 1, where N is the sample
// size), otherwise biased estimator is returned (normalization by N)
// return zero for empty matrix
MATCL_MATREP_EXPORT Matrix  std_vec(const Matrix& v, bool unbiased = true);

// minimum value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
MATCL_MATREP_EXPORT Matrix  min_vec(const Matrix& v);

// maximum value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
MATCL_MATREP_EXPORT Matrix  max_vec(const Matrix& v);

// minimum absolute value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
MATCL_MATREP_EXPORT Matrix  min_abs_vec(const Matrix& v);

// maximum absolute value of elements in a matrix; 
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
MATCL_MATREP_EXPORT Matrix  max_abs_vec(const Matrix& v);

// tuple type
using mat_int2  = tuple<Matrix,Integer,Integer>;

// minimum value of elements in a matrix and positions of minimum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of minimum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
MATCL_MATREP_EXPORT mat_int2   min2_vec(const Matrix& v);

// maximum value of elements in a matrix and positions of maximum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of maximum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
MATCL_MATREP_EXPORT mat_int2   max2_vec(const Matrix& v);

// minimum absolute value of elements in a matrix and positions of minimum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of minimum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
MATCL_MATREP_EXPORT mat_int2   min_abs2_vec(const Matrix& v);

// maximum absolute value of elements in a matrix and positions of maximum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of maximum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
MATCL_MATREP_EXPORT mat_int2   max_abs2_vec(const Matrix& v);

}