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

#include "matcl-scalar/lib_functions/manip.h"
#include "matcl-matrep/matrix/matrix.h"

#include <vector>

namespace matcl
{

// delete rows specified by colon c.
MATCL_MATREP_EXPORT Matrix  delrows(const Matrix& A, const colon& c);
MATCL_MATREP_EXPORT Matrix  delrows(Matrix&& A, const colon& c);

// delete columns specified in colon c.
MATCL_MATREP_EXPORT Matrix  delcols(const Matrix& A, const colon& c);
MATCL_MATREP_EXPORT Matrix  delcols(const Matrix&& A, const colon& c);

// delete rows specified by colon c1 and columns specified by colon c2.
MATCL_MATREP_EXPORT Matrix  delrowscols(const Matrix& A, const colon& c1, const colon& c2);
MATCL_MATREP_EXPORT Matrix  delrowscols(Matrix&& A, const colon& c1, const colon& c2);

// concatenate matrices horizontally (left to right).
MATCL_MATREP_EXPORT Matrix  horzcat(const Matrix& A, const Matrix& B);
MATCL_MATREP_EXPORT Matrix  horzcat(std::initializer_list<Matrix> mat_list);
MATCL_MATREP_EXPORT Matrix  horzcat(const std::vector<Matrix>& mat_list);

// concatenate matrices horizontally (top to bottom).
MATCL_MATREP_EXPORT Matrix  vertcat(const Matrix& A, const Matrix& B);
MATCL_MATREP_EXPORT Matrix  vertcat(std::initializer_list<Matrix> mat_list);
MATCL_MATREP_EXPORT Matrix  vertcat(const std::vector<Matrix>& mat_list);

// create block diagonal matrix with given blocks on the main diagonal
MATCL_MATREP_EXPORT Matrix  blkdiag(const Matrix& A, const Matrix& B);
MATCL_MATREP_EXPORT Matrix  blkdiag(std::initializer_list<Matrix> mat_list);
MATCL_MATREP_EXPORT Matrix  blkdiag(const std::vector<Matrix>& mat_list);

// form a block matrix of size m by n, with a copy of matrix A as each element.
MATCL_MATREP_EXPORT Matrix  repmat(const Matrix &A, Integer m, Integer n);

// convers matrix m into dense matrix.
MATCL_MATREP_EXPORT Matrix  full(const Matrix &m);

// convert matrix m into sparse matrix.
MATCL_MATREP_EXPORT Matrix  sparse(const Matrix &m);

// convert matrix m into band matrix.
MATCL_MATREP_EXPORT Matrix  band(const Matrix &m);

// make independent copy
inline       Matrix         clone(const Matrix& m)      { return m.clone();};

// transpose given matrix m.
MATCL_MATREP_EXPORT Matrix  trans(const Matrix &m);

// makes conjugate transpose of given m.
MATCL_MATREP_EXPORT Matrix  ctrans(const Matrix &m);

// transpose or conjugate transpose of given matrix m.
MATCL_MATREP_EXPORT Matrix  trans(const Matrix &m, trans_type t);

// conjugate, transpose or conjugate transpose of given matrix m.
MATCL_MATREP_EXPORT Matrix  trans(const Matrix &m, trans_type_ext t);

// arranges matrix columns from left to right into vector placing 
// them from top to bottom.
MATCL_MATREP_EXPORT Matrix  vec(const Matrix &m);

// returns a new matrix formed by extracting the lower triangular part
// of the matrix m, and settings all other elements to zero;
// the second argument is optional, and specifies how many diagonals 
// above the diagonal d should also be set to zero (d = 0 is the main diagonal,
// d > 0 are superdiagonals, d < 0 are subdiagonals)
MATCL_MATREP_EXPORT Matrix  tril(const Matrix &m, Integer d = 0);
MATCL_MATREP_EXPORT Matrix  tril(Matrix &&m, Integer d = 0);

// returns a new matrix formed by extracting the upper triangular part of
// the matrix m, and settings all other elements to zero; the second argument
// is optional, and specifies how many diagonals below the diagonal d should 
// also be set to zero (d = 0 is the main diagonal, d > 0 are superdiagonals, 
// d < 0 are subdiagonals)
MATCL_MATREP_EXPORT Matrix  triu(const Matrix &m,Integer d = 0);
MATCL_MATREP_EXPORT Matrix  triu(Matrix &&m,Integer d = 0);

// returns a new matrix formed by extracting the diagonals with index in range
// [fd, ld] of the matrix m, and settings all other elements to zero
MATCL_MATREP_EXPORT Matrix  select_band(const Matrix &m, Integer fd, Integer ld);

// return a copy of m with the order of the rows reversed.
MATCL_MATREP_EXPORT Matrix  flipud(const Matrix &m);

// return a copy of m with the order of the columns reversed.
MATCL_MATREP_EXPORT Matrix  fliplr(const Matrix &m);

// return a matrix with the specified dimensions m x n.
MATCL_MATREP_EXPORT Matrix  reshape(const Matrix &A,Integer m,Integer n);

// return vector of elements on diagonal d; rhe second argument is optional, 
// and specifies which diagonal is to be returned (d = 0 is the main diagonal, 
// d > 0 are superdiagonals, d < 0 are subdiagonals)
MATCL_MATREP_EXPORT Matrix  get_diag(const Matrix& m, Integer d = 0);

// counterclockwise rotation of matrix by n x 90 degree
MATCL_MATREP_EXPORT Matrix  rot90(const Matrix& m, Integer n = 1);

// convert the matrix to a new type given by new_type; warning will be
// printed if this conversion leads to precision lost
MATCL_MATREP_EXPORT Matrix  convert(const Matrix& A, matcl::mat_code new_type);

// convert value type of stored elements not changing matrix structure; 
// warning will be printed if this conversion leads to precision lost
MATCL_MATREP_EXPORT Matrix  convert_value(const Matrix& A, matcl::value_code vc);

// convert a scalar of type From to a scalar of type To;
// this function is enabled if To and From are matcl scalars
// warnings are not printed if this conversion leads to precision lost
// redeclaration of function defined in matcl-scalar
template<class To, class From, class Enable>
To                          convert_scalar(const From& s);

// convert objects stored in the matrix to the type ti; if elements stored
// in the matrix are not of Object type, then original matrix is returned
MATCL_MATREP_EXPORT Matrix  convert_object(const Matrix& A, ti::ti_object ti);

// get number of nonzero subdiagonals in the matrix; if min != -1, then 
// searching for nonzero diagonal is stopped, if number of subdiagonals is
// higher than min
MATCL_MATREP_EXPORT Integer get_ld(const Matrix& A, Integer min = -1);

// get number of nonzero superdiagonals in the matrix; if min != -1, then 
// searching for nonzero diagonal is stopped, if number of superdiagonals is
// higher than min
MATCL_MATREP_EXPORT Integer get_ud(const Matrix& A, Integer min = -1);

// check if a matrix is lower triangular
MATCL_MATREP_EXPORT bool    is_tril(const Matrix& A);

// check if a matrix is upper triangular
MATCL_MATREP_EXPORT bool    is_triu(const Matrix& A);

// check if a matrix is diagonal
MATCL_MATREP_EXPORT bool    is_diag(const Matrix& A);

// check if a matrix is symmetric; tolerance tol is used to
// check equality of two elements, i.e. a and b are treated as
// equal if tol > 0 and |a-b| < eps(a) * tol or tol <= 0 and
// a == b; if a and b are integers, then a,b are treated as equal
// only if a == b
MATCL_MATREP_EXPORT bool    is_sym(const Matrix& A, Real tol);

// check if a matrix is hermitian; tolerance tol is used to
// check equality of two elements, i.e. a and b are treated as
// equal if tol > 0 and |a-b| < eps(a) * tol or tol <= 0 and
// a == b; if a and b are integers, then a,b are treated as equal
// only if a == b
MATCL_MATREP_EXPORT bool    is_her(const Matrix& A, Real tol);

// check if matrix stores real float or integer scalars
MATCL_MATREP_EXPORT bool    is_real_matrix(const Matrix& A);

// check if matrix stores real float scalars
MATCL_MATREP_EXPORT bool    is_real_float_matrix(const Matrix& A);

// check if matrix stores integer scalars
MATCL_MATREP_EXPORT bool    is_integer_matrix(const Matrix& A);

// check if matrix stores complex scalars
MATCL_MATREP_EXPORT bool    is_complex_matrix(const Matrix& A);

// check if matrix stores object scalars
MATCL_MATREP_EXPORT bool    is_object_matrix(const Matrix& A);

// check if matrix has dense representation
MATCL_MATREP_EXPORT bool    is_dense_matrix(const Matrix& A);

// check if matrix has sparse representation
MATCL_MATREP_EXPORT bool    is_sparse_matrix(const Matrix& A);

// check if matrix has band representation
MATCL_MATREP_EXPORT bool    is_band_matrix(const Matrix& A);

// check if matrix has scalar representation
MATCL_MATREP_EXPORT bool    is_scalar_matrix(const Matrix& A);

// check if two matrices A and B point to the same internal representation
MATCL_MATREP_EXPORT bool    is_same_matrix(const Matrix& A, const Matrix& B);

// get number of rows of a matrix
inline Integer              rows(const Matrix& A)               { return A.rows(); };

// get number of columns of a matrix
inline Integer              cols(const Matrix& A)               { return A.cols(); };

// get number of subdiagonals of a matrix based on representation
// and struct flags (if use_flags = true) only
inline Integer              structural_ldiags(const Matrix& A, bool use_flags = true)  
                                                                { return A.structural_ldiags(use_flags);};

// get number of superdiagonals of a matrix based on representation
// and struct flags (if use_flags = true) only
inline Integer              structural_udiags(const Matrix& A, bool use_flags = true)
                                                                { return A.structural_udiags(use_flags); };

// get number of nonzero elements in a matrix based on representation
// and struct flags only
inline Integer              structural_nnz(const Matrix& A)     { return A.structural_nnz(); };

// larger or rows() and cols(); return zero for empty matrices
inline Integer              length(const Matrix& A)             { return A.length(); };

// number of rows times number of columns
inline Real                 numel(const Matrix& A)              { return A.numel(); };

// check if all elements are finite
inline bool                 all_finite(const Matrix& A)         { return A.all_finite(); };

// number of nonzero elements
MATCL_MATREP_EXPORT Integer nnz(const Matrix& A);

// remove elements v in a sparse matrix A satisfying abs(v) <= tol.
// do nothing for dense and band matrices
MATCL_MATREP_EXPORT Matrix  drop_sparse(const Matrix& A, Real tol);

// search for the non-zero elements in v; return indices of nonzero
// elements
MATCL_MATREP_EXPORT Matrix  find(const Matrix& v);

// searches for the non-zero elements in v, using a user-defined test function;
// return indices of nonzero elements
MATCL_MATREP_EXPORT Matrix  find(const Matrix& v,const test_function& t);

// search for the non-zero elements in v; return a 2-tupple containing row and 
// column indices of nonzero elements
MATCL_MATREP_EXPORT mat_tup_2 find2(const Matrix& v);

// search for the non-zero elements in v, using a user-defined test function
// return a 2-tupple containing row and column indices of nonzero elements
MATCL_MATREP_EXPORT mat_tup_2 find2(const Matrix& v,const test_function& t);

// search for the non-zero elements in v; return a 3-tupple containing vectors
// of row (via .get<1>()), column (via .get<2>()) indices of non-zero elements 
// and the non-zero elements (via .get<3>())
MATCL_MATREP_EXPORT mat_tup_3 find3(const Matrix& v);

// search for the non-zero elements in v, using a user-defined test function
// return a 3-tupple containing vectors of row (via .get<1>()), column 
// (via .get<2>()) indices of non-zero elements  and the non-zero elements 
// (via .get<3>())
MATCL_MATREP_EXPORT mat_tup_3 find3(const Matrix& v,const test_function& t);

// sort the matrix v along dimension dim in asceding or desceding order;
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT Matrix  sort(const Matrix& v, int dim = 1, bool asceding = true);

// sort the matrix v along dimension dim in asceding or desceding order;
// return a 2-tupple with the sorted matrix (via .get<1>()) and a matrix 
// with the indices after sorting permutation (via .get<2>())
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT mat_tup_2 sort2(const Matrix& v, int dim = 1, bool asceding = true);

// sort rows of the matrix v as a group;
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT Matrix  sortrows(const Matrix& v);

// sort rows of the matrix v as a group; return a 2-tupple with the sorted 
// matrix (via .get<1>()) and a vector of indices of rows after sorting 
// permutation (via .get<2>()); complex numbers are sorted by real part 
// and then by imaginary part
MATCL_MATREP_EXPORT mat_tup_2 sortrows2(const Matrix& v);

// sort rows of the matrix v as a group considering only elements in columns
// cols; do ascending sort when column number is positive, descending if negative;
// for example if cols = [2, -1] then the function use second column for primary
// comparison (ascending) and first column for secondary comparison (descending);
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT Matrix  sortrows(const Matrix& v, const Matrix& cols);

// sort rows of the matrix v as a group considering only elements in columns
// cols; do ascending sort when column number is positive, descending if negative;
// for example if cols = [2, -1] then the function use second column for primary
// comparison (ascending) and first column for secondary comparison (descending);
// return a 2-tupple with the sorted matrix (via .get<1>()) and a vector of indices
// of rows after sorting permutation (via .get<2>());
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT mat_tup_2 sortrows2(const Matrix& v, const Matrix& cols);

// sort columns of the matrix v as a group;
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT Matrix  sortcols(const Matrix& v);

// sort columns of the matrix v as a group; return a 2-tupple with the sorted 
// matrix (via .get<1>()) and a vector of indices of colums after sorting 
// permutation (via .get<2>()); complex numbers are sorted by real part 
// and then by imaginary part
MATCL_MATREP_EXPORT mat_tup_2 sortcols2(const Matrix& v);

// sort columns of the matrix v as a group considering only elements in rows
// 'rows'; do ascending sort when row number is positive, descending if negative;
// for example if rows = [2, -1] then the function use second row for primary
// comparison (ascending) and first row for secondary comparison (descending);
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT Matrix  sortcols(const Matrix& v, const Matrix& dims);

// sort columns of the matrix v as a group considering only elements in rows
// 'rows'; do ascending sort when row number is positive, descending if negative;
// for example if rows = [2, -1] then the function use second row for primary
// comparison (ascending) and first row for secondary comparison (descending);
// return a 2-tupple with the sorted matrix (via .get<1>()) and a vector of indices
// of columns after sorting permutation (via .get<2>());
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT mat_tup_2 sortcols2(const Matrix& v, const Matrix& dims);

// check if a matrix v is sorted along dimension dim; for dim = 1 issorted(v) 
// returns a boolean vector, specifying if columns of v are sorted; for dim = 2
// function does its job for transposed matrix;
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT Matrix  issorted(const Matrix& v, int dim = 1, bool asceding = true);

// check if rows are sorted as a group;
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT bool    issorted_rows(const Matrix& v);

// check if columns are sorted as a group;
// complex numbers are sorted by real part and then by imaginary part
MATCL_MATREP_EXPORT bool    issorted_cols(const Matrix& v);

}
