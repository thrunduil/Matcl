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

#include "matcl-matrep/objects/object_matrix.h"

namespace matcl
{

// -----------------------------------------------------------------------
//                 IO functions
// -----------------------------------------------------------------------

// display matrix using global disp_stream as default, which prints on global 
// output stream. Options controls how printing is performed, see options_disp
// for details
template<class T>
void                        disp(const object_matrix<T>& m, const disp_stream_ptr& os = default_disp_stream(),
                                const options& opts = options());

// output stream operator for a matrix. This operator is intendent to
// save matrices in text format, not for pretty printing. However output looks
// reasonably good. If you want to get nice looking output, use disp function
// with appropriate output stream.
template<class T>
std::ostream&               operator<<(std::ostream&, const object_matrix<T>&);

// input stream operator for a matrix. This operator is intendent to
// load matrices in text format, saved previously with operator<<..
template<class T>
std::istream&               operator>>(std::istream&, object_matrix<T>&);

// serialize matrices using boost::serialization library
template<class T>
void                        save(oarchive& ar,const object_matrix<T>& mat);

// deserialize matrices
template<class T>
void                        load(iarchive& ar,object_matrix<T>& mat);

// -----------------------------------------------------------------------
//                 manip functions
// -----------------------------------------------------------------------

// delete rows specified by colon c.
template<class T>
object_matrix<T>            delrows(const object_matrix<T>& A, const colon& c);
template<class T>
object_matrix<T>            delrows(object_matrix<T>&& A, const colon& c);

// delete columns specified in colon c.
template<class T>
object_matrix<T>            delcols(const object_matrix<T>& A, const colon& c);
template<class T>
object_matrix<T>            delcols(object_matrix<T>&& A, const colon& c);

// delete rows specified by colon c1 and columns specified by colon c2.
template<class T>
object_matrix<T>            delrowscols(const object_matrix<T>& A, const colon& c1, const colon& c2);
template<class T>
object_matrix<T>            delrowscols(object_matrix<T>&& A, const colon& c1, const colon& c2);

// concatenate matrices horizontally (left to right).
template<class T>
object_matrix<T>            horzcat(const object_matrix<T>& A, const object_matrix<T>& B);
template<class T>
object_matrix<T>            horzcat(std::initializer_list<object_matrix<T>> mat_list);
template<class T>
object_matrix<T>            horzcat(const std::vector<object_matrix<T>>& mat_list);

// concatenate matrices horizontally (top to bottom).
template<class T>
object_matrix<T>            vertcat(const object_matrix<T>& A, const object_matrix<T>& B);
template<class T>
object_matrix<T>            vertcat(std::initializer_list<object_matrix<T>> mat_list);
template<class T>
object_matrix<T>            vertcat(const std::vector<object_matrix<T>>& mat_list);

// create block diagonal matrix with given blocks on the main diagonal
template<class T>
object_matrix<T>            blkdiag(const object_matrix<T>& A, const object_matrix<T>& B);
template<class T>
object_matrix<T>            blkdiag(std::initializer_list<object_matrix<T>> mat_list);
template<class T>
object_matrix<T>            blkdiag(const std::vector<object_matrix<T>>& mat_list);

// form a block matrix of size m by n, with a copy of matrix A as each element.
template<class T>
object_matrix<T>            repmat(const object_matrix<T>& A, Integer m, Integer n);

// convers matrix m into dense matrix.
template<class T>
object_matrix<T>            full(const object_matrix<T>& m);

// convers matrix m into sparse matrix.
template<class T>
object_matrix<T>            sparse(const object_matrix<T>& m);

// convers matrix m into band matrix.
template<class T>
object_matrix<T>            band(const object_matrix<T>& m);

// make independent copy
template<class T>
object_matrix<T>            clone(const object_matrix<T>& m);

// transpose given matrix m.
template<class T>
object_matrix<T>            trans(const object_matrix<T>& m);

// makes conjugate transpose of given m.
template<class T>
object_matrix<T>            ctrans(const object_matrix<T>& m);

// transpose or conjugate transpose of given matrix m.
template<class T>
object_matrix<T>            trans(const object_matrix<T>& m, trans_type t);

// conjugate, transpose or conjugate transpose of given matrix m.
template<class T>
object_matrix<T>            trans(const object_matrix<T>& m, trans_type_ext t);

// arranges matrix columns from left to right into vector placing 
// them from top to bottom.
template<class T>
object_matrix<T>            vec(const object_matrix<T>& m);

// returns a new matrix formed by extracting the lower triangular part
// of the matrix m, and settings all other elements to zero;
// the second argument is optional, and specifies how many diagonals 
// above the diagonal d should also be set to zero (d = 0 is the main diagonal,
// d > 0 are superdiagonals, d < 0 are subdiagonals)
template<class T>
object_matrix<T>            tril(const object_matrix<T>& m, Integer d = 0);
template<class T>
object_matrix<T>            tril(object_matrix<T>&& m, Integer d = 0);

// returns a new matrix formed by extracting the upper triangular part of
// the matrix m, and settings all other elements to zero; the second argument
// is optional, and specifies how many diagonals below the diagonal d should 
// also be set to zero (d = 0 is the main diagonal, d > 0 are superdiagonals, 
// d < 0 are subdiagonals)
template<class T>
object_matrix<T>            triu(const object_matrix<T>& m, Integer d = 0);
template<class T>
object_matrix<T>            triu(object_matrix<T>&& m, Integer d = 0);

// returns a new matrix formed by extracting the upper ud diagonals and lower
// ld diagonals of the matrix m, and settings all other elements to zero
template<class T>
object_matrix<T>            select_band(const object_matrix<T> &m, Integer ld, Integer ud);

// return a copy of m with the order of the rows reversed.
template<class T>
object_matrix<T>            flipud(const object_matrix<T>& m);
template<class T>
object_matrix<T>            flipud(object_matrix<T>&& m);

// return a copy of m with the order of the columns reversed.
template<class T>
object_matrix<T>            fliplr(const object_matrix<T>& m);
template<class T>
object_matrix<T>            fliplr(object_matrix<T>&& m);

// return a matrix with the specified dimensions m x n.
template<class T>
object_matrix<T>            reshape(const object_matrix<T>& A, Integer m, Integer n);

// return vector of elements on diagonal d; rhe second argument is optional, 
// and specifies which diagonal is to be returned (d = 0 is the main diagonal, 
// d > 0 are superdiagonals, d < 0 are subdiagonals)
template<class T>
object_matrix<T>            get_diag(const object_matrix<T>& m, Integer d = 0);

// counterclockwise rotation of matrix by n x 90 degree
template<class T>
object_matrix<T>            rot90(const object_matrix<T>& m, Integer n = 1);

// search for the non-zero elements in v; return indices of nonzero
// elements
template<class T>
dense_matrix<Integer>       find(const object_matrix<T>& v);

// searches for the non-zero elements in v, using a user-defined test function;
// return indices of nonzero elements
template<class T>
dense_matrix<Integer>       find(const object_matrix<T>& v, const test_function& t);
// version with test function, that require implementing test for type T only
template<class T>
dense_matrix<Integer>       find(const object_matrix<T>& v, const test_object_function<T>& t);

// search for the non-zero elements in v; return a 2-tupple containing row and 
//column indices of nonzero elements
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const object_matrix<T>& v);

// search for the non-zero elements in v, using a user-defined test function
// return a 2-tupple containing row and column indices of nonzero elements
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const object_matrix<T>& v,const test_function& t);
// version with test function, that require implementing test for type T only
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const object_matrix<T>& v,const test_object_function<T>& t);

// search for the non-zero elements in v; return a 3-tupple containing vectors
// of row (via .get<1>()), column (via .get<2>()) indices of non-zero elements 
// and the non-zero elements (via .get<3>())
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>
                            find3(const object_matrix<T>& v);

// search for the non-zero elements in v, using a user-defined test function
// return a 3-tupple containing vectors of row (via .get<1>()), column 
// (via .get<2>()) indices of non-zero elements  and the non-zero elements 
// (via .get<3>())
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>
                            find3(const object_matrix<T>& v, const test_function& t);
// version with test function, that require implementing test for type T only
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, object_matrix<T>>
                            find3(const object_matrix<T>& v, const test_object_function<T>& t);

// sort the matrix v along dimension dim in asceding or desceding order;
template<class T>
object_matrix<T>            sort(const object_matrix<T>& v, int dim = 1, bool asceding = true);

// sort the matrix v along dimension dim in asceding or desceding order;
// return a 2-tupple with the sorted matrix (via .get<1>()) and a matrix 
// with the indices after sorting permutation (via .get<2>())
template<class T>
tuple<object_matrix<T>, dense_matrix<Integer>>
                            sort2(const object_matrix<T>& v, int dim = 1, bool asceding = true);

// sort rows of the matrix v as a group;
template<class T>
object_matrix<T>            sortrows(const object_matrix<T>& v);

// sort rows of the matrix v as a group; return a 2-tupple with the sorted 
// matrix (via .get<1>()) and a vector of indices of rows after sorting 
// permutation (via .get<2>());
template<class T>
tuple<object_matrix<T>, dense_matrix<Integer>>
                            sortrows2(const object_matrix<T>& v);

// sort rows of the matrix v as a group considering only elements in columns
// cols; do ascending sort when column number is positive, descending if negative;
// for example if cols = [2, -1] then the function use second column for primary
// comparison (ascending) and first column for secondary comparison (descending);
template<class T>
object_matrix<T>            sortrows(const object_matrix<T>& v, const Matrix& cols);

// sort rows of the matrix v as a group considering only elements in columns
// cols; do ascending sort when column number is positive, descending if negative;
// for example if cols = [2, -1] then the function use second column for primary
// comparison (ascending) and first column for secondary comparison (descending);
// return a 2-tupple with the sorted matrix (via .get<1>()) and a vector of indices
// of rows after sorting permutation (via .get<2>());
template<class T>
tuple<object_matrix<T>, dense_matrix<Integer>>
                            sortrows2(const object_matrix<T>& v, const Matrix& cols);

// sort columns of the matrix v as a group;
template<class T>
object_matrix<T>            sortcols(const object_matrix<T>& v);

// sort columns of the matrix v as a group; return a 2-tupple with the sorted 
// matrix (via .get<1>()) and a vector of indices of colums after sorting 
// permutation (via .get<2>());
template<class T>
tuple<object_matrix<T>, dense_matrix<Integer>>
                            sortcols2(const object_matrix<T>& v);

// sort columns of the matrix v as a group considering only elements in rows
// 'rows'; do ascending sort when row number is positive, descending if negative;
// for example if rows = [2, -1] then the function use second row for primary
// comparison (ascending) and first row for secondary comparison (descending);
template<class T>
object_matrix<T>            sortcols(const object_matrix<T>& v, const Matrix& dims);

// sort columns of the matrix v as a group considering only elements in rows
// 'rows'; do ascending sort when row number is positive, descending if negative;
// for example if rows = [2, -1] then the function use second row for primary
// comparison (ascending) and first row for secondary comparison (descending);
// return a 2-tupple with the sorted matrix (via .get<1>()) and a vector of indices
// of columns after sorting permutation (via .get<2>());
template<class T>
tuple<object_matrix<T>, dense_matrix<Integer>>
                            sortcols2(const object_matrix<T>& v, const Matrix& dims);

// check if a matrix v is sorted along dimension dim; for dim = 1 issorted(v) 
// returns a boolean vector, specifying if columns of v are sorted; for dim = 2
// function does its job for transposed matrix;
template<class T>
dense_matrix<Integer>       issorted(const object_matrix<T>& v, int dim = 1, bool asceding = true);

// remove elements v in a sparse matrix A satisfying abs(v) <= tol.
// do nothing for dense and band matrices
template<class T>
object_matrix<T>            drop_sparse(const object_matrix<T>& A, Real tol);

// convert objects stored in the matrix to object_type<S>; 
template<class S, class T>
object_matrix<S>            convert_object(const object_matrix<T>& A);

//--------------------------------------------------------------
//      functions operating on rows or columns
//--------------------------------------------------------------

// number of nonzero elements along dimension dim
// return Integer matrix of size 1 x A.cols() if dim == 1
// or Integer matrix of size A.rows() x 1 if dim == 2
template<class T>
dense_matrix<Integer>       nnz(const object_matrix<T>& A, Integer dim);

// check if all elements in a matrix are nonzero; if dim == 1 then
// function operates on columns and return Integer matrix of 
// size 1 x A.cols(); if dim == 2 then function operators on rows
// and return Integer matrix of size A.rows() x 1
template<class T>
dense_matrix<Integer>       all(const object_matrix<T>& v,int dim = 1);

// check if all elements in a matrix pass the test function t; 
// if dim == 1 then function operates on columns and return
// Integer matrix of size 1 x A.cols(); if dim == 2 then function
// operators on rows and return Integer matrix of size A.rows() x 1
template<class T>
dense_matrix<Integer>       all(const object_matrix<T>& v,const test_function& t, int dim = 1);

// version with test function, that require implementing test for type T only
template<class T>
dense_matrix<Integer>       all(const object_matrix<T>& v, const test_object_function<T>& t, int dim = 1);

// check if any element in a matrix is nonzero; if dim == 1 then
// function operates on columns and return Integer matrix of 
// size 1 x A.cols(); if dim == 2 then function operators on rows
// and return Integer matrix of size A.rows() x 1
template<class T>
dense_matrix<Integer>       any(const object_matrix<T>& v,int dim = 1);

// check if any element in a matrix pass the test function t; 
// if dim == 1 then function operates on columns and return
// Integer matrix of size 1 x A.cols(); if dim == 2 then function
// operators on rows and return Integer matrix of size A.rows() x 1
template<class T>
dense_matrix<Integer>       any(const object_matrix<T>& v,const test_function& t, int dim = 1);

// version with test function, that require implementing test for type T only
template<class T>
dense_matrix<Integer>       any(const object_matrix<T>& v, const test_object_function<T>& t, int dim = 1);

// sum of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
object_matrix<T>            sum(const object_matrix<T>& v, int dim = 1);

// product of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
object_matrix<T>            prod(const object_matrix<T>& v, int dim = 1);

// cumulative sum of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
object_matrix<T>            cumsum(const object_matrix<T>& v, int dim = 1);

// cumulative product of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
object_matrix<T>            cumprod(const object_matrix<T>& v, int dim = 1);

// mean value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
// this function is available only if S operator/(T,Integer) is defined
// additionally function object_type<S> div(object_type<T>,OInteger)
// must be registered
template<class T>
object_matrix<decltype(T()/Integer())>
                            mean(const object_matrix<T>& v, int dim = 1);

// standard deviation of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// if unbiaded == true, then unbiased estimator of standard deviation
// is returned (i.e. normalization by N-1 if N > 1, where N is the sample
// size), otherwise biased estimator is returned (normalization by N)
// this function is available only if S abs(T) and R operator/(S,Integer)
// are defined; additionally functions object_type<S> abs(object_type<T>)
// abs object_type<R> div(object_type<T>,OInteger) must be registered
template<class T>
object_matrix<decltype(decltype(abs(T()))()/Integer())>
                            std(const object_matrix<T>& v, int dim = 1, bool unbiased = true);

// minimum value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
template<class T>
object_matrix<T>            min_d(const object_matrix<T>& v, int dim = 1);

// maximum value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
template<class T>
object_matrix<T>            max_d(const object_matrix<T>& v, int dim = 1);

// minimum absolute value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
object_matrix<decltype(abs(T()))>
                            min_abs_d(const object_matrix<T>& v, int dim = 1);

// maximum absolute value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
object_matrix<decltype(abs(T()))>
                            max_abs_d(const object_matrix<T>& v, int dim = 1);

// minimum value of elements in a matrix and positions of minimum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// return a 2-tupple of minimum values (via .get<1>()) and a 
// indices (via .get<2>())
template<class T>
tuple<object_matrix<T>, dense_matrix<Integer>>
                            min2(const object_matrix<T>& v, int dim = 1);

// maximum value of elements in a matrix and positions of maximum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// return a 2-tupple of maximum values (via .get<1>()) and a 
// indices (via .get<2>())
template<class T>
tuple<object_matrix<T>, dense_matrix<Integer>>
                            max2(const object_matrix<T>& v, int dim = 1);

// minimum absolute value of elements in a matrix and positions of minimum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// return a 2-tupple of minimum values (via .get<1>()) and a 
// indices (via .get<2>())
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
tuple<object_matrix<decltype(abs(T()))>, dense_matrix<Integer>>
                            min_abs2(const object_matrix<T>& v, int dim = 1);

// maximum absolute value of elements in a matrix and positions of maximum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// return a 2-tupple of maximum values (via .get<1>()) and a 
// indices (via .get<2>())
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
tuple<object_matrix<decltype(abs(T()))>, dense_matrix<Integer>>
                            max_abs2(const object_matrix<T>& v, int dim = 1);

// function equivalent to min_d
template<class T>
inline object_matrix<T>     min(const object_matrix<T>& v)      { return min_d(v,1); };

// function equivalent to max_d
template<class T>
inline object_matrix<T>     max(const object_matrix<T>& v)      { return max_d(v,1); };

// function equivalent to min_abs_d
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
inline object_matrix<decltype(abs(T()))>
                            min_abs(const object_matrix<T>& v)  { return min_abs_d(v,1); };

// function equivalent to max_abs_d
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
inline object_matrix<decltype(abs(T()))>
                            max_abs(const object_matrix<T>& v)  { return max_abs_d(v,1); };

//--------------------------------------------------------------
//      functions operating on all elements
//--------------------------------------------------------------

// number of nonzero elements;
template<class T>
Integer                     nnz_vec(const object_matrix<T>& A);

// check if all elements in a matrix are nonzero; 
// return true for empty matrix
template<class T>
bool                        all_vec(const object_matrix<T>& v);

// check if all elements in a matrix pass the test function t; 
// return true for empty matrix
template<class T>
bool                        all_vec(const object_matrix<T>& v,const test_function& t);

// version with test function, that require implementing test for type T only
template<class T>
bool                        all_vec(const object_matrix<T>& v, const test_object_function<T>& t);

// check if any element in a matrix is nonzero;
// return false for empty matrix
template<class T>
bool                        any_vec(const object_matrix<T>& v);

// check if any element in a matrix pass the test function t; 
// return false for empty matrix
template<class T>
bool                        any_vec(const object_matrix<T>& v,const test_function& t);

// version with test function, that require implementing test for type T only
template<class T>
bool                        any_vec(const object_matrix<T>& v, const test_object_function<T>& t);

// sum of elements in a matrix; return scalar
// return zero for empty matrix
template<class T>
object_type<T>              sum_vec(const object_matrix<T>& v);

// product of elements in a matrix; return scalar
// return zero for empty matrix
template<class T>
object_type<T>              prod_vec(const object_matrix<T>& v);

// mean value of elements in a matrix; return scalar
// return zero for empty matrix
// this function is available only if S operator/(T,Integer) is defined
// additionally function object_type<S> div(object_type<T>,OInteger)
// must be registered
template<class T>
object_type<decltype(T()/Integer())>
                            mean_vec(const object_matrix<T>& v);

// standard deviation of elements in a matrix; return scalar
// if unbiaded == true, then unbiased estimator of standard deviation
// is returned (i.e. normalization by N-1 if N > 1, where N is the sample
// size), otherwise biased estimator is returned (normalization by N)
// this function is available only if S operator/(T,Integer) is defined
// this function is available only if S abs(T) and R operator/(S,Integer)
// are defined; additionally functions object_type<S> abs(object_type<T>)
// abs object_type<R> div(object_type<T>,OInteger) must be registered
template<class T>
object_type<decltype(decltype(abs(T()))()/Integer())>
                            std_vec(const object_matrix<T>& v, bool unbiased = true);

// minimum value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
template<class T>
object_type<T>              min_vec(const object_matrix<T>& v);

// maximum value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
template<class T>
object_type<T>              max_vec(const object_matrix<T>& v);

// minimum absolute value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
object_type<decltype(abs(T()))>
                            min_abs_vec(const object_matrix<T>& v);

// maximum absolute value of elements in a matrix; 
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
object_type<decltype(abs(T()))>
                            max_abs_vec(const object_matrix<T>& v);

// minimum value of elements in a matrix and positions of minimum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of minimum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
template<class T>
tuple<object_type<T>, Integer, Integer >
                            min2_vec(const object_matrix<T>& v);

// maximum value of elements in a matrix and positions of maximum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of maximum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
template<class T>
tuple<object_type<T>, Integer, Integer >
                            max2_vec(const object_matrix<T>& v);

// minimum absolute value of elements in a matrix and positions of minimum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of minimum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
tuple<object_type<decltype(abs(T()))>, Integer, Integer >
                            min_abs2_vec(const object_matrix<T>& v);

// maximum absolute value of elements in a matrix and positions of maximum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of maximum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
// this function is available only if S abs(T) is defined
// additionally function object_type<S> abs(object_type<T>) must be registered
template<class T>
tuple<object_type<decltype(abs(T()))>, Integer, Integer >
                            max_abs2_vec(const object_matrix<T>& v);

};
