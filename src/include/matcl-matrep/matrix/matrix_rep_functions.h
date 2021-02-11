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

#include "matcl-matrep/matrix/matrix_rep_dense.h"

namespace matcl
{

// -----------------------------------------------------------------------
//                 IO functions
// -----------------------------------------------------------------------

// display matrix using global disp_stream as default, which prints on global 
// output stream. Options controls how printing is performed, see options_disp
// for details
template<class T>
void                        disp(const dense_matrix<T>& m, const disp_stream_ptr& os = default_disp_stream(),
                                const options& opts = options());
template<class T>
void                        disp(const sparse_matrix<T>& m, const disp_stream_ptr& os = default_disp_stream(),
                                const options& opts = options());
template<class T>
void                        disp(const band_matrix<T>& m, const disp_stream_ptr& os = default_disp_stream(),
                                const options& opts = options());

// output stream operator for a matrix. This operator is intendent to
// save matrices in text format, not for pretty printing. However output looks
// reasonably good. If you want to get nice looking output, use disp function
// with appropriate output stream.
template<class T>
std::ostream&               operator<<(std::ostream&, const dense_matrix<T>&);
template<class T>
std::ostream&               operator<<(std::ostream&, const sparse_matrix<T>&);
template<class T>
std::ostream&               operator<<(std::ostream&, const band_matrix<T>&);

// input stream operator for a matrix. This operator is intendent to
// load matrices in text format, saved previously with operator<<..
template<class T>
std::istream&               operator>>(std::istream&, dense_matrix<T>&);
template<class T>
std::istream&               operator>>(std::istream&, sparse_matrix<T>&);
template<class T>
std::istream&               operator>>(std::istream&, band_matrix<T>&);

// serialize matrices using boost::serialization library
template<class T>
void                        save(oarchive& ar,const dense_matrix<T>& mat);
template<class T>
void                        save(oarchive& ar,const sparse_matrix<T>& mat);
template<class T>
void                        save(oarchive& ar,const band_matrix<T>& mat);

// deserialize matrices
template<class T>
void                        load(iarchive& ar,dense_matrix<T>& mat);
template<class T>
void                        load(iarchive& ar,sparse_matrix<T>& mat);
template<class T>
void                        load(iarchive& ar,band_matrix<T>& mat);

// -----------------------------------------------------------------------
//                 manip functions
// -----------------------------------------------------------------------

// delete rows specified by colon c.
template<class T>
dense_matrix<T>             delrows(const dense_matrix<T>& A, const colon& c);
template<class T>
sparse_matrix<T>            delrows(const sparse_matrix<T>& A, const colon& c);
template<class T>
sparse_matrix<T>            delrows(const band_matrix<T>& A, const colon& c);

template<class T>
dense_matrix<T>             delrows(dense_matrix<T>&& A, const colon& c);
template<class T>
sparse_matrix<T>            delrows(sparse_matrix<T>&& A, const colon& c);
template<class T>
sparse_matrix<T>            delrows(band_matrix<T>&& A, const colon& c);

// delete columns specified in colon c.
template<class T>
dense_matrix<T>             delcols(const dense_matrix<T>& A, const colon& c);
template<class T>
sparse_matrix<T>            delcols(const sparse_matrix<T>& A, const colon& c);
template<class T>
sparse_matrix<T>            delcols(const band_matrix<T>& A, const colon& c);

template<class T>
dense_matrix<T>             delcols(dense_matrix<T>&& A, const colon& c);
template<class T>
sparse_matrix<T>            delcols(sparse_matrix<T>&& A, const colon& c);
template<class T>
sparse_matrix<T>            delcols(band_matrix<T>&& A, const colon& c);

// delete rows specified by colon c1 and columns specified by colon c2.
template<class T>
dense_matrix<T>             delrowscols(const dense_matrix<T>& A, const colon& c1, const colon& c2);
template<class T>
sparse_matrix<T>            delrowscols(const sparse_matrix<T>& A, const colon& c1, const colon& c2);
template<class T>
sparse_matrix<T>            delrowscols(const band_matrix<T>& A, const colon& c1, const colon& c2);

template<class T>
dense_matrix<T>             delrowscols(dense_matrix<T>&& A, const colon& c1, const colon& c2);
template<class T>
sparse_matrix<T>            delrowscols(sparse_matrix<T>&& A, const colon& c1, const colon& c2);
template<class T>
sparse_matrix<T>            delrowscols(band_matrix<T>&& A, const colon& c1, const colon& c2);

// concatenate matrices horizontally (left to right).
template<class T>
dense_matrix<T>             horzcat(const dense_matrix<T>& A, const dense_matrix<T>& B);
template<class T>
dense_matrix<T>             horzcat(std::initializer_list<dense_matrix<T>> mat_list);
template<class T>
dense_matrix<T>             horzcat(const std::vector<dense_matrix<T>>& mat_list);
template<class T>
sparse_matrix<T>            horzcat(const sparse_matrix<T>& A, const sparse_matrix<T>& B);
template<class T>
sparse_matrix<T>            horzcat(std::initializer_list<sparse_matrix<T>> mat_list);
template<class T>
sparse_matrix<T>            horzcat(const std::vector<sparse_matrix<T>>& mat_list);
template<class T>
sparse_matrix<T>            horzcat(const sparse_matrix<T>& A, const band_matrix<T>& B);
template<class T>
sparse_matrix<T>            horzcat(const band_matrix<T>& A, const sparse_matrix<T>& B);
template<class T>
sparse_matrix<T>            horzcat(const band_matrix<T>& A, const band_matrix<T>& B);
template<class T>
sparse_matrix<T>            horzcat(std::initializer_list<band_matrix<T>> mat_list);
template<class T>
sparse_matrix<T>            horzcat(const std::vector<band_matrix<T>>& mat_list);

// concatenate matrices horizontally (top to bottom).
template<class T>
dense_matrix<T>             vertcat(const dense_matrix<T>& A, const dense_matrix<T>& B);
template<class T>
dense_matrix<T>             vertcat(std::initializer_list<dense_matrix<T>> mat_list);
template<class T>
dense_matrix<T>             vertcat(const std::vector<dense_matrix<T>>& mat_list);
template<class T>
sparse_matrix<T>            vertcat(const sparse_matrix<T>& A, const sparse_matrix<T>& B);
template<class T>
sparse_matrix<T>            vertcat(std::initializer_list<sparse_matrix<T>> mat_list);
template<class T>
sparse_matrix<T>            vertcat(const std::vector<sparse_matrix<T>>& mat_list);
template<class T>
sparse_matrix<T>            vertcat(const sparse_matrix<T>& A, const band_matrix<T>& B);
template<class T>
sparse_matrix<T>            vertcat(const band_matrix<T>& A, const sparse_matrix<T>& B);
template<class T>
sparse_matrix<T>            vertcat(const band_matrix<T>& A, const band_matrix<T>& B);
template<class T>
sparse_matrix<T>            vertcat(std::initializer_list<band_matrix<T>> mat_list);
template<class T>
sparse_matrix<T>            vertcat(const std::vector<band_matrix<T>>& mat_list);

// create block diagonal matrix with given blocks on the main diagonal
template<class T>
dense_matrix<T>             blkdiag(const dense_matrix<T>& A, const dense_matrix<T>& B);
template<class T>
dense_matrix<T>             blkdiag(std::initializer_list<dense_matrix<T>> mat_list);
template<class T>
dense_matrix<T>             blkdiag(const std::vector<dense_matrix<T>>& mat_list);
template<class T>
sparse_matrix<T>            blkdiag(const sparse_matrix<T>& A, const sparse_matrix<T>& B);
template<class T>
sparse_matrix<T>            blkdiag(std::initializer_list<sparse_matrix<T>> mat_list);
template<class T>
sparse_matrix<T>            blkdiag(const std::vector<sparse_matrix<T>>& mat_list);
template<class T>
sparse_matrix<T>            blkdiag(const sparse_matrix<T>& A, const band_matrix<T>& B);
template<class T>
sparse_matrix<T>            blkdiag(const band_matrix<T>& A, const sparse_matrix<T>& B);
template<class T>
sparse_matrix<T>            blkdiag(const band_matrix<T>& A, const band_matrix<T>& B);
template<class T>
sparse_matrix<T>            blkdiag(std::initializer_list<band_matrix<T>> mat_list);
template<class T>
sparse_matrix<T>            blkdiag(const std::vector<band_matrix<T>>& mat_list);

// form a block matrix of size m by n, with a copy of matrix A as each element.
template<class T>
dense_matrix<T>             repmat(const dense_matrix<T>& A, Integer m, Integer n);
template<class T>
sparse_matrix<T>            repmat(const sparse_matrix<T>& A, Integer m, Integer n);
template<class T>
sparse_matrix<T>            repmat(const band_matrix<T>& A, Integer m, Integer n);

// convers matrix m into dense matrix.
template<class T>
dense_matrix<T>             full(const dense_matrix<T>& m);
template<class T>
dense_matrix<T>             full(const sparse_matrix<T>& m);
template<class T>
dense_matrix<T>             full(const band_matrix<T>& m);

// convers matrix m into sparse matrix.
template<class T>
sparse_matrix<T>            sparse(const dense_matrix<T>& m);
template<class T>
sparse_matrix<T>            sparse(const sparse_matrix<T>& m);
template<class T>
sparse_matrix<T>            sparse(const band_matrix<T>& m);

// convers matrix m into band matrix.
template<class T>
band_matrix<T>              band(const dense_matrix<T>& m);
template<class T>
band_matrix<T>              band(const sparse_matrix<T>& m);
template<class T>
band_matrix<T>              band(const band_matrix<T>& m);

// make independent copy
template<class T>
dense_matrix<T>             clone(const dense_matrix<T>& m);
template<class T>
sparse_matrix<T>            clone(const sparse_matrix<T>& m);
template<class T>
band_matrix<T>              clone(const band_matrix<T>& m);

// transpose given matrix m.
template<class T>
dense_matrix<T>             trans(const dense_matrix<T>& m);
template<class T>
sparse_matrix<T>            trans(const sparse_matrix<T>& m);
template<class T>
band_matrix<T>              trans(const band_matrix<T>& m);

// makes conjugate transpose of given m.
template<class T>
dense_matrix<T>             ctrans(const dense_matrix<T>& m);
template<class T>
sparse_matrix<T>            ctrans(const sparse_matrix<T>& m);
template<class T>
band_matrix<T>              ctrans(const band_matrix<T>& m);

// transpose or conjugate transpose of given matrix m.
template<class T>
dense_matrix<T>             trans(const dense_matrix<T>& m, trans_type t);
template<class T>
sparse_matrix<T>            trans(const sparse_matrix<T>& m, trans_type t);
template<class T>
band_matrix<T>              trans(const band_matrix<T>& m, trans_type t);

// transpose or conjugate transpose of given matrix m.
template<class T>
dense_matrix<T>             trans(const dense_matrix<T>& m, trans_type_ext t);
template<class T>
sparse_matrix<T>            trans(const sparse_matrix<T>& m, trans_type_ext t);
template<class T>
band_matrix<T>              trans(const band_matrix<T>& m, trans_type_ext t);

// arranges matrix columns from left to right into vector placing 
// them from top to bottom.
template<class T>
dense_matrix<T>             vec(const dense_matrix<T>& m);
template<class T>
sparse_matrix<T>            vec(const sparse_matrix<T>& m);
template<class T>
sparse_matrix<T>            vec(const band_matrix<T>& m);

// returns a new matrix formed by extracting the lower triangular part
// of the matrix m, and settings all other elements to zero;
// the second argument is optional, and specifies how many diagonals 
// above the diagonal d should also be set to zero (d = 0 is the main diagonal,
// d > 0 are superdiagonals, d < 0 are subdiagonals)
template<class T>
dense_matrix<T>             tril(const dense_matrix<T>& m, Integer d = 0);
template<class T>
sparse_matrix<T>            tril(const sparse_matrix<T>& m, Integer d = 0);
template<class T>
band_matrix<T>              tril(const band_matrix<T>& m, Integer d = 0);

template<class T>
dense_matrix<T>             tril(dense_matrix<T>&& m, Integer d = 0);
template<class T>
sparse_matrix<T>            tril(sparse_matrix<T>&& m, Integer d = 0);
template<class T>
band_matrix<T>              tril(band_matrix<T>&& m, Integer d = 0);

// returns a new matrix formed by extracting the upper triangular part of
// the matrix m, and settings all other elements to zero; the second argument
// is optional, and specifies how many diagonals below the diagonal d should 
// also be set to zero (d = 0 is the main diagonal, d > 0 are superdiagonals, 
// d < 0 are subdiagonals)
template<class T>
dense_matrix<T>             triu(const dense_matrix<T>& m, Integer d = 0);
template<class T>
sparse_matrix<T>            triu(const sparse_matrix<T>& m, Integer d = 0);
template<class T>
band_matrix<T>              triu(const band_matrix<T>& m, Integer d = 0);

template<class T>
dense_matrix<T>             triu(dense_matrix<T>&& m, Integer d = 0);
template<class T>
sparse_matrix<T>            triu(sparse_matrix<T>&& m, Integer d = 0);
template<class T>
band_matrix<T>              triu(band_matrix<T>&& m, Integer d = 0);

// returns a new matrix formed by extracting the upper ud diagonals and lower
// ld diagonals of the matrix m, and settings all other elements to zero
template<class T>
band_matrix<T>              select_band(const dense_matrix<T> &m, Integer ld, Integer ud);
template<class T>
band_matrix<T>              select_band(const sparse_matrix<T> &m, Integer ld, Integer ud);
template<class T>
band_matrix<T>              select_band(const band_matrix<T> &m, Integer ld, Integer ud);

// return a copy of m with the order of the rows reversed.
template<class T>
dense_matrix<T>             flipud(const dense_matrix<T>& m);
template<class T>
sparse_matrix<T>            flipud(const sparse_matrix<T>& m);
template<class T>
sparse_matrix<T>            flipud(const band_matrix<T>& m);

// return a copy of m with the order of the columns reversed.
template<class T>
dense_matrix<T>             fliplr(const dense_matrix<T>& m);
template<class T>
sparse_matrix<T>            fliplr(const sparse_matrix<T>& m);
template<class T>
sparse_matrix<T>            fliplr(const band_matrix<T>& m);

// return a matrix with the specified dimensions m x n.
template<class T>
dense_matrix<T>             reshape(const dense_matrix<T>& A, Integer m, Integer n);
template<class T>
sparse_matrix<T>            reshape(const sparse_matrix<T>& A, Integer m, Integer n);
template<class T>
sparse_matrix<T>            reshape(const band_matrix<T>& A, Integer m, Integer n);

// return vector of elements on diagonal d; rhe second argument is optional, 
// and specifies which diagonal is to be returned (d = 0 is the main diagonal, 
// d > 0 are superdiagonals, d < 0 are subdiagonals)
template<class T>
dense_matrix<T>             get_diag(const dense_matrix<T>& m, Integer d = 0);
template<class T>
dense_matrix<T>             get_diag(const sparse_matrix<T>& m, Integer d = 0);
template<class T>
dense_matrix<T>             get_diag(const band_matrix<T>& m, Integer d = 0);

// counterclockwise rotation of matrix by n x 90 degree
template<class T>
dense_matrix<T>             rot90(const dense_matrix<T>& m, Integer n = 1);
template<class T>
sparse_matrix<T>            rot90(const sparse_matrix<T>& m, Integer n = 1);
template<class T>
sparse_matrix<T>            rot90(const band_matrix<T>& m, Integer n = 1);

// search for the non-zero elements in v; return indices of nonzero
// elements
template<class T>
dense_matrix<Integer>       find(const dense_matrix<T>& v);
template<class T>
dense_matrix<Integer>       find(const sparse_matrix<T>& v);
template<class T>
dense_matrix<Integer>       find(const band_matrix<T>& v);

// searches for the non-zero elements in v, using a user-defined test function;
// return indices of nonzero elements
template<class T>
dense_matrix<Integer>       find(const dense_matrix<T>& v, const test_function& t);
template<class T>
dense_matrix<Integer>       find(const sparse_matrix<T>& v, const test_function& t);
template<class T>
dense_matrix<Integer>       find(const band_matrix<T>& v, const test_function& t);

// version with test function, that require implementing test for type T only
template<class T>
dense_matrix<Integer>       find(const dense_matrix<T>& v, const test_type_function<T>& t);
template<class T>
dense_matrix<Integer>       find(const sparse_matrix<T>& v, const test_type_function<T>& t);
template<class T>
dense_matrix<Integer>       find(const band_matrix<T>& v, const test_type_function<T>& t);

// search for the non-zero elements in v; return a 2-tupple containing row and 
//column indices of nonzero elements
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const dense_matrix<T>& v);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const sparse_matrix<T>& v);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const band_matrix<T>& v);

// search for the non-zero elements in v, using a user-defined test function
// return a 2-tupple containing row and column indices of nonzero elements
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const dense_matrix<T>& v,const test_function& t);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const sparse_matrix<T>& v,const test_function& t);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const band_matrix<T>& v,const test_function& t);

// version with test function, that require implementing test for type T only
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const dense_matrix<T>& v,const test_type_function<T>& t);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const sparse_matrix<T>& v,const test_type_function<T>& t);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>>
                            find2(const band_matrix<T>& v,const test_type_function<T>& t);

// search for the non-zero elements in v; return a 3-tupple containing vectors
// of row (via .get<1>()), column (via .get<2>()) indices of non-zero elements 
// and the non-zero elements (via .get<3>())
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const dense_matrix<T>& v);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const sparse_matrix<T>& v);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const band_matrix<T>& v);

// search for the non-zero elements in v, using a user-defined test function
// return a 3-tupple containing vectors of row (via .get<1>()), column 
// (via .get<2>()) indices of non-zero elements  and the non-zero elements 
// (via .get<3>())
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const dense_matrix<T>& v, const test_function& t);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const sparse_matrix<T>& v, const test_function& t);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const band_matrix<T>& v, const test_function& t);

// version with test function, that require implementing test for type T only
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const dense_matrix<T>& v, const test_type_function<T>& t);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const sparse_matrix<T>& v, const test_type_function<T>& t);
template<class T>
tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
                            find3(const band_matrix<T>& v, const test_type_function<T>& t);

// sort the matrix v along dimension dim in asceding or desceding order;
// complex numbers are sorted by real part and then by imaginary part
template<class T>
dense_matrix<T>             sort(const dense_matrix<T>& v, int dim = 1, bool asceding = true);
template<class T>
sparse_matrix<T>            sort(const sparse_matrix<T>& v, int dim = 1, bool asceding = true);
template<class T>
sparse_matrix<T>            sort(const band_matrix<T>& v, int dim = 1, bool asceding = true);

// sort the matrix v along dimension dim in asceding or desceding order;
// return a 2-tupple with the sorted matrix (via .get<1>()) and a matrix 
// with the indices after sorting permutation (via .get<2>())
// complex numbers are sorted by real part and then by imaginary part
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            sort2(const dense_matrix<T>& v, int dim = 1, bool asceding = true);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sort2(const sparse_matrix<T>& v, int dim = 1, bool asceding = true);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sort2(const band_matrix<T>& v, int dim = 1, bool asceding = true);

// sort rows of the matrix v as a group;
// complex numbers are sorted by real part and then by imaginary part
template<class T>
dense_matrix<T>             sortrows(const dense_matrix<T>& v);
template<class T>
sparse_matrix<T>            sortrows(const sparse_matrix<T>& v);
template<class T>
sparse_matrix<T>            sortrows(const band_matrix<T>& v);

// sort rows of the matrix v as a group; return a 2-tupple with the sorted 
// matrix (via .get<1>()) and a vector of indices of rows after sorting 
// permutation (via .get<2>()); complex numbers are sorted by real part 
// and then by imaginary part
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            sortrows2(const dense_matrix<T>& v);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sortrows2(const sparse_matrix<T>& v);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sortrows2(const band_matrix<T>& v);

// sort rows of the matrix v as a group considering only elements in columns
// cols; do ascending sort when column number is positive, descending if negative;
// for example if cols = [2, -1] then the function use second column for primary
// comparison (ascending) and first column for secondary comparison (descending);
// complex numbers are sorted by real part and then by imaginary part
template<class T>
dense_matrix<T>             sortrows(const dense_matrix<T>& v, const Matrix& cols);
template<class T>
sparse_matrix<T>            sortrows(const sparse_matrix<T>& v, const Matrix& cols);
template<class T>
sparse_matrix<T>            sortrows(const band_matrix<T>& v, const Matrix& cols);

// sort rows of the matrix v as a group considering only elements in columns
// cols; do ascending sort when column number is positive, descending if negative;
// for example if cols = [2, -1] then the function use second column for primary
// comparison (ascending) and first column for secondary comparison (descending);
// return a 2-tupple with the sorted matrix (via .get<1>()) and a vector of indices
// of rows after sorting permutation (via .get<2>());
// complex numbers are sorted by real part and then by imaginary part
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            sortrows2(const dense_matrix<T>& v, const Matrix& cols);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sortrows2(const sparse_matrix<T>& v, const Matrix& cols);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sortrows2(const band_matrix<T>& v, const Matrix& cols);

// sort columns of the matrix v as a group;
// complex numbers are sorted by real part and then by imaginary part
template<class T>
dense_matrix<T>             sortcols(const dense_matrix<T>& v);
template<class T>
sparse_matrix<T>            sortcols(const sparse_matrix<T>& v);
template<class T>
sparse_matrix<T>            sortcols(const band_matrix<T>& v);

// sort columns of the matrix v as a group; return a 2-tupple with the sorted 
// matrix (via .get<1>()) and a vector of indices of colums after sorting 
// permutation (via .get<2>()); complex numbers are sorted by real part 
// and then by imaginary part
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            sortcols2(const dense_matrix<T>& v);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sortcols2(const sparse_matrix<T>& v);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sortcols2(const band_matrix<T>& v);

// sort columns of the matrix v as a group considering only elements in rows
// 'rows'; do ascending sort when row number is positive, descending if negative;
// for example if rows = [2, -1] then the function use second row for primary
// comparison (ascending) and first row for secondary comparison (descending);
// complex numbers are sorted by real part and then by imaginary part
template<class T>
dense_matrix<T>             sortcols(const dense_matrix<T>& v, const Matrix& dims);
template<class T>
sparse_matrix<T>            sortcols(const sparse_matrix<T>& v, const Matrix& dims);
template<class T>
sparse_matrix<T>            sortcols(const band_matrix<T>& v, const Matrix& dims);

// sort columns of the matrix v as a group considering only elements in rows
// 'rows'; do ascending sort when row number is positive, descending if negative;
// for example if rows = [2, -1] then the function use second row for primary
// comparison (ascending) and first row for secondary comparison (descending);
// return a 2-tupple with the sorted matrix (via .get<1>()) and a vector of indices
// of columns after sorting permutation (via .get<2>());
// complex numbers are sorted by real part and then by imaginary part
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            sortcols2(const dense_matrix<T>& v, const Matrix& dims);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sortcols2(const sparse_matrix<T>& v, const Matrix& dims);
template<class T>
tuple<sparse_matrix<T>, dense_matrix<Integer>>
                            sortcols2(const band_matrix<T>& v, const Matrix& dims);

// check if a matrix v is sorted along dimension dim; for dim = 1 issorted(v) 
// returns a boolean vector, specifying if columns of v are sorted; for dim = 2
// function does its job for transposed matrix;
// complex numbers are sorted by real part and then by imaginary part
template<class T>
dense_matrix<Integer>       issorted(const dense_matrix<T>& v, int dim = 1, bool asceding = true);
template<class T>
dense_matrix<Integer>       issorted(const sparse_matrix<T>& v, int dim = 1, bool asceding = true);
template<class T>
dense_matrix<Integer>       issorted(const band_matrix<T>& v, int dim = 1, bool asceding = true);

// remove elements v in a sparse matrix A satisfying abs(v) <= tol.
// do nothing for dense and band matrices
template<class T>
dense_matrix<T>             drop_sparse(const dense_matrix<T>& A, Real tol);
template<class T>
sparse_matrix<T>            drop_sparse(const sparse_matrix<T>& A, Real tol);
template<class T>
band_matrix<T>              drop_sparse(const band_matrix<T>& A, Real tol);

// convert the matrix to a new type given by new_type
template<class New_mat, class T>
New_mat                     convert(const dense_matrix<T>& A);
template<class New_mat, class T>
New_mat                     convert(const sparse_matrix<T>& A);
template<class New_mat, class T>
New_mat                     convert(const band_matrix<T>& A);

// convert value type of stored elements not changing matrix structure
template<class S, class T>
dense_matrix<S>             convert_value(const dense_matrix<T>& A);
template<class S, class T>
sparse_matrix<S>            convert_value(const sparse_matrix<T>& A);
template<class S, class T>
band_matrix<S>              convert_value(const band_matrix<T>& A);

// convert objects stored in the matrix to the type ti; if elements stored
// in the matrix are not of Object type, then original matrix is returned
template<class T>
dense_matrix<T>             convert_object(const dense_matrix<T>& A, ti::ti_object ti);
template<class T>
sparse_matrix<T>            convert_object(const sparse_matrix<T>& A, ti::ti_object ti);
template<class T>
band_matrix<T>              convert_object(const band_matrix<T>& A, ti::ti_object ti);

//--------------------------------------------------------------
//      functions operating on rows or columns
//--------------------------------------------------------------

// number of nonzero elements along dimension dim
// return Integer matrix of size 1 x A.cols() if dim == 1
// or Integer matrix of size A.rows() x 1 if dim == 2
template<class T>
dense_matrix<Integer>       nnz(const dense_matrix<T>& A, Integer dim);
template<class T>
dense_matrix<Integer>       nnz(const sparse_matrix<T>& A, Integer dim);
template<class T>
dense_matrix<Integer>       nnz(const band_matrix<T>& A, Integer dim);

// check if all elements in a matrix are nonzero; if dim == 1 then
// function operates on columns and return Integer matrix of 
// size 1 x A.cols(); if dim == 2 then function operators on rows
// and return Integer matrix of size A.rows() x 1
template<class T>
dense_matrix<Integer>       all(const dense_matrix<T>& v,int dim = 1);
template<class T>
dense_matrix<Integer>       all(const sparse_matrix<T>& v,int dim = 1);
template<class T>
dense_matrix<Integer>       all(const band_matrix<T>& v,int dim = 1);

// check if all elements in a matrix pass the test function t; 
// if dim == 1 then function operates on columns and return
// Integer matrix of size 1 x A.cols(); if dim == 2 then function
// operators on rows and return Integer matrix of size A.rows() x 1
template<class T>
dense_matrix<Integer>       all(const dense_matrix<T>& v,const test_function& t, int dim = 1);
template<class T>
dense_matrix<Integer>       all(const sparse_matrix<T>& v,const test_function& t, int dim = 1);
template<class T>
dense_matrix<Integer>       all(const band_matrix<T>& v,const test_function& t, int dim = 1);

// version with test function, that require implementing test for type T only
template<class T>
dense_matrix<Integer>       all(const dense_matrix<T>& v, const test_type_function<T>& t, int dim = 1);
template<class T>
dense_matrix<Integer>       all(const sparse_matrix<T>& v, const test_type_function<T>& t, int dim = 1);
template<class T>
dense_matrix<Integer>       all(const band_matrix<T>& v, const test_type_function<T>& t, int dim = 1);

// check if any element in a matrix is nonzero; if dim == 1 then
// function operates on columns and return Integer matrix of 
// size 1 x A.cols(); if dim == 2 then function operators on rows
// and return Integer matrix of size A.rows() x 1
template<class T>
dense_matrix<Integer>       any(const dense_matrix<T>& v,int dim = 1);
template<class T>
dense_matrix<Integer>       any(const sparse_matrix<T>& v,int dim = 1);
template<class T>
dense_matrix<Integer>       any(const band_matrix<T>& v,int dim = 1);

// check if any element in a matrix pass the test function t; 
// if dim == 1 then function operates on columns and return
// Integer matrix of size 1 x A.cols(); if dim == 2 then function
// operators on rows and return Integer matrix of size A.rows() x 1
template<class T>
dense_matrix<Integer>       any(const dense_matrix<T>& v,const test_function& t, int dim = 1);
template<class T>
dense_matrix<Integer>       any(const sparse_matrix<T>& v,const test_function& t, int dim = 1);
template<class T>
dense_matrix<Integer>       any(const band_matrix<T>& v,const test_function& t, int dim = 1);

// version with test function, that require implementing test for type T only
template<class T>
dense_matrix<Integer>       any(const dense_matrix<T>& v, const test_type_function<T>& t, int dim = 1);
template<class T>
dense_matrix<Integer>       any(const sparse_matrix<T>& v, const test_type_function<T>& t, int dim = 1);
template<class T>
dense_matrix<Integer>       any(const band_matrix<T>& v, const test_type_function<T>& t, int dim = 1);

// sum of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
dense_matrix<T>             sum(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             sum(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             sum(const band_matrix<T>& v, int dim = 1);

// product of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
dense_matrix<T>             prod(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             prod(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             prod(const band_matrix<T>& v, int dim = 1);

// cumulative sum of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
dense_matrix<T>             cumsum(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             cumsum(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             cumsum(const band_matrix<T>& v, int dim = 1);

// cumulative product of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
dense_matrix<T>             cumprod(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             cumprod(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             cumprod(const band_matrix<T>& v, int dim = 1);

// mean value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1
template<class T>
dense_matrix<typename md::unify_types<T,Float>::type>
                            mean(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<typename md::unify_types<T,Float>::type>
                            mean(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<typename md::unify_types<T,Float>::type>
                            mean(const band_matrix<T>& v, int dim = 1);

// standard deviation of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// if unbiaded == true, then unbiased estimator of standard deviation
// is returned (i.e. normalization by N-1 if N > 1, where N is the sample
// size), otherwise biased estimator is returned (normalization by N)
template<class T>
dense_matrix<typename md::real_type_int_real<T>::type>
                            std(const dense_matrix<T>& v, int dim = 1, bool unbiased = true);
template<class T>
dense_matrix<typename md::real_type_int_real<T>::type>
                            std(const sparse_matrix<T>& v, int dim = 1, bool unbiased = true);
template<class T>
dense_matrix<typename md::real_type_int_real<T>::type>
                            std(const band_matrix<T>& v, int dim = 1, bool unbiased = true);

// minimum value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
template<class T>
dense_matrix<T>             min_d(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             min_d(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             min_d(const band_matrix<T>& v, int dim = 1);

// maximum value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
template<class T>
dense_matrix<T>             max_d(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             max_d(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<T>             max_d(const band_matrix<T>& v, int dim = 1);

// minimum absolute value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
template<class T>
dense_matrix<typename md::real_type<T>::type>
                            min_abs_d(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<typename md::real_type<T>::type>
                            min_abs_d(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<typename md::real_type<T>::type>
                            min_abs_d(const band_matrix<T>& v, int dim = 1);

// maximum absolute value of elements in a matrix; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
template<class T>
dense_matrix<typename md::real_type<T>::type>
                            max_abs_d(const dense_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<typename md::real_type<T>::type>
                            max_abs_d(const sparse_matrix<T>& v, int dim = 1);
template<class T>
dense_matrix<typename md::real_type<T>::type>
                            max_abs_d(const band_matrix<T>& v, int dim = 1);

// minimum value of elements in a matrix and positions of minimum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// return a 2-tupple of minimum values (via .get<1>()) and a 
// indices (via .get<2>())
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            min2(const dense_matrix<T>& v, int dim = 1);
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            min2(const sparse_matrix<T>& v, int dim = 1);
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            min2(const band_matrix<T>& v, int dim = 1);

// maximum value of elements in a matrix and positions of maximum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// return a 2-tupple of maximum values (via .get<1>()) and a 
// indices (via .get<2>())
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            max2(const dense_matrix<T>& v, int dim = 1);
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            max2(const sparse_matrix<T>& v, int dim = 1);
template<class T>
tuple<dense_matrix<T>, dense_matrix<Integer>>
                            max2(const band_matrix<T>& v, int dim = 1);

// minimum absolute value of elements in a matrix and positions of minimum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// return a 2-tupple of minimum values (via .get<1>()) and a 
// indices (via .get<2>())
template<class T>
tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
                            min_abs2(const dense_matrix<T>& v, int dim = 1);
template<class T>
tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
                            min_abs2(const sparse_matrix<T>& v, int dim = 1);
template<class T>
tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
                            min_abs2(const band_matrix<T>& v, int dim = 1);

// maximum absolute value of elements in a matrix and positions of maximum values; 
// if dim == 1 then function operates on columns and return a matrix 
// of size 1 x A.cols(); if dim == 2 then function operators on rows 
// and return a matrix of size A.rows() x 1;
// return a 2-tupple of maximum values (via .get<1>()) and a 
// indices (via .get<2>())
template<class T>
tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
                            max_abs2(const dense_matrix<T>& v, int dim = 1);
template<class T>
tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
                            max_abs2(const sparse_matrix<T>& v, int dim = 1);
template<class T>
tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
                            max_abs2(const band_matrix<T>& v, int dim = 1);

// function equivalent to min_d
template<class T>
inline dense_matrix<T>      min(const dense_matrix<T>& v)      { return min_d(v,1); };
template<class T>
inline dense_matrix<T>      min(const sparse_matrix<T>& v)     { return min_d(v,1); };
template<class T>
inline dense_matrix<T>      min(const band_matrix<T>& v)       { return min_d(v,1); };

// function equivalent to max_d
template<class T>
inline dense_matrix<T>      max(const dense_matrix<T>& v)      { return max_d(v,1); };
template<class T>
inline dense_matrix<T>      max(const sparse_matrix<T>& v)     { return max_d(v,1); };
template<class T>
inline dense_matrix<T>      max(const band_matrix<T>& v)       { return max_d(v,1); };

// function equivalent to min_abs_d
template<class T>
inline dense_matrix<typename md::real_type<T>::type>
                            min_abs(const dense_matrix<T>& v)  { return min_abs_d(v,1); };
template<class T>
inline dense_matrix<typename md::real_type<T>::type>
                            min_abs(const sparse_matrix<T>& v) { return min_abs_d(v,1); };
template<class T>
inline dense_matrix<typename md::real_type<T>::type>
                            min_abs(const band_matrix<T>& v)   { return min_abs_d(v,1); };

// function equivalent to max_abs_d
template<class T>
inline dense_matrix<typename md::real_type<T>::type>
                            max_abs(const dense_matrix<T>& v)  { return max_abs_d(v,1); };
template<class T>
inline dense_matrix<typename md::real_type<T>::type>
                            max_abs(const sparse_matrix<T>& v) { return max_abs_d(v,1); };
template<class T>
inline dense_matrix<typename md::real_type<T>::type>
                            max_abs(const band_matrix<T>& v)   { return max_abs_d(v,1); };

//--------------------------------------------------------------
//      functions operating on all elements
//--------------------------------------------------------------

// number of nonzero elements;
template<class T>
Integer                     nnz_vec(const dense_matrix<T>& A);
template<class T>
Integer                     nnz_vec(const sparse_matrix<T>& A);
template<class T>
Integer                     nnz_vec(const band_matrix<T>& A);

// check if all elements in a matrix are nonzero; 
// return true for empty matrix
template<class T>
bool                        all_vec(const dense_matrix<T>& v);
template<class T>
bool                        all_vec(const sparse_matrix<T>& v);
template<class T>
bool                        all_vec(const band_matrix<T>& v);

// check if all elements in a matrix pass the test function t; 
// return true for empty matrix
template<class T>
bool                        all_vec(const dense_matrix<T>& v,const test_function& t);
template<class T>
bool                        all_vec(const sparse_matrix<T>& v,const test_function& t);
template<class T>
bool                        all_vec(const band_matrix<T>& v,const test_function& t);

// version with test function, that require implementing test for type T only
template<class T>
bool                        all_vec(const dense_matrix<T>& v, const test_type_function<T>& t);
template<class T>
bool                        all_vec(const sparse_matrix<T>& v, const test_type_function<T>& t);
template<class T>
bool                        all_vec(const band_matrix<T>& v, const test_type_function<T>& t);

// check if any element in a matrix is nonzero;
// return false for empty matrix
template<class T>
bool                        any_vec(const dense_matrix<T>& v);
template<class T>
bool                        any_vec(const sparse_matrix<T>& v);
template<class T>
bool                        any_vec(const band_matrix<T>& v);

// check if any element in a matrix pass the test function t; 
// return false for empty matrix
template<class T>
bool                        any_vec(const dense_matrix<T>& v,const test_function& t);
template<class T>
bool                        any_vec(const sparse_matrix<T>& v,const test_function& t);
template<class T>
bool                        any_vec(const band_matrix<T>& v,const test_function& t);

// version with test function, that require implementing test for type T only
template<class T>
bool                        any_vec(const dense_matrix<T>& v, const test_type_function<T>& t);
template<class T>
bool                        any_vec(const sparse_matrix<T>& v, const test_type_function<T>& t);
template<class T>
bool                        any_vec(const band_matrix<T>& v, const test_type_function<T>& t);

// sum of elements in a matrix; return scalar
// return zero for empty matrix
template<class T>
T                           sum_vec(const dense_matrix<T>& v);
template<class T>
T                           sum_vec(const sparse_matrix<T>& v);
template<class T>
T                           sum_vec(const band_matrix<T>& v);

// product of elements in a matrix; return scalar
// return zero for empty matrix
template<class T>
T                           prod_vec(const dense_matrix<T>& v);
template<class T>
T                           prod_vec(const sparse_matrix<T>& v);
template<class T>
T                           prod_vec(const band_matrix<T>& v);

// mean value of elements in a matrix; return scalar
// return zero for empty matrix
template<class T>
typename md::unify_types<T,Float>::type
                            mean_vec(const dense_matrix<T>& v);
template<class T>
typename md::unify_types<T,Float>::type
                            mean_vec(const sparse_matrix<T>& v);
template<class T>
typename md::unify_types<T,Float>::type
                            mean_vec(const band_matrix<T>& v);

// standard deviation of elements in a matrix; return scalar
// if unbiaded == true, then unbiased estimator of standard deviation
// is returned (i.e. normalization by N-1 if N > 1, where N is the sample
// size), otherwise biased estimator is returned (normalization by N)
template<class T>
typename md::real_type_int_real<T>::type
                            std_vec(const dense_matrix<T>& v, bool unbiased = true);
template<class T>
typename md::real_type_int_real<T>::type
                            std_vec(const sparse_matrix<T>& v, bool unbiased = true);
template<class T>
typename md::real_type_int_real<T>::type
                            std_vec(const band_matrix<T>& v, bool unbiased = true);

// minimum value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
template<class T>
T                           min_vec(const dense_matrix<T>& v);
template<class T>
T                           min_vec(const sparse_matrix<T>& v);
template<class T>
T                           min_vec(const band_matrix<T>& v);

// maximum value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
template<class T>
T                           max_vec(const dense_matrix<T>& v);
template<class T>
T                           max_vec(const sparse_matrix<T>& v);
template<class T>
T                           max_vec(const band_matrix<T>& v);

// minimum absolute value of elements in a matrix; return scalar
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
template<class T>
typename md::real_type<T>::type
                            min_abs_vec(const dense_matrix<T>& v);
template<class T>
typename md::real_type<T>::type
                            min_abs_vec(const sparse_matrix<T>& v);
template<class T>
typename md::real_type<T>::type
                            min_abs_vec(const band_matrix<T>& v);

// maximum absolute value of elements in a matrix; 
// complex numbers are compared by real part and then by imaginary part;
// return zero for empty matrix
template<class T>
typename md::real_type<T>::type
                            max_abs_vec(const dense_matrix<T>& v);
template<class T>
typename md::real_type<T>::type
                            max_abs_vec(const sparse_matrix<T>& v);
template<class T>
typename md::real_type<T>::type
                            max_abs_vec(const band_matrix<T>& v);

// minimum value of elements in a matrix and positions of minimum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of minimum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
template<class T>
tuple<T, Integer, Integer > min2_vec(const dense_matrix<T>& v);
template<class T>
tuple<T, Integer, Integer > min2_vec(const sparse_matrix<T>& v);
template<class T>
tuple<T, Integer, Integer > min2_vec(const band_matrix<T>& v);

// maximum value of elements in a matrix and positions of maximum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of maximum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
template<class T>
tuple<T, Integer, Integer > max2_vec(const dense_matrix<T>& v);
template<class T>
tuple<T, Integer, Integer > max2_vec(const sparse_matrix<T>& v);
template<class T>
tuple<T, Integer, Integer > max2_vec(const band_matrix<T>& v);

// minimum absolute value of elements in a matrix and positions of minimum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of minimum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
template<class T>
tuple<typename md::real_type<T>::type, Integer, Integer >
                            min_abs2_vec(const dense_matrix<T>& v);
template<class T>
tuple<typename md::real_type<T>::type, Integer, Integer >
                            min_abs2_vec(const sparse_matrix<T>& v);
template<class T>
tuple<typename md::real_type<T>::type, Integer, Integer >
                            min_abs2_vec(const band_matrix<T>& v);

// maximum absolute value of elements in a matrix and positions of maximum values; 
// complex numbers are compared by real part and then by imaginary part;
// return a 3-tupple of maximum value (via .get<1>()) and a 
// row and a column index (via .get<2>() and .get<3>())
// return mat_int2(0,0,0) for empty matrix
template<class T>
tuple<typename md::real_type<T>::type, Integer, Integer >
                            max_abs2_vec(const dense_matrix<T>& v);
template<class T>
tuple<typename md::real_type<T>::type, Integer, Integer >
                            max_abs2_vec(const sparse_matrix<T>& v);
template<class T>
tuple<typename md::real_type<T>::type, Integer, Integer >
                            max_abs2_vec(const band_matrix<T>& v);
};

#include "matcl-matrep/details/matrix_rep_functions.inl"
