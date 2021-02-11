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

#include "matcl-matrep/matrix/matrix_rep_functions.h"
#include "matcl-matrep/matrix/matrix_rep_dense.h"
#include "matcl-matrep/IO/matrix_io.h"
#include "matcl-matrep/lib_functions/vecfunc.h"

namespace matcl
{

//-----------------------------------------------------------------
//                  IO functions
//-----------------------------------------------------------------
template<class T>
inline void
matcl::disp(const dense_matrix<T>& m, const disp_stream_ptr& os, const options& opts)
{
    return matcl::disp(matcl::Matrix(m), os, opts);
};

template<class T>
inline std::ostream& matcl::operator<<(std::ostream& os, const dense_matrix<T>& m)
{
    os << Matrix(m);
    return os;
};

template<class T>
inline std::istream& matcl::operator>>(std::istream& is, dense_matrix<T>& m)
{
    Matrix tmp;
    is >> tmp;
    m = dense_matrix<T>(tmp);
    return is;
};

template<class T>
inline void matcl::save(oarchive& ar,const dense_matrix<T>& mat)
{
    save(ar, matcl::Matrix(mat));
};

template<class T>
inline void matcl::load(iarchive& ar,dense_matrix<T>& mat)
{
    Matrix tmp;
    load(ar,tmp);
    mat = dense_matrix<T>(tmp);
};

template<class T>
inline void matcl::disp(const sparse_matrix<T>& m, const disp_stream_ptr& os, const options& opts)
{
    return disp(matcl::Matrix(m), os, opts);
};

template<class T>
inline std::ostream& matcl::operator<<(std::ostream& os, const sparse_matrix<T>& m)
{
    os << Matrix(m);
    return os;
};

template<class T>
inline std::istream& matcl::operator>>(std::istream& is, sparse_matrix<T>& m)
{
    Matrix tmp;
    is >> tmp;
    m = sparse_matrix<T>(tmp);
    return is;
};

template<class T>
inline void matcl::save(oarchive& ar,const sparse_matrix<T>& mat)
{
    save(ar, matcl::Matrix(mat));
};

template<class T>
inline void matcl::load(iarchive& ar,sparse_matrix<T>& mat)
{
    Matrix tmp;
    load(ar,tmp);
    mat = sparse_matrix<T>(tmp);
};

template<class T>
inline void matcl::disp(const band_matrix<T>& m, const disp_stream_ptr& os, const options& opts)
{
    return disp(matcl::Matrix(m), os, opts);
};

template<class T>
inline std::ostream& matcl::operator<<(std::ostream& os, const band_matrix<T>& m)
{
    os << Matrix(m);
    return os;
};

template<class T>
inline std::istream& matcl::operator>>(std::istream& is, band_matrix<T>& m)
{
    Matrix tmp;
    is >> tmp;
    m = band_matrix<T>(tmp);
    return is;
};

template<class T>
inline void matcl::save(oarchive& ar,const band_matrix<T>& mat)
{
    save(ar, matcl::Matrix(mat));
};

template<class T>
inline void matcl::load(iarchive& ar,band_matrix<T>& mat)
{
    Matrix tmp;
    load(ar,tmp);
    mat = band_matrix<T>(tmp);
};

// -----------------------------------------------------------------------
//                 manip functions
// -----------------------------------------------------------------------
template<class T>
inline dense_matrix<T> matcl::delrows(const dense_matrix<T>& A, const colon& c)
{
    return A.delrows(c);
};

template<class T>
inline dense_matrix<T> matcl::delrows(dense_matrix<T>&& A, const colon& c)
{
    return std::move(A).delrows(c);
};

template<class T>
inline dense_matrix<T> matcl::delcols(const dense_matrix<T>& A, const colon& c)
{
    return A.delcols(c);
};

template<class T>
inline dense_matrix<T> matcl::delcols(dense_matrix<T>&& A, const colon& c)
{
    return std::move(A).delcols(c);
};

template<class T>
inline dense_matrix<T> matcl::delrowscols(const dense_matrix<T>& A, const colon& c1, const colon& c2)
{
    return A.delrowscols(c1,c2);
};

template<class T>
inline dense_matrix<T> matcl::delrowscols(dense_matrix<T>&& A, const colon& c1, const colon& c2)
{
    return std::move(A).delrowscols(c1,c2);
};

template<class T>
inline dense_matrix<T> matcl::horzcat(const dense_matrix<T>& A, const dense_matrix<T>& B)
{
    return dense_matrix<T>(horzcat(Matrix(A), Matrix(B)));
};

template<class T>
inline dense_matrix<T> matcl::horzcat(const std::vector<dense_matrix<T>>& mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return dense_matrix<T>(ret);
};

template<class T>
inline dense_matrix<T> matcl::horzcat(std::initializer_list<dense_matrix<T>> mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return dense_matrix<T>(ret);
};

template<class T>
inline dense_matrix<T> matcl::vertcat(const dense_matrix<T>& A, const dense_matrix<T>& B)
{
    return dense_matrix<T>(vertcat(Matrix(A), Matrix(B)));
};

template<class T>
inline dense_matrix<T> matcl::vertcat(const std::vector<dense_matrix<T>>& mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return dense_matrix<T>(ret);
};

template<class T>
inline dense_matrix<T> matcl::vertcat(std::initializer_list<dense_matrix<T>> mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return dense_matrix<T>(ret);
};

template<class T>
inline dense_matrix<T> matcl::blkdiag(const dense_matrix<T>& A, const dense_matrix<T>& B)
{
    return dense_matrix<T>(blkdiag(Matrix(A), Matrix(B)));
};

template<class T>
inline dense_matrix<T> matcl::blkdiag(const std::vector<dense_matrix<T>>& mat_list)
{
    std::vector<Matrix> vec_mat;
    for( auto it: mat_list)
        vec_mat.push_back(Matrix(it));

    Matrix ret = blkdiag(vec_mat);
    return dense_matrix<T>(ret);
};

template<class T>
inline dense_matrix<T> matcl::blkdiag(std::initializer_list<dense_matrix<T>> mat_list)
{
    std::vector<Matrix> vec_mat;
    for( auto it: mat_list)
        vec_mat.push_back(Matrix(it));

    Matrix ret = blkdiag(vec_mat);
    return dense_matrix<T>(ret);
};

template<class T>
inline dense_matrix<T> matcl::repmat(const dense_matrix<T>& A, Integer m, Integer n)
{
    return dense_matrix<T>(repmat(Matrix(A), m, n));
};

template<class T>
inline dense_matrix<T> matcl::full(const dense_matrix<T>& m)
{
    return m;
};

template<class T>
inline dense_matrix<T> matcl::full(const sparse_matrix<T>& m)
{
    return dense_matrix<T>(full(Matrix(m)));
};

template<class T>
inline dense_matrix<T> matcl::full(const band_matrix<T>& m)
{
    return dense_matrix<T>(full(Matrix(m)));
};

template<class T>
inline dense_matrix<T> matcl::clone(const dense_matrix<T>& m)
{
    return dense_matrix<T>(clone(Matrix(m)));
};

template<class T>
inline dense_matrix<T> matcl::trans(const dense_matrix<T>& m)
{
    return dense_matrix<T>(trans(Matrix(m)));
};

template<class T>
inline dense_matrix<T> matcl::ctrans(const dense_matrix<T>& m)
{
    return dense_matrix<T>(ctrans(Matrix(m)));
};

template<class T>
inline dense_matrix<T> matcl::trans(const dense_matrix<T>& m, trans_type t)
{
    return dense_matrix<T>(trans(Matrix(m), t));
};

template<class T>
inline dense_matrix<T> matcl::trans(const dense_matrix<T>& m, trans_type_ext t)
{
    return dense_matrix<T>(trans(Matrix(m), t));
};

template<class T>
inline dense_matrix<T> matcl::vec(const dense_matrix<T>& m)
{
    return dense_matrix<T>(vec(Matrix(m)));
};

template<class T>
inline dense_matrix<T> matcl::tril(const dense_matrix<T>& m, Integer d)
{
    return dense_matrix<T>(tril(Matrix(m),d));
};
template<class T>
inline dense_matrix<T> matcl::tril(dense_matrix<T>&& m, Integer d)
{
    return dense_matrix<T>(tril(Matrix(std::move(m)),d));
};

template<class T>
inline dense_matrix<T> matcl::triu(const dense_matrix<T>& m, Integer d)
{
    return dense_matrix<T>(triu(Matrix(m),d));
};

template<class T>
inline dense_matrix<T> matcl::triu(dense_matrix<T>&& m, Integer d)
{
    return dense_matrix<T>(triu(Matrix(std::move(m)),d));
};

template<class T>
inline band_matrix<T> matcl::select_band(const dense_matrix<T>& m, Integer ld, Integer ud)
{
    return band_matrix<T>(select_band(Matrix(m),ld, ud));
};

template<class T>
inline dense_matrix<T> matcl::flipud(const dense_matrix<T>& m)
{
    return dense_matrix<T>(flipud(Matrix(m)));
};

template<class T>
inline dense_matrix<T> matcl::fliplr(const dense_matrix<T>& m)
{
    return dense_matrix<T>(fliplr(Matrix(m)));
};

template<class T>
inline dense_matrix<T> matcl::reshape(const dense_matrix<T>& A, Integer m, Integer n)
{
    return dense_matrix<T>(reshape(Matrix(A),m,n));
};

template<class T>
inline dense_matrix<T> matcl::get_diag(const dense_matrix<T>& m, Integer d)
{
    return dense_matrix<T>(get_diag(Matrix(m),d));
};

template<class T>
inline dense_matrix<T> matcl::rot90(const dense_matrix<T>& m, Integer n)
{
    return dense_matrix<T>(rot90(Matrix(m),n));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const dense_matrix<T>& v)
{
    return dense_matrix<Integer>(find(Matrix(v)));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const dense_matrix<T>& v, const test_function& t)
{
    return dense_matrix<Integer>(find(Matrix(v), t));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const dense_matrix<T>& v, const test_type_function<T>& t)
{
    return dense_matrix<Integer>(find(Matrix(v), t.to_test_function()));
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const dense_matrix<T>& v)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v));

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const dense_matrix<T>& v,const test_function& t)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v),t);

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const dense_matrix<T>& v,const test_type_function<T>& t)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v),t.to_test_function());

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const dense_matrix<T>& v)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v));

    dense_matrix<Integer>   ret1(R);
    dense_matrix<Integer>   ret2(C);
    dense_matrix<T>         ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const dense_matrix<T>& v, const test_function& t)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v), t);

    dense_matrix<Integer>   ret1(R);
    dense_matrix<Integer>   ret2(C);
    dense_matrix<T>         ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const dense_matrix<T>& v, const test_type_function<T>& t)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v), t.to_test_function());

    dense_matrix<Integer>   ret1(R);
    dense_matrix<Integer>   ret2(C);
    dense_matrix<T>         ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline dense_matrix<T> matcl::sort(const dense_matrix<T>& v, int dim, bool asceding)
{
    return dense_matrix<T>(sort(Matrix(v), dim, asceding));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::sort2(const dense_matrix<T>& v, int dim, bool asceding)
{
    Matrix S, I;
    tie(S,I) = sort2(Matrix(v),dim,asceding);

    dense_matrix<T>         ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<dense_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline dense_matrix<T> matcl::sortrows(const dense_matrix<T>& v)
{
    return dense_matrix<T>(sortrows(Matrix(v)));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::sortrows2(const dense_matrix<T>& v)
{
    Matrix S, I;
    tie(S,I) = sort2(Matrix(v));

    dense_matrix<T>         ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<dense_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline dense_matrix<T> matcl::sortrows(const dense_matrix<T>& v, const Matrix& cols)
{
    return dense_matrix<T>(sortrows(Matrix(v), cols));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::sortrows2(const dense_matrix<T>& v, const Matrix& cols)
{
    Matrix S, I;
    tie(S,I) = sortrows2(Matrix(v),cols);

    dense_matrix<T>         ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<dense_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
}

template<class T>
inline dense_matrix<T> matcl::sortcols(const dense_matrix<T>& v)
{
    return dense_matrix<T>(sortcols(Matrix(v)));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::sortcols2(const dense_matrix<T>& v)
{
    Matrix S, I;
    tie(S,I) = sortcols2(Matrix(v));

    dense_matrix<T>         ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<dense_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline dense_matrix<T> matcl::sortcols(const dense_matrix<T>& v, const Matrix& dims)
{
    return dense_matrix<T>(sortcols(Matrix(v), dims));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::sortcols2(const dense_matrix<T>& v, const Matrix& dims)
{
    Matrix S, I;
    tie(S,I) = sortcols2(Matrix(v),dims);

    dense_matrix<T>         ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<dense_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline dense_matrix<Integer> matcl::issorted(const dense_matrix<T>& v, int dim, bool asceding)
{
    return dense_matrix<Integer>(issorted(Matrix(v),dim,asceding));
};

template<class T>
inline dense_matrix<T> matcl::drop_sparse(const dense_matrix<T>& A, Real tol)
{
    (void)tol;
    return A;
};

// convert the matrix to a new type given by new_type
template<class New_mat, class T>
inline New_mat matcl::convert(const dense_matrix<T>& A)
{
    return New_mat(A.to_matrix());
};

template<class S, class T>
inline dense_matrix<S> matcl::convert_value(const dense_matrix<T>& m)
{
    return dense_matrix<S>(convert_value(Matrix(m), matrix_traits::value_code<S>::value));
};

template<class T>
inline dense_matrix<T> matcl::convert_object(const dense_matrix<T>& A, ti::ti_object ti)
{
    return dense_matrix<T>(convert_object(A.to_matrix(), ti));
};

template<class T>
inline sparse_matrix<T> matcl::delrows(const sparse_matrix<T>& A, const colon& c)
{
    return A.delrows(c);
};

template<class T>
inline sparse_matrix<T> matcl::delrows(sparse_matrix<T>&& A, const colon& c)
{
    return std::move(A).delrows(c);
};

template<class T>
inline sparse_matrix<T> matcl::delcols(const sparse_matrix<T>& A, const colon& c)
{
    return A.delcols(c);
};

template<class T>
inline sparse_matrix<T> matcl::delcols(sparse_matrix<T>&& A, const colon& c)
{
    return std::move(A).delcols(c);
};

template<class T>
inline sparse_matrix<T> matcl::delrowscols(const sparse_matrix<T>& A, const colon& c1, const colon& c2)
{
    return A.delrowscols(c1,c2);
};

template<class T>
inline sparse_matrix<T> matcl::delrowscols(sparse_matrix<T>&& A, const colon& c1, const colon& c2)
{
    return std::move(A).delrowscols(c1,c2);
};

template<class T>
inline sparse_matrix<T> matcl::horzcat(const sparse_matrix<T>& A, const sparse_matrix<T>& B)
{
    return sparse_matrix<T>(horzcat(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::horzcat(const std::vector<sparse_matrix<T>>& mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::horzcat(std::initializer_list<sparse_matrix<T>> mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::horzcat(const sparse_matrix<T>& A, const band_matrix<T>& B)
{
    return sparse_matrix<T>(horzcat(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::horzcat(const band_matrix<T>& A, const sparse_matrix<T>& B)
{
    return sparse_matrix<T>(horzcat(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::horzcat(const band_matrix<T>& A, const band_matrix<T>& B)
{
    return sparse_matrix<T>(horzcat(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::horzcat(const std::vector<band_matrix<T>>& mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::horzcat(std::initializer_list<band_matrix<T>> mat_list)
{
    mat_row ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::vertcat(const sparse_matrix<T>& A, const sparse_matrix<T>& B)
{
    return sparse_matrix<T>(vertcat(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::vertcat(const std::vector<sparse_matrix<T>>& mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::vertcat(std::initializer_list<sparse_matrix<T>> mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::blkdiag(const sparse_matrix<T>& A, const sparse_matrix<T>& B)
{
    return sparse_matrix<T>(blkdiag(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::blkdiag(const std::vector<sparse_matrix<T>>& mat_list)
{
    std::vector<Matrix> vec_mat;
    for( auto it: mat_list)
        vec_mat.push_back(Matrix(it));

    Matrix ret = blkdiag(vec_mat);
    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::blkdiag(std::initializer_list<sparse_matrix<T>> mat_list)
{
    std::vector<Matrix> vec_mat;
    for( auto it: mat_list)
        vec_mat.push_back(Matrix(it));

    Matrix ret = blkdiag(vec_mat);
    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::vertcat(const band_matrix<T>& A, const sparse_matrix<T>& B)
{
    return sparse_matrix<T>(vertcat(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::blkdiag(const band_matrix<T>& A, const sparse_matrix<T>& B)
{
    return sparse_matrix<T>(blkdiag(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::vertcat(const sparse_matrix<T>& A, const band_matrix<T>& B)
{
    return sparse_matrix<T>(vertcat(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::blkdiag(const sparse_matrix<T>& A, const band_matrix<T>& B)
{
    return sparse_matrix<T>(blkdiag(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::vertcat(const band_matrix<T>& A, const band_matrix<T>& B)
{
    return sparse_matrix<T>(vertcat(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::vertcat(const std::vector<band_matrix<T>>& mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::vertcat(std::initializer_list<band_matrix<T>> mat_list)
{
    mat_col ret;
    for( auto it: mat_list)
        ret.add(Matrix(it));

    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::blkdiag(const band_matrix<T>& A, const band_matrix<T>& B)
{
    return sparse_matrix<T>(blkdiag(Matrix(A), Matrix(B)));
};

template<class T>
inline sparse_matrix<T> matcl::blkdiag(const std::vector<band_matrix<T>>& mat_list)
{
    std::vector<Matrix> vec_mat;
    for( auto it: mat_list)
        vec_mat.push_back(Matrix(it));

    Matrix ret = blkdiag(vec_mat);
    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::blkdiag(std::initializer_list<band_matrix<T>> mat_list)
{
    std::vector<Matrix> vec_mat;
    for( auto it: mat_list)
        vec_mat.push_back(Matrix(it));

    Matrix ret = blkdiag(vec_mat);
    return sparse_matrix<T>(ret);
};

template<class T>
inline sparse_matrix<T> matcl::repmat(const sparse_matrix<T>& A, Integer m, Integer n)
{
    return sparse_matrix<T>(repmat(Matrix(A), m, n));
};

template<class T>
inline sparse_matrix<T> matcl::sparse(const sparse_matrix<T>& m)
{
    return m;
};

template<class T>
inline sparse_matrix<T> matcl::sparse(const dense_matrix<T>& m)
{
    return sparse_matrix<T>(sparse(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::sparse(const band_matrix<T>& m)
{
    return sparse_matrix<T>(sparse(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::clone(const sparse_matrix<T>& m)
{
    return sparse_matrix<T>(clone(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::trans(const sparse_matrix<T>& m)
{
    return sparse_matrix<T>(trans(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::ctrans(const sparse_matrix<T>& m)
{
    return sparse_matrix<T>(ctrans(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::trans(const sparse_matrix<T>& m, trans_type t)
{
    return sparse_matrix<T>(trans(Matrix(m), t));
};

template<class T>
inline sparse_matrix<T> matcl::trans(const sparse_matrix<T>& m, trans_type_ext t)
{
    return sparse_matrix<T>(trans(Matrix(m), t));
};

template<class T>
inline sparse_matrix<T> matcl::vec(const sparse_matrix<T>& m)
{
    return sparse_matrix<T>(vec(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::tril(const sparse_matrix<T>& m, Integer d)
{
    return sparse_matrix<T>(tril(Matrix(m),d));
};

template<class T>
inline sparse_matrix<T> matcl::tril(sparse_matrix<T>&& m, Integer d)
{
    return sparse_matrix<T>(tril(Matrix(std::move(m)),d));
};

template<class T>
inline band_matrix<T> matcl::select_band(const sparse_matrix<T>& m, Integer ld, Integer ud)
{
    return band_matrix<T>(select_band(Matrix(m),ld, ud));
};

template<class T>
inline sparse_matrix<T> matcl::triu(const sparse_matrix<T>& m, Integer d)
{
    return sparse_matrix<T>(triu(Matrix(m),d));
};

template<class T>
inline sparse_matrix<T> matcl::triu(sparse_matrix<T>&& m, Integer d)
{
    return sparse_matrix<T>(triu(Matrix(std::move(m)),d));
};

template<class T>
inline sparse_matrix<T> matcl::flipud(const sparse_matrix<T>& m)
{
    return sparse_matrix<T>(flipud(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::fliplr(const sparse_matrix<T>& m)
{
    return sparse_matrix<T>(fliplr(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::reshape(const sparse_matrix<T>& A, Integer m, Integer n)
{
    return sparse_matrix<T>(reshape(Matrix(A),m,n));
};

template<class T>
inline dense_matrix<T> matcl::get_diag(const sparse_matrix<T>& m, Integer d)
{
    return dense_matrix<T>(get_diag(Matrix(m),d));
};

template<class T>
inline sparse_matrix<T> matcl::rot90(const sparse_matrix<T>& m, Integer n)
{
    return sparse_matrix<T>(rot90(Matrix(m),n));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const sparse_matrix<T>& v)
{
    return dense_matrix<Integer>(find(Matrix(v)));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const sparse_matrix<T>& v, const test_function& t)
{
    return dense_matrix<Integer>(find(Matrix(v), t));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const sparse_matrix<T>& v, const test_type_function<T>& t)
{
    return dense_matrix<Integer>(find(Matrix(v), t.to_test_function()));
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const sparse_matrix<T>& v)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v));

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const sparse_matrix<T>& v,const test_function& t)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v),t);

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const sparse_matrix<T>& v,const test_type_function<T>& t)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v),t.to_test_function());

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const sparse_matrix<T>& v)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v));

    dense_matrix<Integer>   ret1(R);
    dense_matrix<Integer>   ret2(C);
    dense_matrix<T>         ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const sparse_matrix<T>& v, const test_function& t)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v), t);

    dense_matrix<Integer>       ret1(R);
    dense_matrix<Integer>       ret2(C);
    dense_matrix<T>             ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const sparse_matrix<T>& v, const test_type_function<T>& t)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v), t.to_test_function());

    dense_matrix<Integer>       ret1(R);
    dense_matrix<Integer>       ret2(C);
    dense_matrix<T>             ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline sparse_matrix<T> matcl::sort(const sparse_matrix<T>& v, int dim, bool asceding)
{
    return sparse_matrix<T>(sort(Matrix(v), dim, asceding));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sort2(const sparse_matrix<T>& v, int dim, bool asceding)
{
    Matrix S, I;
    tie(S,I) = sort2(Matrix(v),dim,asceding);

    sparse_matrix<T>         ret1(S);
    dense_matrix<Integer>       ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline sparse_matrix<T> matcl::sortrows(const sparse_matrix<T>& v)
{
    return sparse_matrix<T>(sortrows(Matrix(v)));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sortrows2(const sparse_matrix<T>& v)
{
    Matrix S, I;
    tie(S,I) = sort2(Matrix(v));

    sparse_matrix<T>        ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline sparse_matrix<T> matcl::sortrows(const sparse_matrix<T>& v, const Matrix& cols)
{
    return sparse_matrix<T>(sortrows(Matrix(v), cols));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sortrows2(const sparse_matrix<T>& v, const Matrix& cols)
{
    Matrix S, I;
    tie(S,I) = sortrows2(Matrix(v),cols);

    sparse_matrix<T>        ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
}

template<class T>
inline sparse_matrix<T> matcl::sortcols(const sparse_matrix<T>& v)
{
    return sparse_matrix<T>(sortcols(Matrix(v)));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sortcols2(const sparse_matrix<T>& v)
{
    Matrix S, I;
    tie(S,I) = sortcols2(Matrix(v));

    sparse_matrix<T>         ret1(S);
    dense_matrix<Integer>    ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline sparse_matrix<T> matcl::sortcols(const sparse_matrix<T>& v, const Matrix& dims)
{
    return sparse_matrix<T>(sortcols(Matrix(v), dims));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sortcols2(const sparse_matrix<T>& v, const Matrix& dims)
{
    Matrix S, I;
    tie(S,I) = sortcols2(Matrix(v),dims);

    sparse_matrix<T>         ret1(S);
    dense_matrix<Integer>    ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline dense_matrix<Integer> matcl::issorted(const sparse_matrix<T>& v, int dim, bool asceding)
{
    return dense_matrix<Integer>(issorted(Matrix(v),dim,asceding));
};

template<class T>
inline sparse_matrix<T> matcl::drop_sparse(const sparse_matrix<T>& A, Real tol)
{
    return sparse_matrix<T>(drop_sparse(A.to_matrix(), tol));
};

// convert the matrix to a new type given by new_type
template<class New_mat, class T>
inline New_mat matcl::convert(const sparse_matrix<T>& A)
{
    return New_mat(A.to_matrix());
};

template<class T>
inline sparse_matrix<T> matcl::convert_object(const sparse_matrix<T>& A, ti::ti_object ti)
{
    return sparse_matrix<T>(convert_object(A.to_matrix(), ti));
};

template<class S, class T>
inline sparse_matrix<S> matcl::convert_value(const sparse_matrix<T>& m)
{
    return sparse_matrix<S>(convert_value(Matrix(m), matrix_traits::value_code<S>::value));
};

template<class T>
inline sparse_matrix<T> matcl::delrows(const band_matrix<T>& A, const colon& c)
{
    return A.delrows(c);
};

template<class T>
inline sparse_matrix<T> matcl::delrows(band_matrix<T>&& A, const colon& c)
{
    return std::move(A).delrows(c);
};

template<class T>
inline sparse_matrix<T> matcl::delcols(const band_matrix<T>& A, const colon& c)
{
    return A.delcols(c);
};

template<class T>
inline sparse_matrix<T> matcl::delcols(band_matrix<T>&& A, const colon& c)
{
    return std::move(A).delcols(c);
};

template<class T>
inline sparse_matrix<T> matcl::delrowscols(const band_matrix<T>& A, const colon& c1, const colon& c2)
{
    return A.delrowscols(c1,c2);
};

template<class T>
inline sparse_matrix<T> matcl::delrowscols(band_matrix<T>&& A, const colon& c1, const colon& c2)
{
    return std::move(A).delrowscols(c1,c2);
};

template<class T>
inline sparse_matrix<T> matcl::repmat(const band_matrix<T>& A, Integer m, Integer n)
{
    return sparse_matrix<T>(repmat(Matrix(A), m, n));
};

template<class T>
inline band_matrix<T> matcl::band(const sparse_matrix<T>& m)
{
    return band_matrix<T>(band(Matrix(m)));
};

template<class T>
inline band_matrix<T> matcl::band(const dense_matrix<T>& m)
{
    return band_matrix<T>(band(Matrix(m)));
};

template<class T>
inline band_matrix<T> matcl::band(const band_matrix<T>& m)
{
    return m;
};

template<class T>
inline band_matrix<T> matcl::clone(const band_matrix<T>& m)
{
    return band_matrix<T>(clone(Matrix(m)));
};

template<class T>
inline band_matrix<T> matcl::trans(const band_matrix<T>& m)
{
    return band_matrix<T>(trans(Matrix(m)));
};

template<class T>
inline band_matrix<T> matcl::ctrans(const band_matrix<T>& m)
{
    return band_matrix<T>(ctrans(Matrix(m)));
};

template<class T>
inline band_matrix<T> matcl::trans(const band_matrix<T>& m, trans_type t)
{
    return band_matrix<T>(trans(Matrix(m), t));
};

template<class T>
inline band_matrix<T> matcl::trans(const band_matrix<T>& m, trans_type_ext t)
{
    return band_matrix<T>(trans(Matrix(m), t));
};

template<class T>
inline sparse_matrix<T> matcl::vec(const band_matrix<T>& m)
{
    return sparse_matrix<T>(vec(Matrix(m)));
};

template<class T>
inline band_matrix<T> matcl::tril(const band_matrix<T>& m, Integer d)
{
    return band_matrix<T>(tril(Matrix(m),d));
};

template<class T>
inline band_matrix<T> matcl::tril(band_matrix<T>&& m, Integer d)
{
    return band_matrix<T>(tril(Matrix(std::move(m)),d));
};

template<class T>
inline band_matrix<T> matcl::select_band(const band_matrix<T>& m, Integer ld, Integer ud)
{
    return band_matrix<T>(select_band(Matrix(m),ld, ud));
};

template<class T>
inline band_matrix<T> matcl::triu(const band_matrix<T>& m, Integer d)
{
    return band_matrix<T>(triu(Matrix(m),d));
};

template<class T>
inline band_matrix<T> matcl::triu(band_matrix<T>&& m, Integer d)
{
    return band_matrix<T>(triu(Matrix(std::move(m)),d));
};

template<class T>
inline sparse_matrix<T> matcl::flipud(const band_matrix<T>& m)
{
    return sparse_matrix<T>(flipud(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::fliplr(const band_matrix<T>& m)
{
    return sparse_matrix<T>(fliplr(Matrix(m)));
};

template<class T>
inline sparse_matrix<T> matcl::reshape(const band_matrix<T>& A, Integer m, Integer n)
{
    return sparse_matrix<T>(reshape(Matrix(A),m,n));
};

template<class T>
inline dense_matrix<T> matcl::get_diag(const band_matrix<T>& m, Integer d)
{
    return dense_matrix<T>(get_diag(Matrix(m),d));
};

template<class T>
inline sparse_matrix<T> matcl::rot90(const band_matrix<T>& m, Integer n)
{
    return sparse_matrix<T>(rot90(Matrix(m),n));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const band_matrix<T>& v)
{
    return dense_matrix<Integer>(find(Matrix(v)));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const band_matrix<T>& v, const test_function& t)
{
    return dense_matrix<Integer>(find(Matrix(v), t));
};

template<class T>
inline dense_matrix<Integer> matcl::find(const band_matrix<T>& v, const test_type_function<T>& t)
{
    return dense_matrix<Integer>(find(Matrix(v), t.to_test_function()));
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const band_matrix<T>& v)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v));

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const band_matrix<T>& v,const test_function& t)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v),t);

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>>
matcl::find2(const band_matrix<T>& v,const test_type_function<T>& t)
{
    Matrix R, C;
    tie(R,C) = find2(Matrix(v),t.to_test_function());

    dense_matrix<Integer> ret1(R);
    dense_matrix<Integer> ret2(C);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const band_matrix<T>& v)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v));

    dense_matrix<Integer>   ret1(R);
    dense_matrix<Integer>   ret2(C);
    dense_matrix<T>         ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const band_matrix<T>& v, const test_function& t)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v), t);

    dense_matrix<Integer>       ret1(R);
    dense_matrix<Integer>       ret2(C);
    dense_matrix<T>             ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>
matcl::find3(const band_matrix<T>& v, const test_type_function<T>& t)
{
    Matrix R, C, X;
    tie(R,C, X) = find3(Matrix(v), t.to_test_function());

    dense_matrix<Integer>       ret1(R);
    dense_matrix<Integer>       ret2(C);
    dense_matrix<T>             ret3(X);

    return tuple<dense_matrix<Integer>, dense_matrix<Integer>, dense_matrix<T>>(ret1, ret2, ret3);
};

template<class T>
inline sparse_matrix<T> matcl::sort(const band_matrix<T>& v, int dim, bool asceding)
{
    return sparse_matrix<T>(sort(Matrix(v), dim, asceding));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sort2(const band_matrix<T>& v, int dim, bool asceding)
{
    Matrix S, I;
    tie(S,I) = sort2(Matrix(v),dim,asceding);

    sparse_matrix<T>         ret1(S);
    dense_matrix<Integer>       ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline sparse_matrix<T> matcl::sortrows(const band_matrix<T>& v)
{
    return sparse_matrix<T>(sortrows(Matrix(v)));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sortrows2(const band_matrix<T>& v)
{
    Matrix S, I;
    tie(S,I) = sort2(Matrix(v));

    sparse_matrix<T>        ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline sparse_matrix<T> matcl::sortrows(const band_matrix<T>& v, const Matrix& cols)
{
    return sparse_matrix<T>(sortrows(Matrix(v), cols));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sortrows2(const band_matrix<T>& v, const Matrix& cols)
{
    Matrix S, I;
    tie(S,I) = sortrows2(Matrix(v),cols);

    sparse_matrix<T>        ret1(S);
    dense_matrix<Integer>   ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
}

template<class T>
inline sparse_matrix<T> matcl::sortcols(const band_matrix<T>& v)
{
    return sparse_matrix<T>(sortcols(Matrix(v)));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sortcols2(const band_matrix<T>& v)
{
    Matrix S, I;
    tie(S,I) = sortcols2(Matrix(v));

    sparse_matrix<T>         ret1(S);
    dense_matrix<Integer>    ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline sparse_matrix<T> matcl::sortcols(const band_matrix<T>& v, const Matrix& dims)
{
    return sparse_matrix<T>(sortcols(Matrix(v), dims));
};

template<class T>
inline tuple<sparse_matrix<T>, dense_matrix<Integer>>
matcl::sortcols2(const band_matrix<T>& v, const Matrix& dims)
{
    Matrix S, I;
    tie(S,I) = sortcols2(Matrix(v),dims);

    sparse_matrix<T>         ret1(S);
    dense_matrix<Integer>    ret2(I);

    return tuple<sparse_matrix<T>, dense_matrix<Integer>>(ret1, ret2);
};

template<class T>
inline dense_matrix<Integer> matcl::issorted(const band_matrix<T>& v, int dim, bool asceding)
{
    return dense_matrix<Integer>(issorted(Matrix(v),dim,asceding));
};

template<class T>
inline band_matrix<T> matcl::drop_sparse(const band_matrix<T>& A, Real tol)
{
    (void)tol;
    return A;
};

// convert the matrix to a new type given by new_type
template<class New_mat, class T>
inline New_mat matcl::convert(const band_matrix<T>& A)
{
    return New_mat(A.to_matrix());
};

template<class S, class T>
inline band_matrix<S> matcl::convert_value(const band_matrix<T>& m)
{
    return band_matrix<S>(convert_value(Matrix(m), matrix_traits::value_code<S>::value));
};

template<class T>
inline band_matrix<T> matcl::convert_object(const band_matrix<T>& A, ti::ti_object ti)
{
    return band_matrix<T>(convert_object(A.to_matrix(), ti));
};

//--------------------------------------------------------------
//      functions operating on rows or columns
//--------------------------------------------------------------

template<class T>
inline dense_matrix<Integer> matcl::nnz(const dense_matrix<T>& A, Integer dim)
{
    return dense_matrix<Integer>(matcl::nnz(matcl::Matrix(A),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const dense_matrix<T>& v, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const dense_matrix<T>& v,const test_function& t, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),t,dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const dense_matrix<T>& v, const test_type_function<T>& t, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),t.to_test_function(),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const dense_matrix<T>& v, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const dense_matrix<T>& v,const test_function& t, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),t,dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const dense_matrix<T>& v, const test_type_function<T>& t, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),t.to_test_function(),dim));
};

template<class T>
inline dense_matrix<T> matcl::sum(const dense_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::sum(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::prod(const dense_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::prod(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::cumsum(const dense_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::cumsum(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::cumprod(const dense_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::cumprod(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<typename md::unify_types<T,Float>::type>
matcl::mean(const dense_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::unify_types<T,Float>::type>;
    return ret_type(matcl::mean(matcl::Matrix(v),dim));
};

template<class T>
dense_matrix<typename md::real_type_int_real<T>::type>
inline matcl::std(const dense_matrix<T>& v, int dim, bool unbiased)
{
    using ret_type = dense_matrix<typename md::real_type_int_real<T>::type>;
    return ret_type(matcl::std(matcl::Matrix(v),dim,unbiased));
};

template<class T>
inline dense_matrix<T> matcl::min_d(const dense_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::min_d(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::max_d(const dense_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::max_d(matcl::Matrix(v),dim));
};

template<class T>
dense_matrix<typename md::real_type<T>::type>
inline matcl::min_abs_d(const dense_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::real_type<T>::type>;
    return ret_type(matcl::min_abs_d(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<typename md::real_type<T>::type>
matcl::max_abs_d(const dense_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::real_type<T>::type>;
    return ret_type(matcl::max_abs_d(matcl::Matrix(v),dim));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::min2(const dense_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<T>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::min2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::max2(const dense_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<T>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::max2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
matcl::min_abs2(const dense_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<typename md::real_type<T>::type>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::min_abs2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
matcl::max_abs2(const dense_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<typename md::real_type<T>::type>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::max_abs2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline dense_matrix<Integer> matcl::nnz(const sparse_matrix<T>& A, Integer dim)
{
    return dense_matrix<Integer>(matcl::nnz(matcl::Matrix(A),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const sparse_matrix<T>& v, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const sparse_matrix<T>& v,const test_function& t, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),t,dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const sparse_matrix<T>& v, const test_type_function<T>& t, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),t.to_test_function(),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const sparse_matrix<T>& v, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const sparse_matrix<T>& v,const test_function& t, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),t,dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const sparse_matrix<T>& v, const test_type_function<T>& t, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),t.to_test_function(),dim));
};

template<class T>
inline dense_matrix<T> matcl::sum(const sparse_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::sum(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::prod(const sparse_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::prod(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::cumsum(const sparse_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::cumsum(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::cumprod(const sparse_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::cumprod(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<typename md::unify_types<T,Float>::type>
matcl::mean(const sparse_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::unify_types<T,Float>::type>;
    return ret_type(matcl::mean(matcl::Matrix(v),dim));
};

template<class T>
dense_matrix<typename md::real_type_int_real<T>::type>
inline matcl::std(const sparse_matrix<T>& v, int dim, bool unbiased)
{
    using ret_type = dense_matrix<typename md::real_type_int_real<T>::type>;
    return ret_type(matcl::std(matcl::Matrix(v),dim,unbiased));
};

template<class T>
inline dense_matrix<T> matcl::min_d(const sparse_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::min_d(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::max_d(const sparse_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::max_d(matcl::Matrix(v),dim));
};

template<class T>
dense_matrix<typename md::real_type<T>::type>
inline matcl::min_abs_d(const sparse_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::real_type<T>::type>;
    return ret_type(matcl::min_abs_d(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<typename md::real_type<T>::type>
matcl::max_abs_d(const sparse_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::real_type<T>::type>;
    return ret_type(matcl::max_abs_d(matcl::Matrix(v),dim));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::min2(const sparse_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<T>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::min2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::max2(const sparse_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<T>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::max2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
matcl::min_abs2(const sparse_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<typename md::real_type<T>::type>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::min_abs2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
matcl::max_abs2(const sparse_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<typename md::real_type<T>::type>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::max_abs2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline dense_matrix<Integer> matcl::nnz(const band_matrix<T>& A, Integer dim)
{
    return dense_matrix<Integer>(matcl::nnz(matcl::Matrix(A),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const band_matrix<T>& v, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const band_matrix<T>& v,const test_function& t, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),t,dim));
};

template<class T>
inline dense_matrix<Integer> matcl::all(const band_matrix<T>& v, const test_type_function<T>& t, int dim)
{
    return dense_matrix<Integer>(matcl::all(matcl::Matrix(v),t.to_test_function(),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const band_matrix<T>& v, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const band_matrix<T>& v,const test_function& t, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),t,dim));
};

template<class T>
inline dense_matrix<Integer> matcl::any(const band_matrix<T>& v, const test_type_function<T>& t, int dim)
{
    return dense_matrix<Integer>(matcl::any(matcl::Matrix(v),t.to_test_function(),dim));
};

template<class T>
inline dense_matrix<T> matcl::sum(const band_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::sum(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::prod(const band_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::prod(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::cumsum(const band_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::cumsum(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::cumprod(const band_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::cumprod(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<typename md::unify_types<T,Float>::type>
matcl::mean(const band_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::unify_types<T,Float>::type>;
    return ret_type(matcl::mean(matcl::Matrix(v),dim));
};

template<class T>
dense_matrix<typename md::real_type_int_real<T>::type>
inline matcl::std(const band_matrix<T>& v, int dim, bool unbiased)
{
    using ret_type = dense_matrix<typename md::real_type_int_real<T>::type>;
    return ret_type(matcl::std(matcl::Matrix(v),dim,unbiased));
};

template<class T>
inline dense_matrix<T> matcl::min_d(const band_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::min_d(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<T> matcl::max_d(const band_matrix<T>& v, int dim)
{
    return dense_matrix<T>(matcl::max_d(matcl::Matrix(v),dim));
};

template<class T>
dense_matrix<typename md::real_type<T>::type>
inline matcl::min_abs_d(const band_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::real_type<T>::type>;
    return ret_type(matcl::min_abs_d(matcl::Matrix(v),dim));
};

template<class T>
inline dense_matrix<typename md::real_type<T>::type>
matcl::max_abs_d(const band_matrix<T>& v, int dim)
{
    using ret_type = dense_matrix<typename md::real_type<T>::type>;
    return ret_type(matcl::max_abs_d(matcl::Matrix(v),dim));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::min2(const band_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<T>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::min2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<T>, dense_matrix<Integer>>
matcl::max2(const band_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<T>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::max2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
matcl::min_abs2(const band_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<typename md::real_type<T>::type>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::min_abs2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

template<class T>
inline tuple<dense_matrix<typename md::real_type<T>::type>, dense_matrix<Integer>>
matcl::max_abs2(const band_matrix<T>& v, int dim)
{
    using first_type    = dense_matrix<typename md::real_type<T>::type>;
    using second_type   = dense_matrix<Integer>;
    using ret_type      = tuple<first_type, second_type>;

    auto ret = matcl::max_abs2(matcl::Matrix(v),dim);
    return ret_type(first_type(ret.get<1>()), second_type(ret.get<2>()));
};

//--------------------------------------------------------------
//      functions operating on all elements
//--------------------------------------------------------------
template<class T>
inline Integer matcl::nnz_vec(const dense_matrix<T>& A)
{
    return matcl::nnz_vec(matcl::Matrix(A));
};

template<class T>
inline bool matcl::all_vec(const dense_matrix<T>& v)
{
    return matcl::all_vec(matcl::Matrix(v));
};

template<class T>
inline bool matcl::all_vec(const dense_matrix<T>& v,const test_function& t)
{
    return matcl::all_vec(matcl::Matrix(v),t);
};

template<class T>
inline bool matcl::all_vec(const dense_matrix<T>& v, const test_type_function<T>& t)
{
    return matcl::all_vec(matcl::Matrix(v),t.to_test_function());
};

template<class T>
inline bool matcl::any_vec(const dense_matrix<T>& v)
{
    return matcl::any_vec(matcl::Matrix(v));
};

template<class T>
inline bool matcl::any_vec(const dense_matrix<T>& v,const test_function& t)
{
    return matcl::any_vec(matcl::Matrix(v),t);
};

template<class T>
inline bool matcl::any_vec(const dense_matrix<T>& v, const test_type_function<T>& t)
{
    return matcl::any_vec(matcl::Matrix(v),t.to_test_function());
};

template<class T>
inline T matcl::sum_vec(const dense_matrix<T>& v)
{
    return matcl::sum_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline T matcl::prod_vec(const dense_matrix<T>& v)
{
    return matcl::prod_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline typename md::unify_types<T,Float>::type
matcl::mean_vec(const dense_matrix<T>& v)
{
    using ret_type = typename md::unify_types<T,Float>::type;
    return matcl::mean_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline typename md::real_type_int_real<T>::type
matcl::std_vec(const dense_matrix<T>& v, bool unbiased)
{
    using ret_type = typename md::real_type_int_real<T>::type;
    return matcl::std_vec(matcl::Matrix(v),unbiased).get_scalar<ret_type>();
};

template<class T>
inline T matcl::min_vec(const dense_matrix<T>& v)
{
    return matcl::min_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline T matcl::max_vec(const dense_matrix<T>& v)
{
    return matcl::max_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline typename md::real_type<T>::type
matcl::min_abs_vec(const dense_matrix<T>& v)
{
    using ret_type = typename md::real_type<T>::type;
    return matcl::min_abs_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline typename md::real_type<T>::type
matcl::max_abs_vec(const dense_matrix<T>& v)
{
    using ret_type = typename md::real_type<T>::type;
    return matcl::max_abs_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline tuple<T, Integer, Integer >
matcl::min2_vec(const dense_matrix<T>& v)
{
    using first_type    = T;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::min2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<T, Integer, Integer >
matcl::max2_vec(const dense_matrix<T>& v)
{
    using first_type    = T;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::max2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<typename md::real_type<T>::type, Integer, Integer >
matcl::min_abs2_vec(const dense_matrix<T>& v)
{
    using first_type    = typename md::real_type<T>::type;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::min_abs2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<typename md::real_type<T>::type, Integer, Integer >
matcl::max_abs2_vec(const dense_matrix<T>& v)
{
    using first_type    = typename md::real_type<T>::type;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::max_abs2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline Integer matcl::nnz_vec(const sparse_matrix<T>& A)
{
    return matcl::nnz_vec(matcl::Matrix(A));
};

template<class T>
inline bool matcl::all_vec(const sparse_matrix<T>& v)
{
    return matcl::all_vec(matcl::Matrix(v));
};

template<class T>
inline bool matcl::all_vec(const sparse_matrix<T>& v,const test_function& t)
{
    return matcl::all_vec(matcl::Matrix(v),t);
};

template<class T>
inline bool matcl::all_vec(const sparse_matrix<T>& v, const test_type_function<T>& t)
{
    return matcl::all_vec(matcl::Matrix(v),t.to_test_function());
};

template<class T>
inline bool matcl::any_vec(const sparse_matrix<T>& v)
{
    return matcl::any_vec(matcl::Matrix(v));
};

template<class T>
inline bool matcl::any_vec(const sparse_matrix<T>& v,const test_function& t)
{
    return matcl::any_vec(matcl::Matrix(v),t);
};

template<class T>
inline bool matcl::any_vec(const sparse_matrix<T>& v, const test_type_function<T>& t)
{
    return matcl::any_vec(matcl::Matrix(v),t.to_test_function());
};

template<class T>
inline T matcl::sum_vec(const sparse_matrix<T>& v)
{
    return matcl::sum_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline T matcl::prod_vec(const sparse_matrix<T>& v)
{
    return matcl::prod_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline typename md::unify_types<T,Float>::type
matcl::mean_vec(const sparse_matrix<T>& v)
{
    using ret_type = typename md::unify_types<T,Float>::type;
    return matcl::mean_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline typename md::real_type_int_real<T>::type
matcl::std_vec(const sparse_matrix<T>& v, bool unbiased)
{
    using ret_type = typename md::real_type_int_real<T>::type;
    return matcl::std_vec(matcl::Matrix(v),unbiased).get_scalar<ret_type>();
};

template<class T>
inline T matcl::min_vec(const sparse_matrix<T>& v)
{
    return matcl::min_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline T matcl::max_vec(const sparse_matrix<T>& v)
{
    return matcl::max_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline typename md::real_type<T>::type
matcl::min_abs_vec(const sparse_matrix<T>& v)
{
    using ret_type = typename md::real_type<T>::type;
    return matcl::min_abs_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline typename md::real_type<T>::type
matcl::max_abs_vec(const sparse_matrix<T>& v)
{
    using ret_type = typename md::real_type<T>::type;
    return matcl::max_abs_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline tuple<T, Integer, Integer >
matcl::min2_vec(const sparse_matrix<T>& v)
{
    using first_type    = T;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::min2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<T, Integer, Integer >
matcl::max2_vec(const sparse_matrix<T>& v)
{
    using first_type    = T;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::max2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<typename md::real_type<T>::type, Integer, Integer >
matcl::min_abs2_vec(const sparse_matrix<T>& v)
{
    using first_type    = typename md::real_type<T>::type;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::min_abs2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<typename md::real_type<T>::type, Integer, Integer >
matcl::max_abs2_vec(const sparse_matrix<T>& v)
{
    using first_type    = typename md::real_type<T>::type;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::max_abs2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline Integer matcl::nnz_vec(const band_matrix<T>& A)
{
    return matcl::nnz_vec(matcl::Matrix(A));
};

template<class T>
inline bool matcl::all_vec(const band_matrix<T>& v)
{
    return matcl::all_vec(matcl::Matrix(v));
};

template<class T>
inline bool matcl::all_vec(const band_matrix<T>& v,const test_function& t)
{
    return matcl::all_vec(matcl::Matrix(v),t);
};

template<class T>
inline bool matcl::all_vec(const band_matrix<T>& v, const test_type_function<T>& t)
{
    return matcl::all_vec(matcl::Matrix(v),t.to_test_function());
};

template<class T>
inline bool matcl::any_vec(const band_matrix<T>& v)
{
    return matcl::any_vec(matcl::Matrix(v));
};

template<class T>
inline bool matcl::any_vec(const band_matrix<T>& v,const test_function& t)
{
    return matcl::any_vec(matcl::Matrix(v),t);
};

template<class T>
inline bool matcl::any_vec(const band_matrix<T>& v, const test_type_function<T>& t)
{
    return matcl::any_vec(matcl::Matrix(v),t.to_test_function());
};

template<class T>
inline T matcl::sum_vec(const band_matrix<T>& v)
{
    return matcl::sum_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline T matcl::prod_vec(const band_matrix<T>& v)
{
    return matcl::prod_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline typename md::unify_types<T,Float>::type
matcl::mean_vec(const band_matrix<T>& v)
{
    using ret_type = typename md::unify_types<T,Float>::type;
    return matcl::mean_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline typename md::real_type_int_real<T>::type
matcl::std_vec(const band_matrix<T>& v, bool unbiased)
{
    using ret_type = typename md::real_type_int_real<T>::type;
    return matcl::std_vec(matcl::Matrix(v),unbiased).get_scalar<ret_type>();
};

template<class T>
inline T matcl::min_vec(const band_matrix<T>& v)
{
    return matcl::min_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline T matcl::max_vec(const band_matrix<T>& v)
{
    return matcl::max_vec(matcl::Matrix(v)).get_scalar<T>();
};

template<class T>
inline typename md::real_type<T>::type
matcl::min_abs_vec(const band_matrix<T>& v)
{
    using ret_type = typename md::real_type<T>::type;
    return matcl::min_abs_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline typename md::real_type<T>::type
matcl::max_abs_vec(const band_matrix<T>& v)
{
    using ret_type = typename md::real_type<T>::type;
    return matcl::max_abs_vec(matcl::Matrix(v)).get_scalar<ret_type>();
};

template<class T>
inline tuple<T, Integer, Integer >
matcl::min2_vec(const band_matrix<T>& v)
{
    using first_type    = T;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::min2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<T, Integer, Integer >
matcl::max2_vec(const band_matrix<T>& v)
{
    using first_type    = T;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::max2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<typename md::real_type<T>::type, Integer, Integer >
matcl::min_abs2_vec(const band_matrix<T>& v)
{
    using first_type    = typename md::real_type<T>::type;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::min_abs2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

template<class T>
inline tuple<typename md::real_type<T>::type, Integer, Integer >
matcl::max_abs2_vec(const band_matrix<T>& v)
{
    using first_type    = typename md::real_type<T>::type;
    using ret_type      = tuple<first_type, Integer, Integer >;

    auto ret = matcl::max_abs2_vec(matcl::Matrix(v));

    return ret_type(ret.get<1>().get_scalar<first_type>(),ret.get<2>(), ret.get<3>());
};

};
