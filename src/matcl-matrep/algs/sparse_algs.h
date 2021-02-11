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

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/base/colon_info.h"

namespace matcl { namespace algorithm
{
    namespace details
    {
        template<class SM>
        struct drop_entries_functor
        {
            static SM eval(SM& mat,const matcl::details::colon_info&, Real tol);
        };

        template<class SM>
        struct drop_entries_functor_2
        {
            static SM eval(SM& mat,const matcl::details::colon_info&, Real tol);
        };

        template<class SM>
        struct zero_entries_functor
        {
            static SM eval(SM& mat,const matcl::details::colon_info&);
        };

        template<class SM>
        struct zero_entries_functor_2
        {
            static SM eval(SM& mat,const matcl::details::colon_info&);
        };

        template<class SM>
        struct change_entries_functor
        {
            using value_type = typename SM::value_type;
            static SM eval(SM& mat,const matcl::details::colon_info& ci, const value_type& val);
        };

        template<class SM>
        struct change_entries_functor_2
        {
            using value_type = typename SM::value_type;
            static SM eval(SM& mat,const matcl::details::colon_info& ci, const value_type& val);
        };

        template<class SM>
        struct add_entries_functor
        {
            using value_type = typename SM::value_type;
            static SM eval(SM& mat,const matcl::details::colon_info& ci);
        };

        template<class SM>
        struct add_entries_functor_2
        {
            using value_type = typename SM::value_type;
            static SM eval(SM& mat,const matcl::details::colon_info& ci);
        };

        template<class SM>
        struct get_submatrix_functor
        {
            public:
                using Vector = raw::integer_dense;

                static SM eval(const SM& mat,const matcl::details::colon_info& ci);

            private:
                static SM eval_00(const SM& mat,const Vector& ri,const Vector& ci);
                static SM eval_01(const SM& mat,const matcl::details::colon_info& ci);
                static SM eval_10(const SM& mat,const matcl::details::colon_info& ci);
                static SM eval_11(const SM& mat,const matcl::details::colon_info& ci);

                static SM get_cols_0(const SM& mat,const matcl::details::colon_info& ci);
                static SM get_cols_1(const SM& mat,const matcl::details::colon_info& ci);
        };

        template<class SM>
        struct get_submatrix_functor_2
        {
            public:
                static void eval(matcl::Matrix& ret, const SM& mat,const matcl::details::colon_info& ci);

            private:
                using Vector = raw::integer_dense;

                static SM eval_0(const SM& mat,const Vector& ri);
                static SM eval_0_dc(const SM& mat,const Vector& ri, const Vector& ci);
                static SM eval_1(const SM& mat,const matcl::details::colon_info& ci);

                static SM eval_0_increasing(const SM& mat,const Vector& ri, Integer n_rep);
                static SM eval_0_dc_increasing(const SM& mat,const Vector& r,const Vector& ci, Integer n_rep);
                static SM eval_0_decreasing(const SM& mat,const Vector& ri, Integer n_rep);
                static SM eval_1_increasing(const SM& mat,const matcl::details::colon_info& ci);
                static SM eval_1_decreasing(const SM& mat,const matcl::details::colon_info& ci);
        };

        template<class SM1,class SM2>
        struct change_submatrix_functor
        {
            public:
                static SM1 eval(const SM1& mat, const matcl::details::colon_info& ci,const SM2& mat2);
        };

        template<class SM1,class SM2>
        struct change_submatrix_functor_2
        {
            public:
                static SM1 eval(const SM1& mat,const matcl::details::colon_info& ci,const SM2& mat2);
        };

        template<class SM1,class DM2>
        struct change_submatrix_dense_functor
        {
            public:
                static SM1 eval(const SM1& mat, const matcl::details::colon_info& ci,const DM2& mat2);
        };

        template<class SM1,class DM2>
        struct change_submatrix_dense_functor_2
        {
            public:
                static SM1 eval(const SM1& mat,const matcl::details::colon_info& ci,const DM2& mat2);
        };

        template<class SM1,class DM2>
        struct change_submatrix_band_functor
        {
            public:
                static SM1 eval(const SM1& mat, const matcl::details::colon_info& ci,const DM2& mat2);
        };

        template<class SM1,class DM2>
        struct change_submatrix_band_functor_2
        {
            public:
                static SM1 eval(const SM1& mat,const matcl::details::colon_info& ci,const DM2& mat2);
        };

        template<class SM>
        struct del_rows_sparse_functor
        {
            static void eval(Matrix& ret, const SM& mat,const matcl::details::colon_info& ci, bool rvalue);
        };

        template<class SM>
        struct del_cols_sparse_functor
        {
            static void eval(Matrix& ret, const SM& mat,const matcl::details::colon_info& ci, bool rvalue);
        };

        template<class SM>
        struct del_rowscols_sparse_functor
        {
            static void eval(Matrix& ret, const SM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_00(Matrix& ret, const SM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_01(Matrix& ret, const SM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_10(Matrix& ret, const SM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_11(Matrix& ret, const SM& mat,const matcl::details::colon_info& ci, bool rvalue);
        };

        template<class value_type>
        struct sparse_change_diag_functor
        {
            using DM = raw::Matrix<value_type,struct_dense>;
            using SM = raw::Matrix<value_type,struct_sparse>;
            static SM eval(SM& mat, Integer d, const DM& val);
            static SM eval(SM& mat, Integer d, const value_type& val);
        };

        template<class value_type>
        struct sparse_drop_diag_functor
        {
            using SM = raw::Matrix<value_type,struct_sparse>;
            static SM eval(SM& mat, Integer d, Real tol);
            static SM eval_zero(SM& mat, Integer d);
        };

        template<class value_type>
        struct sparse_add_diag_functor
        {
            using SM = raw::Matrix<value_type,struct_sparse>;
            static SM eval(SM& mat, Integer d);
        };
    };

//removes entries A(r,c) given by vectors of rows and columns
template<class value_type>
raw::Matrix<value_type,struct_sparse>
drop_entries(raw::Matrix<value_type,struct_sparse>& mat,
                  const matcl::details::colon_info& ci, Real tol)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::drop_entries_functor<SparseMatrix>::eval(mat,ci, tol);
};

template<class value_type>
raw::Matrix<value_type,struct_sparse>
zero_entries(raw::Matrix<value_type,struct_sparse>& mat, const matcl::details::colon_info& ci)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::zero_entries_functor<SparseMatrix>::eval(mat,ci);
};

//removes entries A(r_i,c_i), i=1..s, where r, c are vectors of row and column indices with number of elements s
template<class value_type>
raw::Matrix<value_type,struct_sparse> 
drop_entries_2(raw::Matrix<value_type,struct_sparse>& mat, const matcl::details::colon_info& ci,
               Real tol)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::drop_entries_functor_2<SparseMatrix>::eval(mat,ci,tol);
};

template<class value_type>
raw::Matrix<value_type,struct_sparse> 
zero_entries_2(raw::Matrix<value_type,struct_sparse>& mat, const matcl::details::colon_info& ci)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::zero_entries_functor_2<SparseMatrix>::eval(mat,ci);
};

//eval A(r,c) = val for given vectors of rows and columns
template<class value_type>
raw::Matrix<value_type,struct_sparse>
change_entries(raw::Matrix<value_type,struct_sparse>& mat,
                  const matcl::details::colon_info& ci, const value_type& val)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::change_entries_functor<SparseMatrix>::eval(mat,ci,val);
};

//eval A(r_i,c_i) = val, i=1..s, where r, c are vectors of row and column indices with number of elements s
template<class value_type>
raw::Matrix<value_type,struct_sparse>
change_entries_2(raw::Matrix<value_type,struct_sparse>& mat,
                  const matcl::details::colon_info& ci, const value_type& val)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::change_entries_functor_2<SparseMatrix>::eval(mat,ci,val);
};

//add structure A(r,c) for given vectors of rows and columns
template<class value_type>
raw::Matrix<value_type,struct_sparse>
add_entries(raw::Matrix<value_type,struct_sparse>& mat, const matcl::details::colon_info& ci)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::add_entries_functor<SparseMatrix>::eval(mat,ci);
};

//add structure A(r_i,c_i), i=1..s, where r, c are vectors of row and column indices with number of elements s
template<class value_type>
raw::Matrix<value_type,struct_sparse>
add_entries_2(raw::Matrix<value_type,struct_sparse>& mat, const matcl::details::colon_info& ci)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::add_entries_functor_2<SparseMatrix>::eval(mat,ci);
};

//eval A(r,c) for given vectors of row and column indices
template<class value_type>
raw::Matrix<value_type,struct_sparse>
get_submatrix(const raw::Matrix<value_type,struct_sparse>& mat, const matcl::details::colon_info& ci)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::get_submatrix_functor<SparseMatrix>::eval(mat,ci);
};

//eval B(i,1) = A(r_i,c_i), i=1..s, where r, c are vectors of row and column indices with number of elements s
template<class value_type>
void get_submatrix_2(matcl::Matrix& ret, const raw::Matrix<value_type,struct_sparse>& mat, 
                const matcl::details::colon_info& ci)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    details::get_submatrix_functor_2<SparseMatrix>::eval(ret, mat,ci);
};

//eval A(r,c) = B, where r, c are vectors of row and column indices
template<class value_type_1,class value_type_2>
raw::Matrix<value_type_1,struct_sparse>
change_submatrix(const raw::Matrix<value_type_1,struct_sparse>& mat,
                   const matcl::details::colon_info& ci,
                   const raw::Matrix<value_type_2,struct_sparse>& mat_2)
{
    using SparseMatrix_1 = raw::Matrix<value_type_1,struct_sparse>;
    using SparseMatrix_2 = raw::Matrix<value_type_2,struct_sparse>;
    return details::change_submatrix_functor<SparseMatrix_1,SparseMatrix_2>::eval(mat,ci,mat_2);
};

//eval A(r_i,c_i) = B(i,1), i=1..s, where r, c are vectors of row and column indices with number of elements s
template<class value_type_1,class value_type_2>
raw::Matrix<value_type_1,struct_sparse>
change_submatrix_2(const raw::Matrix<value_type_1,struct_sparse>& mat,
                  const matcl::details::colon_info& ci,
                  const raw::Matrix<value_type_2,struct_sparse>& mat_2)
{
    using SparseMatrix_1 = raw::Matrix<value_type_1,struct_sparse>;
    using SparseMatrix_2 = raw::Matrix<value_type_2,struct_sparse>;
    return details::change_submatrix_functor_2<SparseMatrix_1,SparseMatrix_2>::eval(mat,ci,mat_2);
};

//eval A(r,c) = B, where r, c are vectors of row and column indices, B is a dense matrix
template<class value_type_1,class value_type_2>
raw::Matrix<value_type_1,struct_sparse>
change_submatrix_dense(const raw::Matrix<value_type_1,struct_sparse>& mat,
                   const matcl::details::colon_info& ci,
                   const raw::Matrix<value_type_2,struct_dense>& mat_2)
{
    using SparseMatrix_1    = raw::Matrix<value_type_1,struct_sparse>;
    using DenseMatrix_2     = raw::Matrix<value_type_2,struct_dense>;
    return details::change_submatrix_dense_functor<SparseMatrix_1,DenseMatrix_2>::eval(mat,ci,mat_2);
};

//eval A(r_i,c_i) = B(i,1), i=1..s, where r, c are vectors of row and column indices with number of elements s
//B is a dense matrix
template<class value_type_1,class value_type_2>
raw::Matrix<value_type_1,struct_sparse> 
change_submatrix_dense_2(const raw::Matrix<value_type_1,struct_sparse>& mat,
                  const matcl::details::colon_info& ci,
                  const raw::Matrix<value_type_2,struct_dense>& mat_2)
{
    using SparseMatrix_1    = raw::Matrix<value_type_1,struct_sparse>;
    using DenseMatrix_2     = raw::Matrix<value_type_2,struct_dense>;
    return details::change_submatrix_dense_functor_2<SparseMatrix_1,DenseMatrix_2>::eval(mat,ci,mat_2);
};

//eval A(r,c) = B, where r, c are vectors of row and column indices, B is a banded matrix
template<class value_type_1,class value_type_2>
raw::Matrix<value_type_1,struct_sparse>
change_submatrix_band(const raw::Matrix<value_type_1,struct_sparse>& mat,
                   const matcl::details::colon_info& ci,
                   const raw::Matrix<value_type_2,struct_banded>& mat_2)
{
    using SparseMatrix_1 = raw::Matrix<value_type_1,struct_sparse>;
    using DenseMatrix_2  = raw::Matrix<value_type_2,struct_banded>;
    return details::change_submatrix_band_functor<SparseMatrix_1,DenseMatrix_2>::eval(mat,ci,mat_2);
};

//eval A(r_i,c_i) = B(i,1), i=1..s, where r, c are vectors of row and column indices with number of elements s
//B is a band matrix
template<class value_type_1,class value_type_2>
raw::Matrix<value_type_1,struct_sparse>
change_submatrix_band_2(const raw::Matrix<value_type_1,struct_sparse>& mat,
                  const matcl::details::colon_info& ci,
                  const raw::Matrix<value_type_2,struct_banded>& mat_2)
{
    using SparseMatrix_1    = raw::Matrix<value_type_1,struct_sparse>;
    using DenseMatrix_2     = raw::Matrix<value_type_2,struct_banded>;
    return details::change_submatrix_band_functor_2<SparseMatrix_1,DenseMatrix_2>::eval(mat,ci,mat_2);
};

//eval A(c_i,:) = []
template<class value_type>
void del_rows_sparse(Matrix& ret, const raw::Matrix<value_type,struct_sparse>& mat,
                  const matcl::details::colon_info& ci, bool rvalue)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::del_rows_sparse_functor<SparseMatrix>::eval(ret,mat,ci,rvalue);
};

//eval A(:,c_i) = []
template<class value_type>
void del_cols_sparse(Matrix& ret, const raw::Matrix<value_type,struct_sparse>& mat,
                  const matcl::details::colon_info& ci, bool rvalue)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::del_cols_sparse_functor<SparseMatrix>::eval(ret,mat,ci,rvalue);
};

//eval A(c_i,c_i) = []
template<class value_type>
void del_rowscols_sparse(Matrix& ret, const raw::Matrix<value_type,struct_sparse>& mat,
                  const matcl::details::colon_info& ci, bool rvalue)
{
    using SparseMatrix = raw::Matrix<value_type,struct_sparse>;
    return details::del_rowscols_sparse_functor<SparseMatrix>::eval(ret,mat,ci,rvalue);
};

//eval A.diag(d) = val for given diagonal
template<class value_type>
raw::Matrix<value_type,struct_sparse>
sparse_change_diag(raw::Matrix<value_type,struct_sparse>& mat, Integer d, const value_type& val)
{
    return details::sparse_change_diag_functor<value_type>::eval(mat,d,val);
};

//eval drop(A.diag(d)) for given diagonal
template<class value_type>
raw::Matrix<value_type,struct_sparse>
sparse_drop_diag(raw::Matrix<value_type,struct_sparse>& mat, Integer d, Real tol)
{
    return details::sparse_drop_diag_functor<value_type>::eval(mat,d,tol);
};

//A.diag(d) = 0
template<class value_type>
raw::Matrix<value_type,struct_sparse>
zero_entries_diag(raw::Matrix<value_type,struct_sparse>& mat, Integer d)
{
    return details::sparse_drop_diag_functor<value_type>::eval_zero(mat,d);
};

//eval A.diag(d) = mat for given diagonal
template<class value_type>
raw::Matrix<value_type,struct_sparse>
sparse_change_diag(raw::Matrix<value_type,struct_sparse>& mat, Integer d, 
                       const raw::Matrix<value_type,struct_dense>& val)
{
    return details::sparse_change_diag_functor<value_type>::eval(mat,d,val);
};

//add diagonal d to structure
template<class value_type>
raw::Matrix<value_type,struct_sparse>
sparse_add_diag(raw::Matrix<value_type,struct_sparse>& mat, Integer d)
{
    return details::sparse_add_diag_functor<value_type>::eval(mat,d);
};

};};