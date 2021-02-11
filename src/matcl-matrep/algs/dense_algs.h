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

#include "matcl-matrep/base/colon_info.h"

namespace matcl { namespace algorithm
{
    namespace details
    {
        template<class DM>
        struct dense_change_entries_functor
        {
            using value_type = typename DM::value_type;
            static void eval(DM& mat,const matcl::details::colon_info& ci, const value_type& val);
        };

        template<class DM>
        struct dense_change_entries_functor_2
        {
            using value_type = typename DM::value_type;
            static void eval(DM& mat,const matcl::details::colon_info& ci, const value_type& val);
        };

        template<class DM>
        struct del_rows_dense_functor
        {
            static void eval(Matrix& ret, const DM& mat,const matcl::details::colon_info& ci, bool rvalue);
        };

        template<class DM>
        struct del_cols_dense_functor
        {
            static void eval(Matrix& ret, const DM& mat,const matcl::details::colon_info& ci, bool rvalue);
        };

        template<class DM>
        struct del_rowscols_dense_functor
        {
            static void eval(Matrix& ret, const DM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_00(Matrix& ret, const DM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_01(Matrix& ret, const DM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_10(Matrix& ret, const DM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_11(Matrix& ret, const DM& mat,const matcl::details::colon_info& ci, bool rvalue);
        };

        template<class DM>
        struct dense_change_diag_functor
        {
            using value_type = typename DM::value_type;
            static void eval(DM& mat, Integer d, const DM& val);
            static void eval(DM& mat, Integer d, const value_type& val);
        };
    };


//eval A(c) = val for given colon c
template<class value_type>
void dense_change_entries(raw::Matrix<value_type,struct_dense>& mat,
                    const matcl::details::colon_info& ci, const value_type& val)
{
    using DenseMatrix = raw::Matrix<value_type,struct_dense>;
    return details::dense_change_entries_functor<DenseMatrix>::eval(mat,ci,val);
};

//eval A(c1,c2) = val
template<class value_type>
void dense_change_entries_2(raw::Matrix<value_type,struct_dense>& mat,
                    const matcl::details::colon_info& ci, const value_type& val)
{
    using DenseMatrix = raw::Matrix<value_type,struct_dense>;
    return details::dense_change_entries_functor_2<DenseMatrix>::eval(mat,ci,val);
};

//eval A(c_i,:) = []
template<class value_type>
void del_rows_dense(Matrix& ret, const raw::Matrix<value_type,struct_dense>& mat,
                        const matcl::details::colon_info& ci, bool rvalue)
{
    using DenseMatrix = raw::Matrix<value_type,struct_dense>;
    return details::del_rows_dense_functor<DenseMatrix>::eval(ret,mat,ci,rvalue);
};

//eval A(:,c_i) = []
template<class value_type>
void del_cols_dense(Matrix& ret, const raw::Matrix<value_type,struct_dense>& mat,
                  const matcl::details::colon_info& ci, bool rvalue)
{
    using DenseMatrix = raw::Matrix<value_type,struct_dense>;
    return details::del_cols_dense_functor<DenseMatrix>::eval(ret,mat,ci,rvalue);
};

//eval A(c_i,c_i) = []
template<class value_type>
void del_rowscols_dense(Matrix& ret, const raw::Matrix<value_type,struct_dense>& mat,
                  const matcl::details::colon_info& ci, bool rvalue)
{
    using DenseMatrix = raw::Matrix<value_type,struct_dense>;
    return details::del_rowscols_dense_functor<DenseMatrix>::eval(ret,mat,ci,rvalue);
};

//eval A(d) = val for given diagonal
template<class value_type>
void dense_change_diag(raw::Matrix<value_type,struct_dense>& mat, Integer d, const value_type& val)
{
    using DenseMatrix = raw::Matrix<value_type,struct_dense>;
    return details::dense_change_diag_functor<DenseMatrix>::eval(mat,d,val);
};

//eval A(d) = mat for given diagonal
template<class value_type>
void dense_change_diag(raw::Matrix<value_type,struct_dense>& mat, Integer d, 
                       const raw::Matrix<value_type,struct_dense>& val)
{
    using DenseMatrix = raw::Matrix<value_type,struct_dense>;
    return details::dense_change_diag_functor<DenseMatrix>::eval(mat,d,val);
};

};};