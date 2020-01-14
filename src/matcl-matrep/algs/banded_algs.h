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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/base/colon_info.h"

namespace matcl { namespace algorithm
{

    namespace details
    {
        template<class value_type>
        struct del_rows_banded_functor
        {
            using SM = raw::Matrix<value_type,struct_sparse>;
            using BM = raw::Matrix<value_type,struct_banded>;

            static void eval(Matrix& ret, const BM& mat,const matcl::details::colon_info& ci, bool rvalue);
        };

        template<class value_type>
        struct del_cols_banded_functor
        {
            using SM = raw::Matrix<value_type,struct_sparse>;
            using BM = raw::Matrix<value_type,struct_banded>;
            
            static void eval(Matrix& ret, const BM& mat,const matcl::details::colon_info& ci, bool rvalue);
        };

        template<class value_type>
        struct del_rowscols_banded_functor
        {
            using SM = raw::Matrix<value_type,struct_sparse>;
            using BM = raw::Matrix<value_type,struct_banded>;
            
            static void eval(Matrix& ret, const BM& mat,const matcl::details::colon_info& ci, bool rvalue);
            static void eval_00(Matrix& ret, const BM& mat,const matcl::details::colon_info& ci);
            static void eval_01(Matrix& ret, const BM& mat,const matcl::details::colon_info& ci);
            static void eval_10(Matrix& ret, const BM& mat,const matcl::details::colon_info& ci);
            static void eval_11(Matrix& ret, const BM& mat,const matcl::details::colon_info& ci);
        };

        template<class value_type>
        struct band_change_diag_functor
        {
            using DM = raw::Matrix<value_type,struct_dense>;
            using BM = raw::Matrix<value_type,struct_banded>;

            static void eval(Matrix& ret, const BM& mat, Integer d, const DM& val);
            static void eval(Matrix& ret, const BM& mat, Integer d, const value_type& val);
        };
    };

//eval A(c_i,:) = []
template<class value_type>
void del_rows_banded(Matrix& ret, const raw::Matrix<value_type,struct_banded>& mat,
                  const matcl::details::colon_info& ci, bool rvalue)
{
    return details::del_rows_banded_functor<value_type>::eval(ret,mat,ci,rvalue);
};

//eval A(:,c_i) = []
template<class value_type>
void del_cols_banded(Matrix& ret, const raw::Matrix<value_type,struct_banded>& mat,
                  const matcl::details::colon_info& ci, bool rvalue)
{
    return details::del_cols_banded_functor<value_type>::eval(ret,mat,ci,rvalue);
};

//eval A(c_i,c_i) = []
template<class value_type>
void del_rowscols_banded(Matrix& ret, const raw::Matrix<value_type,struct_banded>& mat,
                  const matcl::details::colon_info& ci, bool rvalue)
{
    return details::del_rowscols_banded_functor<value_type>::eval(ret,mat,ci,rvalue);
};

//eval A(d) = val for given diagonal
template<class value_type>
void
band_change_diag(Matrix& ret, const raw::Matrix<value_type,struct_banded>& mat, 
                 Integer d, const value_type& val)
{
    return details::band_change_diag_functor<value_type>::eval(ret, mat, d, val);
};

//eval A(d) = mat for given diagonal
template<class value_type>
void
band_change_diag(Matrix& ret, const raw::Matrix<value_type,struct_banded>& mat, Integer d, 
                       const raw::Matrix<value_type,struct_dense>& val)
{
    return details::band_change_diag_functor<value_type>::eval(ret, mat, d, val);
};

};};