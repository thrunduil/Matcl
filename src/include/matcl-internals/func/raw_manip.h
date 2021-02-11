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

#include "matcl-matrep/details/mpl.h"
#include "matcl-internals/base/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;

namespace details
{

template<class MP>
struct MATCL_MATREP_EXPORT manip_tr_helper
{
    using ret_type_tril = MP;                
    using ret_type_triu = MP;

    static void eval_tril(matcl::Matrix& ret, const MP& m, Integer d, bool rvalue);
    static void eval_triu(matcl::Matrix& ret, const MP& m, Integer d, bool rvalue);
    static void eval_select_band(matcl::Matrix& ret, const MP& m, Integer fd, Integer ld);
};

template<class struct_type, class Val_in, class Val_ret>
struct MATCL_MATREP_EXPORT manip_trans_converter_helper
{
    using in_type   = raw::Matrix<Val_in, struct_type>;
    using ret_type  = raw::Matrix<Val_ret, struct_type>;

    static void     eval_trans(ret_type& ret, const in_type& m);
    static void     eval_ctrans(ret_type& ret, const in_type& m);
};

template<class MP>
struct MATCL_MATREP_EXPORT manip_trans_helper
{
    using ret_type = MP;

    //ret is empty matrix
    static void     eval_trans(MP& ret, const MP& m);
    static void     eval_ctrans(MP& ret, const MP& m);
};

template<class Val, class Val_ret>
struct MATCL_MATREP_EXPORT manip_trans_reshaper_helper
{
    using in_type   = raw::Matrix<Val, struct_sparse>;
    using ret_type  = raw::Matrix<Val_ret, struct_sparse>;

    // drop rows with index higher than max_row, and columns with index higher than
    // max_col, create matrix of size ret_rows x ret_cols
    static void     eval_trans(ret_type& ret, const in_type& m, Integer max_row, Integer max_col, 
                               Integer ret_rows, Integer ret_cols);
    static void     eval_ctrans(ret_type& ret, const in_type& m, Integer max_row, Integer max_col, 
                                Integer ret_rows, Integer ret_cols);
};

template<class Val>
struct MATCL_MATREP_EXPORT manip_trans_helper_mat
{
    using MP        = raw::Matrix<Val,struct_dense>;

    //store trans in ret matrix; ret must be properly initialized
    static void     eval_trans(Val* ret_ptr, Integer ld, const MP& m);
    static void     eval_ctrans(Val* ret_ptr, Integer ld, const MP& m);
};

};

template<class V>
Integer MATCL_MATREP_EXPORT
get_ld(const Matrix<V,struct_dense>& A, Integer min, bool use_flag = true);

template<class V>
Integer MATCL_MATREP_EXPORT
get_ld(const Matrix<V,struct_sparse>& A, Integer min, bool use_flag = true);

template<class V>
Integer MATCL_MATREP_EXPORT
get_ld(const Matrix<V,struct_banded>& A, Integer min, bool use_flag = true);

template<class V>
Integer MATCL_MATREP_EXPORT
get_ud(const Matrix<V,struct_dense>& A, Integer min, bool use_flag = true);

template<class V>
Integer MATCL_MATREP_EXPORT
get_ud(const Matrix<V,struct_sparse>& A, Integer min, bool use_flag = true);

template<class V>
Integer MATCL_MATREP_EXPORT
get_ud(const Matrix<V,struct_banded>& A, Integer min, bool use_flag = true);

template<class V, class S>
bool is_diag(const Matrix<V,S>& A)
{
    return (get_ld(A,0) == 0) && (get_ud(A,0) == 0);
};

};};
