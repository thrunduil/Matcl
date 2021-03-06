/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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
#include "matcl-internals/func/raw_manip.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;

namespace details
{

template<class MP>
struct manip_reshape_helper
{
    using value_type        = typename MP::value_type;
    using struct_type       = typename MP::struct_type;
    using struct_type_ret   = typename matcl::details::select_if
                            <
                                std::is_same<struct_type,struct_banded>::value,
                                struct_sparse,
                                struct_type
                            >::type;

    using ret_type_vec      = Matrix<value_type,struct_type_ret>;
    using ret_type_flip     = Matrix<value_type,struct_type_ret>;
    using ret_type_reshape  = Matrix<value_type,struct_type_ret>;
    using ret_type_repmat   = Matrix<value_type,struct_type_ret>;

    static void     eval_vec(matcl::Matrix& ret, const MP& m);
    static void     eval_reshape(matcl::Matrix& ret, const MP& A, Integer m, Integer n);
    static void     eval_flipud(matcl::Matrix& ret, const MP& m);
    static void     eval_fliplr(matcl::Matrix& ret, const MP& m);        
    static void     eval_repmat(matcl::Matrix& ret, const MP& A, Integer m, Integer n);
};

};

template<class val_type>
Matrix<val_type,struct_dense> full(const Matrix<val_type,struct_sparse>& A);

template<class val_type>
Matrix<val_type,struct_sparse> sparse(const Matrix<val_type,struct_dense>& A);

template<class val_type>
Matrix<val_type,struct_sparse> sparse(const Matrix<val_type,struct_banded>& A);

template<class val_type>
Matrix<val_type,struct_banded> band(const Matrix<val_type,struct_dense>& A);

template<class val_type>
Matrix<val_type,struct_banded> band(const Matrix<val_type,struct_sparse>& A);

template<class V, class S>
bool is_tril(const Matrix<V,S>& A)
{
    return (get_ud(A,0) == 0);
};

template<class V, class S>
bool is_triu(const Matrix<V,S>& A)
{
    return (get_ld(A,0) == 0);
};

template<class V>
Integer nnz_total(const Matrix<V,struct_dense>& A);

template<class V>
Integer nnz_total(const Matrix<V,struct_sparse>& A);

template<class V>
Integer nnz_total(const Matrix<V,struct_banded>& A);

// tol is used to check equality of a and b, i.e.
// a ~== b if |a-b| < eps(a) * tol for tol > 0 and
// a ~== b if a == b for tol <= 0
template<class V>
bool is_sym(const Matrix<V,struct_dense>& A, Real tol, bool use_flag);

template<class V>
bool is_sym(const Matrix<V,struct_sparse>& A, Real tol, bool use_flag);

template<class V>
bool is_sym(const Matrix<V,struct_banded>& A, Real tol, bool use_flag);

template<class V>
bool is_her(const Matrix<V,struct_dense>& A, Real tol, bool use_flag);

template<class V>
bool is_her(const Matrix<V,struct_sparse>& A, Real tol, bool use_flag);

template<class V>
bool is_her(const Matrix<V,struct_banded>& A, Real tol, bool use_flag);

template<class Mat>
struct all_finite_helper
{
    static bool eval(const Mat& mat);
};

};};
