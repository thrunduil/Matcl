/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include "matcl-linalg/decompositions/lu.h"
#include "matcl-internals/func/converter.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-matrep/matrix/permvec.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-linalg/lusol/lusol.h"
#include "matcl-internals/func/lapack_utils.h"

namespace matcl {  namespace details
{

template<class V> struct lu_dense   
{ 
    static void eval(lu_return_type& ret, const raw::Matrix<V,struct_dense>& mat, 
                     const matcl::options& opts); 
};
template<class V> struct lu_sparse  
{ 
    using VL    = typename lusol_value_type<V>::type;

    static void eval(lu_return_type& ret, const raw::Matrix<V,struct_sparse>& mat, 
                     const matcl::options& opts); 

    static void eval_lusol(lu_return_type& ret, const raw::Matrix<V,struct_sparse>& mat, 
                     const matcl::options& opts); 
    static void eval_superlu(lu_return_type& ret, const raw::Matrix<V,struct_sparse>& mat, 
                     const matcl::options& opts); 

    static void extract_factors(lusol::LUSOLrec<VL>*LUSOL, Matrix& out_L, Matrix& out_U,
                    const permvec& p, const permvec& q);
};

template<class V> struct lu_band    
{ 
    static void eval(lu_return_type& ret, const raw::Matrix<V,struct_banded>& mat, 
                     const matcl::options& opts); 
};

template<class V, class S>
struct lu_impl{};

template<class V>
struct lu_impl<V,struct_sparse>
{
    using matrix_type   = raw::Matrix<V,struct_sparse>;

    static void eval(lu_return_type& ret, const matrix_type& mat, const matcl::options& opts)
    {
        return lu_sparse<V>::eval(ret, mat,opts);
    };
};
template<class V>
struct lu_impl<V,struct_dense>
{
    using matrix_type   = raw::Matrix<V,struct_dense>;

    static void eval(lu_return_type& ret, const matrix_type& mat, const matcl::options& opts)
    {
        return lu_dense<V>::eval(ret, mat,opts);
    };
};
template<class V>
struct lu_impl<V,struct_banded>
{
    using matrix_type   = raw::Matrix<V,struct_banded>;

    static void eval(lu_return_type& ret, const matrix_type& mat, const matcl::options& opts)
    {
        return lu_band<V>::eval(ret, mat,opts);
    };
};

template<class V, class S>
struct lu
{
    using matrix_type   = raw::Matrix<V,S>;

    static void eval(lu_return_type& ret, const matrix_type& mat, const matcl::options& opts)
    {
        return details::lu_impl<V,S>::eval(ret, mat,opts);
    };
};
template<class S>
struct lu<Object,S>
{
    using matrix_type   = raw::Matrix<Object,S>;

    static void eval(lu_return_type&, const matrix_type&, const matcl::options&)
    {
        throw error::object_value_type_not_allowed("lu");
    };
};

template<class S>
struct lu<Integer,S>
{
    using matrix_type       = raw::Matrix<Integer,S>;
    using real_matrix_type  = raw::Matrix<Real,S>;

    static void eval(lu_return_type& ret, const matrix_type& mat, const matcl::options& opts)
    {
        real_matrix_type mat_real = raw::converter<real_matrix_type,matrix_type>::eval(mat);
        return lu<Real,S>::eval(ret, mat_real,opts);
    };
};

};};
