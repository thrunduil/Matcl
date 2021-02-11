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
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/details/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"

namespace matcl { namespace raw { namespace details
{

template<class Val_ret, class M1, class M2,class str_1, class str_2>
struct eval_mult { };

template<class Val_ret, class M1, class M2, class str_1>
struct eval_mult_abs { };

}}};

namespace matcl { namespace raw { namespace details
{

namespace md = matcl::details;

// if refcount is 1 then inplace version can be used; increase refcount of
// nontemporary matrices first
template<class M1,class M2>
struct mult_helper
{
    using val_type_1    = typename M1::value_type;
    using val_type_2    = typename M2::value_type;
    using str_type_1    = typename M1::struct_type;
    using str_type_2    = typename M2::struct_type;

    using Val_ret       = typename md::unify_types<val_type_1,val_type_2>::type;

    static void eval(matcl::Matrix& ret, const M1& A, const M2& B, trans_type t_A, trans_type t_B);
    static void eval_abs(matcl::Matrix& ret, const M1& A, const M2& X, trans_type t_A, const matcl::Matrix& C);
};

template<class T1, class T2>
struct mult_helper_mat_scal
{
    static void eval(matcl::Matrix& ret, const T1& mat, const T2& scal, trans_type t_A, trans_type t_B);
    static void eval_abs(matcl::Matrix& ret, const T1& mat, const T2& scal, trans_type t_A, const matcl::Matrix& C);
};

template<class T1, class T2>
struct mult_helper_scal_mat
{
    static void eval(matcl::Matrix& ret, const T1& scal, const T2& mat, trans_type t_A, trans_type t_B);
    static void eval_abs(matcl::Matrix& ret, const T1& scal, const T2& mat, const matcl::Matrix& C);
};

struct mult_helper_scal_scal
{
    static void eval_abs(matcl::Matrix& ret, const matcl::Matrix& A, const matcl::Matrix& B,
                         const matcl::Matrix& C);
};

template<class M1,class M2>
struct gemm_helper
{
    static void eval(const M1& A, const M2& B, trans_type t_A, trans_type t_B,
                     const matcl::Matrix& alpha, const matcl::Matrix& beta, matcl::Matrix& C,
                     Integer fr, Integer rows, Integer fc, Integer cols);
};

template<class T1, class T2>
struct gemm_helper_mat_scal
{
    static void eval(const T1& mat, const T2& scal, trans_type t_A, trans_type t_B,
                     const matcl::Matrix& alpha, const matcl::Matrix& beta, matcl::Matrix& C,
                     Integer fr, Integer rows, Integer fc, Integer cols);
};

template<class T1, class T2>
struct gemm_helper_scal_mat
{
    static void eval(const T1& scal, const T2& mat, trans_type t_A, trans_type t_B,
                     const matcl::Matrix& alpha, const matcl::Matrix& beta, matcl::Matrix& C,
                     Integer fr, Integer rows, Integer fc, Integer cols);
};

// C = a + beta*C
struct gemm_helper_scal_scal
{
    static void eval(const matcl::Matrix& a, const matcl::Matrix& beta, matcl::Matrix& C, 
                     Integer fr, Integer rows, Integer fc, Integer cols);
};

}}}
