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

#include "matcl-linalg/decompositions/eig/svd_range.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-core/utils/workspace.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/decompositions/hess.h"
#include "matcl-linalg/decompositions/eig/schur_utils.h"

#if 0

namespace matcl { namespace details
{

// commented out; gesvdx is generally broken and not faster than standard svd

template<class V>
struct singsel_str
{
    using VR            = typename md::real_type<V>::type;
    using Mat           = raw::Matrix<V,struct_dense>;
    using Mat_R         = raw::Matrix<VR,struct_dense>;

    static void eval(bool need_UV, Matrix& U, Matrix& S, Matrix& VV, const Matrix& A, 
                        bool range, Real VL, Real VU, Integer IF, Integer IL)    
    {
        Integer M           = A.rows();
        Integer N           = A.cols();
        Integer K           = std::min(M, N);

        Mat mat_A           = A.impl<Mat>().make_unique();
        const char* JOBU    = need_UV ? "V" : "N";
        const char* JOBVT   = need_UV ? "V" : "N";
        const char* RANGE   = range ? "V" : "I";

        Integer KU          = range ? K : IL - IF + 1;
        Integer MU          = M;
        Integer NU          = need_UV ? KU : 0;
        Integer MV          = need_UV ? K : 0;
        Integer NV          = N;

        V* ptr_A            = mat_A.ptr();
        Integer LDA         = mat_A.ld();
        Integer NS;

        //gesvdx write out of bound for ptr_S array
        Integer mult        = 2;

        Mat_R mat_S(ti::ti_empty(), K * mult, 1);
        Mat mat_U(ti::ti_empty(), MU, NU);
        Mat mat_V(ti::ti_empty(), MV, NV);

        VR* ptr_S           = mat_S.ptr();
        V* ptr_U            = mat_U.ptr();
        V* ptr_V            = mat_V.ptr();
        Integer LDU         = mat_U.ld();
        Integer LDV         = mat_V.ld();

        using iworkspace    = matcl::pod_workspace<Integer>;
        iworkspace IWORK    = iworkspace(12 * K * mult);
        Integer* IWORK_ptr  = IWORK.ptr();

        bool is_complex     = md::is_complex<V>::value;
        Integer lrwork      = is_complex ? K * (K*2+15*K) : 0;
        using rworkspace    = matcl::pod_workspace<VR>;
        rworkspace RWORK    = rworkspace(lrwork);
        VR* RWORK_ptr       = RWORK.ptr();

        V work_query;
        Integer info;

        lapack::gesvdx(JOBU, JOBVT, RANGE, M, N, lap(ptr_A), LDA, (VR)VL, (VR)VU, IF, IL, NS, 
                       lap(ptr_S), lap(ptr_U), LDU, lap(ptr_V), LDV, lap(&work_query), -1, RWORK_ptr, 
                       IWORK_ptr, info);

        Integer lwork       = (Integer)real(work_query) * mult;

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = (V*)WORK.ptr();

        lapack::gesvdx(JOBU, JOBVT, RANGE, M, N, lap(ptr_A), LDA, (VR)VL, (VR)VU, IF, IL, NS, 
                       lap(ptr_S), lap(ptr_U), LDU, lap(ptr_V), LDV, lap(ptr_WORK), lwork, RWORK_ptr,
                       IWORK_ptr, info);

        if (info < 0)
            throw error::error_general("invalid argument passed to gesvdx");

        else if (info > 0)
            throw error::error_svd();

        S                   = Matrix(mat_S, false)(colon(1, NS));

        if (need_UV == false)
            return;

        S                   = bdiag(S);
        U                   = Matrix(mat_U, true)(colon(), colon(1, NS));
        VV                  = Matrix(mat_V, true)(colon(1, NS), colon());
        VV                  = ctrans(VV);

        if (U.rows() == U.cols())
        {
            struct_flag str_unit;
            str_unit.set_user(unitary_flag());

            U.add_struct(str_unit);
        };
        if (VV.rows() == VV.cols())
        {
            struct_flag str_unit;
            str_unit.set_user(unitary_flag());

            VV.add_struct(str_unit);
        };
    };    
};

struct visit_singsel : public md::extract_type_switch<void, visit_singsel,true>
{
    using base_type = md::extract_type_switch<void, visit_singsel,true>;

    template<class Mat, class ... Arg>
    static void eval(const matcl::Matrix&, const Mat&, Arg&& ... args)
    {
        using V0    = typename Mat::value_type;
        using V     = typename md::unify_types<V0, Float>::type;
        return singsel_str<V>::eval(std::forward<Arg>(args)...);
    };

    template<class Mat, class ... Arg>
    static void eval_scalar(const matcl::Matrix&, const Mat&, Arg&& ... args)
    {
        using V = typename md::unify_types<Mat, Float>::type;
        return singsel_str<V>::eval(std::forward<Arg>(args)...);
    };

    template<class Str, class ... Arg>
    static void eval(const matcl::Matrix&, const raw::Matrix<Object, Str>&, Arg&& ...)
    {
        throw error::object_value_type_not_allowed("singsel");
    };

    template<class ... Arg>
    static void eval_scalar(const matcl::Matrix&, const Object&, Arg&& ...)
    {
        throw error::object_value_type_not_allowed("singsel");
    };
};

static void eval_impl(Matrix& U, Matrix& S, Matrix& V, bool need_UV, const Matrix& A, 
                bool range, Real VL, Real VU, Integer IF, Integer IL)
{
    Integer M           = A.rows();
    Integer N           = A.cols();
    Integer K           = std::min(M, N);
    value_code vc       = A.get_value_code();
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (K == 0)
    {
        vc              = matrix_traits::real_value_type(vc);

        if (need_UV == true)
        {
            U           = zeros(M,0, vc);
            V           = zeros(N,0, vc);
            S           = zeros(0,0,vc);
        }
        else
        {
            S           = zeros(0,1,vc);
        };
        return;
    }

    if (range == false)
    {
        IF              = std::max(IF, 1);
        IL              = std::max(IL, 1);
        IF              = std::min(IF, K);
        IL              = std::min(IL, K);
        IL              = std::max(IF, IL);
    };

    bool isv            = A.all_finite();

    if (range == true)
        isv             = isv && is_finite(VL) && is_finite(VU);

    if (isv == false)
    {
        Integer N_sel   = 1;
        if (range == false)
            N_sel       = IL - IF + 1;

        if (need_UV == true)
        {
            U           = make_nan_matrix(M, N_sel, vc);
            V           = make_nan_matrix(N, N_sel, vc);
            S           = make_nan_matrix(N_sel, N_sel,vc);
        }
        else
        {
            S           = make_nan_matrix(N_sel,1,vc);
        };

        return;
    };

    if (is_diag(A) == true)
        return eval_diag(need_UV, U, S, V, A, range, VL, VU, IF, IL);

    Matrix sel = zeros(0,0,vc);
    visit_singsel::make<const Matrix&>(sel, need_UV, U, S, V, A, range, VL, VU, IF, IL);
};

}};

namespace matcl
{

void svd_range::eval_range(Matrix& ret, const Matrix& A, Real VL, Real VU)
{
    Matrix U, V;
    return details::eval_impl(U, ret, V, false, A, true, VL, VU, -1, -1);
};

void svd_range::eval_index(Matrix& ret, const Matrix& A, Integer IF, Integer IL)
{
    Matrix U, V;
    return details::eval_impl(U, ret, V, false, A, false, -1, -1, IF, IL);
};

void svd_range::eval2_range(mat_tup_3& ret, const Matrix& A, Real VL, Real VU)
{
    Matrix U, S, V;
    details::eval_impl(U, S, V, true, A, true, VL, VU, -1, -1);
    ret = mat_tup_3(U, S, V);
}
void svd_range::eval2_index(mat_tup_3& ret, const Matrix& A, Integer IF, Integer IL)
{
    Matrix U, S, V;
    details::eval_impl(U, S, V, true, A, false, -1, -1, IF, IL);
    ret = mat_tup_3(U, S, V);
}

// compute selected singular values eigenvalues
//     E = singsel_range(A, VL, VU);
//     E = singsel_index(A, IF, IL);
// A       - an M x N matrix
// VL, VU  - find singular values in half-open interval (VL,VU], VL >= 0
// IF, IL  - find singular values with index IF through IL in set of
//           all singular values sorted decreasingly
// E       - selected singular values sorted decreasingly 
//
// not available for sparse and band matrices
MATCL_LINALG_EXPORT Matrix  singsel_range(const Matrix& A, Real VL, Real VU);
MATCL_LINALG_EXPORT Matrix  singsel_range(Matrix&& A, Real VL, Real VU);
MATCL_LINALG_EXPORT Matrix  singsel_index(const Matrix& A, Integer IF, Integer IL);
MATCL_LINALG_EXPORT Matrix  singsel_index(Matrix&& A, Integer IF, Integer IL);

// compute selected singular values and singular vectors
//     [U, S, V] = singsel_range2(A, VL, VU);
//     [U, S, V] = singsel_index2(A, IF, IL);
// A       - a M x N matrix
// VL, VU  - find singular values in half-open interval (VL,VU], VL >= 0
// IF, IL  - find singular values with index IF through IL in set of
//           all singular values sorted decreasingly
// U       - unitary matrix of size M x K, where K is number of selected
//           singular values
// V       - unitary matrix of size N x K
// S       - diagonal matrix of size K x K with singular values sorted
//           decreasingly
// not available for sparse and band matrices
MATCL_LINALG_EXPORT mat_tup_3  singsel_range2(const Matrix& A, Real VL, Real VU);
MATCL_LINALG_EXPORT mat_tup_3  singsel_range2(Matrix&& A, Real VL, Real VU);
MATCL_LINALG_EXPORT mat_tup_3  singsel_index2(const Matrix& A, Integer IF, Integer IL);
MATCL_LINALG_EXPORT mat_tup_3  singsel_index2(Matrix&& A, Integer IF, Integer IL);

}

#endif

