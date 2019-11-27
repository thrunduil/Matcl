/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-linalg/decompositions/eig/schur_range.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-core/utils/workspace.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/decompositions/hess.h"
#include "matcl-linalg/decompositions/eig/schur_utils.h"

namespace matcl { namespace details
{

template<class V, bool Is_comples = md::is_complex<V>::value>
struct eigsel_str
{
    using VR            = typename md::real_type<V>::type;

    static void eval(Matrix& ret_E, const Matrix& diag, const Matrix& subdiag, 
                        bool range, Real VL, Real VU, Integer IF, Integer IL)
    {
        Integer N           = diag.length();
        const V* ptr_D      = diag.get_array<V>();
        const V* ptr_E      = subdiag.get_array<V>();

        using Mat_D         = raw::Matrix<VR,struct_dense>;

        Mat_D W(ti::ti_empty(), N, 1);
        VR* ptr_W           = W.ptr();

        using workspace     = matcl::pod_workspace<VR>;
        workspace WORK      = workspace(5*N);
        VR* ptr_WORK        = WORK.ptr();

        using iworkspace    = matcl::pod_workspace<Integer>;
        iworkspace IWORK    = iworkspace(2*N + 3*N);
        Integer* IWORK_ptr  = IWORK.ptr();

        Integer n_eig;
        const char* RANGE   = range ? "V" : "I";
        const char* ORDER   = "E";
        
        VR ABSTOL           = 0;

        Integer n_split;        
        
        Integer* IBLOCK     = IWORK_ptr;
        Integer* ISPLIT     = IWORK_ptr + N;
        Integer* ptr_IWORK  = IWORK_ptr + 2*N;

        Integer info;

        lapack::stebz<V>(RANGE, ORDER, N, (V)VL, (V)VU, IF, IL, ABSTOL, lap(ptr_D), lap(ptr_E), n_eig, 
              n_split, lap(ptr_W), IBLOCK, ISPLIT, ptr_WORK, ptr_IWORK, info);

        if (info < 0)
            throw error::error_general("invalid argument passed to stebz");
        else if (info > 0)
            throw error::error_schur();

        ret_E               = Matrix(W.resize(n_eig, 1), false);
    };
};

template<class V>
struct eigsel_str<V,true>
{
    using VR            = typename md::real_type<V>::type;

    static void eval(Matrix& ret_E, const Matrix& diag, const Matrix& subdiag, 
                        bool range, Real VL, Real VU, Integer IF, Integer IL)
    {
        Integer N           = diag.length();
        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK_U    = workspace(N);

        V* ptr_U            = reinterpret_cast<V*>(WORK_U.ptr());

        Matrix D1r          = make_tridiag_subdiag_real(subdiag, ptr_U, true);

        return eigsel_str<VR>::eval(ret_E, diag, D1r, range, VL, VU, IF, IL);
    };
};

struct visit_eigsel : public md::extract_type_switch<void, visit_eigsel,true>
{
    using base_type = md::extract_type_switch<void, visit_eigsel,true>;

    template<class Mat, class ... Arg>
    static void eval(const matcl::Matrix&, const Mat&, Arg&& ... args)
    {
        using V0    = typename Mat::value_type;
        using V     = typename md::unify_types<V0, Float>::type;
        return eigsel_str<V>::eval(std::forward<Arg>(args)...);
    };

    template<class Mat, class ... Arg>
    static void eval_scalar(const matcl::Matrix&, const Mat&, Arg&& ... args)
    {
        using V = typename md::unify_types<Mat, Float>::type;
        return eigsel_str<V>::eval(std::forward<Arg>(args)...);
    };

    template<class Str, class ... Arg>
    static void eval(const matcl::Matrix&, const raw::Matrix<Object, Str>&, Arg&& ...)
    {
        throw error::object_value_type_not_allowed("eigsel_tridiag");
    };

    template<class ... Arg>
    static void eval_scalar(const matcl::Matrix&, const Object&, Arg&& ...)
    {
        throw error::object_value_type_not_allowed("eigsel_tridiag");
    };
};

static void eval_tridiag_impl(Matrix& E, Matrix& V, bool need_V, const Matrix& diag, const Matrix& subdiag, 
                bool range, Real VL, Real VU, Integer IF, Integer IL)
{
    check_tridiag(diag, subdiag);

    Integer N           = diag.length();
    value_code vc_A     = diag.get_value_code();
    value_code vc_B     = subdiag.get_value_code();
    value_code vc       = matrix_traits::unify_value_types(vc_A, vc_B);
    vc                  = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (N == 0)
    {
        vc              = matrix_traits::real_value_type(vc);

        Matrix S        = zeros(0,0,vc);
        E               = S;
        V               = S;
        return;
    }

    if (N == 1)
    {
        Matrix I        = eye(1,1,vc);
        E               = diag;
        V               = I;
        return;
    }

    if (range == false)
    {
        IF              = std::max(IF, 1);
        IL              = std::max(IL, 1);
        IF              = std::min(IF, N);
        IL              = std::min(IL, N);
        IL              = std::max(IF, IL);
    };

    bool isv            = diag.all_finite() && subdiag.all_finite();

    if (range == true)
        isv             = isv && is_finite(VL) && is_finite(VU);

    if (isv == false)
    {
        Integer N_sel   = 1;
        if (range == false)
            N_sel       = IL - IF + 1;

        E               = make_nan_matrix(N_sel, 1, vc);
        if (need_V)
            V           = make_nan_matrix(N, N_sel, vc);

        return;
    };

    Matrix sel = zeros(0,0,vc);
    visit_eigsel::make<const Matrix&>(sel, E, diag, subdiag, range, VL, VU, IF, IL);

    if (need_V == false)
        return;

    Matrix H            = make_band_noinit(N,N,-1,1,vc);
    H.diag(0)           = diag;
    H.diag(-1)          = subdiag(colon(1,N-1));
    H.diag(1)           = conj(subdiag(colon(1,N-1)));    

    H.add_struct(predefined_struct_type::her);

    V = hess_right_eig(H, E);
};

}};

namespace matcl
{

void schur_range::eval_tridiag_range(Matrix& ret, const Matrix& diag, const Matrix& subdiag, Real VL, Real VU)
{
    Matrix V;
    return details::eval_tridiag_impl(ret, V, false, diag, subdiag, true, VL, VU, -1, -1);
};

void schur_range::eval_tridiag_index(Matrix& ret, const Matrix& diag, const Matrix& subdiag, Integer IF, Integer IL)
{
    Matrix V;
    return details::eval_tridiag_impl(ret, V, false, diag, subdiag, false, -1, -1, IF, IL);
};

void schur_range::eval_tridiag2_range(mat_tup_2& ret, const Matrix& diag, const Matrix& subdiag, Real VL, Real VU)
{
    Matrix E, V;
    details::eval_tridiag_impl(E, V, true, diag, subdiag, true, VL, VU, -1, -1);
    ret = mat_tup_2(E,V);
}
void schur_range::eval_tridiag2_index(mat_tup_2& ret, const Matrix& diag, const Matrix& subdiag, Integer IF, Integer IL)
{
    Matrix E, V;
    details::eval_tridiag_impl(E, V, true, diag, subdiag, false, -1, -1, IF, IL);
    ret = mat_tup_2(E,V);
}

void schur_range::eval_range(Matrix& ret, const Matrix& A, Real VL, Real VU)
{
    A.add_struct(predefined_struct_type::her);
    
    Matrix H        = hess(A);
    Matrix diag     = H.diag(0);
    Matrix subdiag  = H.diag(-1);
    return eval_tridiag_range(ret, diag, subdiag, VL, VU);
};
void schur_range::eval_index(Matrix& ret, const Matrix& A, Integer IF, Integer IL)
{
    A.add_struct(predefined_struct_type::her);

    Matrix H        = hess(A);
    Matrix diag     = H.diag(0);
    Matrix subdiag  = H.diag(-1);
    return eval_tridiag_index(ret, diag, subdiag, IF, IL);
};

void schur_range::eval2_range(mat_tup_2& ret, const Matrix& A, Real VL, Real VU)
{
    A.add_struct(predefined_struct_type::her);

    unitary_matrix U;
    Matrix H;
    tie(U,H) = hess2(A);
    
    Matrix E;
    Matrix diag     = H.diag(0);
    Matrix subdiag  = H.diag(-1);

    eval_tridiag_range(E, diag, subdiag, VL, VU);

    Matrix V        = hess_right_eig(H, E);
    V               = U * std::move(V);

    ret             = mat_tup_2(E,V);
    return;
};

void schur_range::eval2_index(mat_tup_2& ret, const Matrix& A, Integer IF, Integer IL)
{
    A.add_struct(predefined_struct_type::her);

    unitary_matrix U;
    Matrix H;
    tie(U,H) = hess2(A);

    Matrix E;
    Matrix diag     = H.diag(0);
    Matrix subdiag  = H.diag(-1);

    eval_tridiag_index(E, diag, subdiag, IF, IL);

    Matrix V        = hess_right_eig(H, E);
    V               = U * std::move(V);

    ret             = mat_tup_2(E,V);
    return;
};

}

