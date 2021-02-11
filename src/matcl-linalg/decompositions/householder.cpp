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

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/general/config_linalg.h"
#include "matcl-linalg/decompositions/householder.h"
#include "matcl-linalg/general/linalg_exception.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/inplace.h"
#include "matcl-linalg/decompositions/householder_q.h"

namespace matcl { namespace details
{

template<class V, class S>
struct construct_str
{};

template<class V>
struct construct_str<V,struct_dense>
{
    using Mat = raw::Matrix<V,struct_dense>;
    static void eval(ret_construct_householder& ret, const Mat& mat0)
    {
        using VR    = typename md::real_type<V>::type;

        Integer M   = mat0.rows();
        Integer N   = mat0.cols();
        Integer K   = std::max(M,N);

        Mat A       = mat0.make_unique();
        A.get_struct().reset();

        V* ptr_A        = A.ptr();
        Integer A_ld    = A.ld();
        Integer inc_A   = (N == 1)? 1 : A_ld;
        V tau;

        using VL = typename details::lapack_value_type<V>::type;

        lapack::larfg<VL>(K, lap(ptr_A), lap(ptr_A + inc_A), inc_A, lap(&tau));

        V beta          = ptr_A[0];
        ptr_A[0]        = VR(1.0);

        ret = ret_construct_householder(Matrix(A,false), tau, beta);
    };

    static void eval_house_to_unit(unitary_matrix& ret, const Mat& T, const Mat& tau, Integer cols)
    {
        using unitary_matrix_data_ptr = unitary_matrix::unitary_matrix_data_ptr;

        Integer off = 0;

        bool isv    = tau.all_finite() && T.all_finite();
        if (isv == false)
        {
            ret = unitary_matrix::from_nan(T.rows(), cols, matrix_traits::value_code<V>::value);
            return;
        }
        else
        {
            unitary_matrix_data_ptr data(new householder_q<V>(cols,T,tau,T.ld(),off) );
            ret = unitary_matrix(data);
            return;
        };
    };
    static void eval_house_WY(Matrix& ret, const Mat& mat_V, const Mat& tau)
    {
        using VR    = typename md::real_type<V>::type;

        Integer N   = mat_V.rows();
        Integer K   = mat_V.cols();
        
        Mat T(mat_V.get_type(), K, K);

        const V* ptr_V  = mat_V.ptr();
        Integer V_ld    = mat_V.ld();
        V* ptr_T        = T.ptr();
        Integer T_ld    = T.ld();

        const V* ptr_tau= tau.ptr();

        using VL = typename details::lapack_value_type<V>::type;        

        lapack::larft<VL>("F", "C", N, K, lap(ptr_V), V_ld, lap(ptr_tau), lap(ptr_T), T_ld );

        Matrix ret_T;
        matcl::raw::inplace::make_triu<Mat>::eval(ret_T, T, 0);

        ret = ret_T;
        return;
    };
};

template<class V>
struct construct_str<V,struct_banded>
{
    using Mat   = raw::Matrix<V,struct_banded>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval(ret_construct_householder& ret, const Mat& mat)
    {
        using DM = raw::Matrix<V,struct_dense>;
        DM Ac = raw::converter<DM,Mat>::eval(mat);
        return construct_str<V,struct_dense>::eval(ret,Ac);
    };

    static void eval_house_to_unit(unitary_matrix& ret, const Mat& T, const Mat_D& tau, Integer cols)
    {
        using unitary_matrix_data_ptr = unitary_matrix::unitary_matrix_data_ptr;

        bool isv = T.all_finite() && tau.all_finite();

        if (isv == false)
        {
            ret = unitary_matrix::from_nan(T.rows(), cols, matrix_traits::value_code<V>::value);
            return;
        }
        else
        {
            unitary_matrix_data_ptr data(new householder_band_q<V>(cols,T,tau) );
            ret = unitary_matrix(data);
            return;
        };
    };
    static void eval_house_WY(Matrix& ret, const Mat& T, const Mat_D& tau)
    {
        using DM    = raw::Matrix<V,struct_dense>;

        if (T.has_diag(0) == false)
            throw error::band_matrix_with_main_diag_required(T.first_diag(), T.last_diag());

        Integer ldiags  = T.number_subdiagonals() + 1;        
        DM Ac           = raw::converter<DM,Mat>::eval(T.make_view(1, ldiags, T.cols()));
        return construct_str<V,struct_dense>::eval_house_WY(ret,Ac,tau);
    };
};

template<class V>
struct construct_str<V,struct_sparse>
{
    using Mat   = raw::Matrix<V,struct_sparse>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval(ret_construct_householder& ret, const Mat& mat)
    {
        using DM = raw::Matrix<V,struct_dense>;
        DM Ac = raw::converter<DM,Mat>::eval(mat);
        return construct_str<V,struct_dense>::eval(ret,Ac);
    };

    static void eval_house_to_unit(unitary_matrix& ret, const Mat& T, const Mat_D& tau, Integer cols)
    {
        using DM = raw::Matrix<V,struct_dense>;
        DM Ac = raw::converter<DM,Mat>::eval(T);
        return construct_str<V,struct_dense>::eval_house_to_unit(ret,Ac,tau, cols);
    };
    static void eval_house_WY(Matrix& ret, const Mat& T, const Mat_D& tau)
    {
        using DM = raw::Matrix<V,struct_dense>;
        DM Ac = raw::converter<DM,Mat>::eval(T);
        return construct_str<V,struct_dense>::eval_house_WY(ret,Ac,tau);
    };
};

template<class V, class S1, class S2>
struct construct_str2
{
    using Mat1  = raw::Matrix<V,S1>;
    using Mat2  = raw::Matrix<V,S2>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval_house_to_unit(unitary_matrix& ret, const Mat1& T, const Mat2& tau, Integer cols)
    {
        const Mat_D& tau_c  = raw::converter<Mat_D,Mat2>::eval(tau);
        return construct_str<V,S1>::eval_house_to_unit(ret,T,tau_c, cols);
    };
    static void eval_house_WY(Matrix& ret, const Mat1& T, const Mat2& tau)
    {
        const Mat_D& tau_c  = raw::converter<Mat_D,Mat2>::eval(tau);
        return construct_str<V,S1>::eval_house_WY(ret,T,tau_c);
    };
};


template<class V, class S>
struct construct_val
{
    using Mat = raw::Matrix<V,S>;
    static void eval(ret_construct_householder& ret, const Mat& mat)
    {
        return construct_str<V,S>::eval(ret,mat);
    };
};
template<class V, class S1, class S2>
struct construct_val2
{
    using Mat1 = raw::Matrix<V,S1>;
    using Mat2 = raw::Matrix<V,S2>;

    static void eval_house_to_unit(unitary_matrix& ret, const Mat1& A, const Mat2& B, Integer cols)
    {
        return construct_str2<V,S1,S2>::eval_house_to_unit(ret,A,B,cols);
    };
    static void eval_house_WY(Matrix& ret, const Mat1& A, const Mat2& B)
    {
        return construct_str2<V,S1,S2>::eval_house_WY(ret,A,B);
    };
};

template<class S>
struct construct_val<Integer,S>
{
    using Mat = raw::Matrix<Integer,S>;
    static void eval(ret_construct_householder& ret, const Mat& mat)
    {
        using DM = raw::Matrix<Real,struct_dense>;
        DM Ac = raw::converter<DM,Mat>::eval(mat);
        return construct_val<Real,struct_dense>::eval(ret,Ac);
    };
};

template<class S1, class S2>
struct construct_val2<Integer,S1, S2>
{
    using Mat1 = raw::Matrix<Integer,S1>;
    using Mat2 = raw::Matrix<Integer,S2>;

    static void eval_house_to_unit(unitary_matrix& ret, const Mat1& A, const Mat2& B, Integer cols)
    {
        using DM = raw::Matrix<Real,S1>;
        DM Ac   = raw::converter<DM,Mat1>::eval(A);
        DM Bc   = raw::converter<DM,Mat2>::eval(B);
        return construct_val2<Real,S1,S1>::eval_house_to_unit(ret,Ac,Bc,cols);
    };
    static void eval_house_WY(Matrix& ret, const Mat1& A, const Mat2& B)
    {
        using DM = raw::Matrix<Real,S1>;
        DM Ac   = raw::converter<DM,Mat1>::eval(A);
        DM Bc   = raw::converter<DM,Mat2>::eval(B);
        return construct_val2<Real,S1,S1>::eval_house_WY(ret,Ac,Bc);
    };
};

template<class S>
struct construct_val<Object,S>
{
    using Mat = raw::Matrix<Object,S>;

    static void eval(ret_construct_householder&, const Mat&)
    {
        throw error::object_value_type_not_allowed("construct_householder");
    };
};
template<class S1, class S2>
struct construct_val2<Object,S1, S2>
{
    using Mat1 = raw::Matrix<Object,S1>;
    using Mat2 = raw::Matrix<Object,S2>;

    static void eval_house_to_unit(unitary_matrix&, const Mat1&, const Mat2&, Integer)
    {
        throw error::object_value_type_not_allowed("householder_to_unitary");
    };
    static void eval_house_WY(Matrix&, const Mat1&, const Mat2&)
    {
        throw error::object_value_type_not_allowed("householder_WY");
    };
};
struct construct_householder_vis : public extract_type_switch<void, construct_householder_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, ret_construct_householder& ret)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return construct_val<V,S>::eval(ret, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& s, ret_construct_householder& ret)
    {
        using VT    = typename md::unify_types<T,Integer>::type;
        ret = ret_construct_householder(VT(s), VT(0.0), VT(s));
    };

    static void eval_scalar(const Matrix&, const Object&, ret_construct_householder&)
    {
        throw error::object_value_type_not_allowed("construct_householder");
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, ret_construct_householder&)
    {
        throw error::object_value_type_not_allowed("construct_householder");
    };
};

struct householder_to_unitary_vis : public extract_type2_switch<void,householder_to_unitary_vis, 
                                            mr::val_type_corrector_diag>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, unitary_matrix& ret, Integer cols)
    {
        using V     = typename T1::value_type;
        using S1    = typename T1::struct_type;
        using S2    = typename T2::struct_type;
        return construct_val2<V,S1,S2>::eval_house_to_unit(ret, A, B, cols);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, unitary_matrix& ret, Integer cols)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<Dense_1,Dense_2>(Ac,Bc, ret, cols);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& A, const T2& B, unitary_matrix& ret, Integer cols)
    {
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<T1,Dense_2>(A,Bc, ret,cols);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& A, const T2& B, unitary_matrix& ret, Integer cols)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<Dense_1,T2>(Ac,B, ret,cols);
    };
};

struct householder_WY_vis : public extract_type2_switch<void,householder_WY_vis, 
                                            mr::val_type_corrector_diag>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, Matrix& ret)
    {
        using V     = typename T1::value_type;
        using S1    = typename T1::struct_type;
        using S2    = typename T2::struct_type;
        return construct_val2<V,S1,S2>::eval_house_WY(ret, A, B);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, Matrix& ret)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<Dense_1,Dense_2>(Ac,Bc, ret);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& A, const T2& B, Matrix& ret)
    {
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<T1,Dense_2>(A,Bc, ret);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& A, const T2& B, Matrix& ret)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<Dense_1,T2>(Ac,B, ret);
    };
};

void construct_householder_impl(ret_construct_householder& ret, const Matrix& y)
{
    Integer M   = y.rows();
    Integer N   = y.cols();

    if ( (M == 1 && N >= 1 || N == 1 && M >= 1) == false)
        throw error::invalid_size(M,N);

    return construct_householder_vis::make<const Matrix&>(y,ret);
};

void householder_to_unitary_impl(unitary_matrix& ret, const Matrix& T, const Matrix& tau, Integer cols)
{
    Integer N   = T.rows();
    Integer K   = T.cols();

    if (N < 1 || K > N)
        throw error::invalid_size(N,K);

    if ( tau.length() != K)
        throw error::invalid_size(tau.rows(),tau.cols());

    if (cols < 1 || cols > N)
        throw error::invalid_holseholder_col(cols,N);

    return householder_to_unitary_vis::make(T, tau, ret, cols);
};

void householder_WY_impl(Matrix& ret, const Matrix& T, const Matrix& tau)
{
    Integer N   = T.rows();
    Integer K   = T.cols();

    if (N < 1 || K > N)
        throw error::invalid_size(N,K);

    if ( tau.length() != K)
        throw error::invalid_size(tau.rows(),tau.cols());

    return householder_WY_vis::make(T, tau, ret);
};

}};
namespace matcl
{

namespace md = matcl::details;

ret_construct_householder matcl::construct_householder(const Matrix& y0)
{
    Matrix y(y0);

    ret_construct_householder ret;
    details::construct_householder_impl(ret,y);
    return ret;
};

ret_construct_householder matcl::construct_householder(Matrix&& y0)
{
    Matrix y(std::move(y0));

    ret_construct_householder ret;
    details::construct_householder_impl(ret,y);
    return ret;
};

unitary_matrix matcl::householder_to_unitary(const Matrix& T0, const Matrix& tau)
{
    Matrix T(T0);

    unitary_matrix ret;
    details::householder_to_unitary_impl(ret, T, tau, T.rows());
    return ret;
};
unitary_matrix matcl::householder_to_unitary(const Matrix& T0, const Matrix& tau, Integer cols)
{
    Matrix T(T0);

    unitary_matrix ret;
    details::householder_to_unitary_impl(ret, T, tau, cols);
    return ret;
};

unitary_matrix matcl::householder_to_unitary(Matrix&& T0, const Matrix& tau)
{
    Matrix T(std::move(T0));

    unitary_matrix ret;
    details::householder_to_unitary_impl(ret, T, tau, T.rows());
    return ret;
};
unitary_matrix matcl::householder_to_unitary(Matrix&& T0, const Matrix& tau, Integer cols)
{
    Matrix T(std::move(T0));

    unitary_matrix ret;
    details::householder_to_unitary_impl(ret, T, tau, cols);
    return ret;
};
Matrix matcl::householder_WY(const Matrix& T0, const Matrix& tau)
{
    Matrix T(T0);

    Matrix ret;
    details::householder_WY_impl(ret, T, tau);
    return ret;
};
Matrix matcl::householder_WY(Matrix&& T0, const Matrix& tau)
{
    Matrix T(std::move(T0));

    Matrix ret;
    details::householder_WY_impl(ret, T, tau);
    return ret;
};

};
