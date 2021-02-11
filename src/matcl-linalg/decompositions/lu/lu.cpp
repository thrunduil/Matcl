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

#include "matcl-linalg/decompositions/lu.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/linear_eq/linsolve.h"

#include "matcl-linalg/decompositions/lu/lu_impl.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include "matcl-linalg/linear_eq/linsolve_objects_decomp.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "lu_superlu.h"

namespace matcl { namespace details
{

//---------------------------------------------------------------
//                  ilu_impl
//---------------------------------------------------------------

template<class V>
struct ilu_impl
{
    using Mat   = raw::Matrix<V,struct_sparse>;
    static void eval(const Mat& A, lu_return_type& ret, const options& opts)
    {
        superlu_wrap<V>::eval_ilu(A, ret, opts);
    };
};

template<>
struct ilu_impl<Object>
{
    using Mat   = raw::Matrix<Object,struct_sparse>;
    static void eval(const Mat&, lu_return_type&, const options&)
    {
        throw error::object_value_type_not_allowed("ilu");
    };
};

template<>
struct ilu_impl<Integer>
{
    using Mat   = raw::Matrix<Integer,struct_sparse>;
    static void eval(const Mat& A, lu_return_type& ret, const options& opts)
    {
        using Mat_R = raw::Matrix<Real,struct_sparse>;

        Mat_R Ar    = raw::converter<Mat_R,Mat>::eval(A);
        return ilu_impl<Real>::eval(Ar, ret, opts);
    };
};

//---------------------------------------------------------------
//                  linsolve_tridiag_impl
//---------------------------------------------------------------
template<class V>
struct linsolve_tridiag_impl
{
    using VR    = typename md::real_type<V>::type;

    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_B = raw::Matrix<V,struct_banded>;
    using Mat_I = raw::Matrix<Integer,struct_dense>;

    static void eval(linsolve_obj& ret, const Mat& Dm1, const Mat& D0, const Mat& Dp1, const options& opts)
    { 
        Integer N           = D0.length();
        const V* ptr_D      = D0.ptr();
        const V* ptr_Em1    = Dm1.ptr();
        const V* ptr_Ep1    = Dp1.ptr();

        Mat_B A(D0.get_type(), N, N, -1, 1);

        V* ptr_A0       = A.rep_ptr() + A.first_elem_diag(0);
        V* ptr_AL       = A.rep_ptr() + A.first_elem_diag(-1);
        V* ptr_AU       = A.rep_ptr() + A.first_elem_diag(1);
        Integer A_ld    = A.ld();

        for (Integer i = 0; i < N-1; ++i)
        {
            ptr_A0[0]   = ptr_D[i];
            ptr_AL[0]   = ptr_Em1[i];
            ptr_AU[0]   = ptr_Ep1[i];

            ptr_A0      += A_ld;
            ptr_AL      += A_ld;
            ptr_AU      += A_ld;
        };
        ptr_A0[0]       = ptr_D[N-1];

        Mat D0_f        = D0.make_unique();
        Mat Dm1_f       = Dm1.make_unique();
        Mat Dp1_f       = Dp1.make_unique();
        Mat Dp2_f(D0.get_type(), std::max(N-2,1), 1);
        Mat_I ipiv(D0.get_type(), N, 1);

        V* ptr_Df       = D0_f.ptr();
        V* ptr_Em1f     = Dm1_f.ptr();
        V* ptr_Ep1f     = Dp1_f.ptr();
        V* ptr_Ep2f     = Dp2_f.ptr();
        Integer* ptr_piv = ipiv.ptr();

        Integer info;
        lapack::gttrf(N, lap(ptr_Em1f), lap(ptr_Df), lap(ptr_Ep1f), lap(ptr_Ep2f), ptr_piv, info);

        if (info < 0)
            throw error::error_general("invalid argument passed to gttrf");

        Real tol        = opts.get_option<Real>(opt::linsolve::tol_sing());

        if (info > 0 && tol == 0.0)
            throw error::error_singular();

        bool isv    = D0_f.all_finite() && Dm1_f.all_finite() && Dp1_f.all_finite() && Dp2_f.all_finite();

        if (isv == false)
        {          
            value_code vc   = matrix_traits::value_code<VR>::value;
            using data_ptr  = linsolve_obj::linsolve_data_ptr;

            ret = linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, 
                        ti::convert_ti_object<VR>(D0.get_type()))));
            return;
        };

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        ret = linsolve_obj(data_ptr(new details::linsolve_obj_tridiag_fac<V>
                                    (A, D0_f, Dm1_f, Dp1_f, Dp2_f, ipiv, opts)));
    };
};

struct unary_visitor_lu : public extract_type_switch<void,unary_visitor_lu,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, const matcl::options& opts, lu_return_type& ret)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return details::lu<V,S>::eval(ret, mat,opts);
    };
    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, const matcl::options& opts, 
                            lu_return_type& ret)
    {
        using dense_mat = raw::Matrix<T,struct_dense>;
        return eval<dense_mat>(handle,dense_mat(ti::get_ti(mat),mat,1,1),opts, ret);
    };
};

struct linsolve_lu_tridiag_vis : public extract_type_switch<void, linsolve_lu_tridiag_vis,true>
{
    template<class V>
    static void eval_impl(linsolve_obj& ret, const Matrix& Dm1_0, const Matrix& D0_0, const Matrix& Dp1_0,
                          const options& opts)
    {
        using VT    = typename md::unify_types<V, Float>::type;
        using VR    = typename real_type_int_real<VT>::type;
        using S     = struct_dense;        
        using Mat   = raw::Matrix<VT,struct_dense>;

        const Mat& Dm1  = Dm1_0.impl<Mat>().make_explicit();
        const Mat& D0   = D0_0.impl<Mat>().make_explicit();
        const Mat& Dp1  = Dp1_0.impl<Mat>().make_explicit();

        bool isv    = Dm1.all_finite() && D0.all_finite() && Dp1.all_finite();

        Integer N   = D0.length();

        if (isv == false)
        {          
            value_code vc   = matrix_traits::value_code<V>::value;

            using data_ptr = linsolve_obj::linsolve_data_ptr;
            ret = linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, 
                        ti::convert_ti_object<VR>(D0.get_type()))));
            return;
        };

        return details::linsolve_tridiag_impl<VT>::eval(ret, Dm1, D0, Dp1, opts);
    };

    template<class T>
    static void eval(const Matrix&, const T&, linsolve_obj& ret, const Matrix& Dm1, const Matrix& D0, 
                     const Matrix& Dp1, const options& opts)
    {
        using V     = typename T::value_type;
        return eval_impl<V>(ret, Dm1, D0, Dp1, opts);
    };
    
    template<class T>
    static void eval_scalar(const Matrix&, const T&, linsolve_obj& ret, const Matrix& Dm1, const Matrix& D0, 
                     const Matrix& Dp1, const options& opts)
    {
        return eval_impl<T>(ret, Dm1, D0, Dp1, opts);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, linsolve_obj&, const Matrix&, 
                     const Matrix&, const Matrix&, const options&)
    {
        throw error::object_value_type_not_allowed("linsolve_tridiag");
    };

    static void eval_scalar(const Matrix&, const Object&, linsolve_obj&, const Matrix&, const Matrix&, const Matrix&,
                            const options&)
    {
        throw error::object_value_type_not_allowed("linsolve_tridiag");
    };
};

struct visit_ilu : public md::extract_type_switch<void, visit_ilu,true>
{
    using base_type = md::extract_type_switch<void, visit_ilu,true>;

    template<class V>
    static void eval(const matcl::Matrix&, const raw::Matrix<V,struct_sparse>& mat, lu_return_type& ret,
                     const options& opts)
    {
        return ilu_impl<V>::eval(mat, ret, opts);
    };

    template<class Mat, class ... Arg>
    static void eval(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };

    template<class Mat, class ... Arg>
    static void eval_scalar(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };
};

};};

namespace matcl
{

lu_return_type matcl::lu(const Matrix& A0, const matcl::options& opts)
{
    Matrix A(A0);

    lu_return_type ret;
    details::unary_visitor_lu::make<const Matrix&>(A,opts, ret);
    return ret;
};
lu_return_type matcl::lu(Matrix&& A0, const matcl::options& opts)
{
    Matrix A(std::move(A0));

    lu_return_type ret;
    details::unary_visitor_lu::make<const Matrix&>(A,opts, ret);
    return ret;
};

lu_return_type matcl::ilu(const Matrix& A0, const options& opts)
{
    Matrix A(A0);
        
    struct_code st  = A.get_struct_code();
    if (st == struct_code::struct_dense || st == struct_code::struct_banded)
    {
        return lu(A, opts);
    };
    
    lu_return_type ret;
    details::visit_ilu::make<const Matrix&>(A, ret, opts);
    return ret;
};

lu_return_type matcl::ilu(Matrix&& A0, const options& opts)
{
    Matrix A(std::move(A0));    

    struct_code st  = A.get_struct_code();

    if (st == struct_code::struct_dense || st == struct_code::struct_banded)
        return lu(std::move(A), opts);

    lu_return_type ret;
    details::visit_ilu::make<const Matrix&>(A, ret, opts);
    return ret;
};

linsolve_obj matcl::linsolve_lu(const Matrix& A, const Matrix& L, const Matrix& U, const permvec& p, 
                                   const permvec& q, const options& opts)
{
    if (matcl::is_tril(L) == false)
        throw error::invalid_lu_factors();

    if (matcl::is_triu(U) == false)
        throw error::invalid_lu_factors();

    if (L.rows() != p.length())
        throw error::invalid_lu_factors();

    if (U.cols() != q.length())
        throw error::invalid_lu_factors();
    if (A.rows() != L.rows() )
        throw error::invalid_lu_factors();
    if (A.cols() != U.cols())
        throw error::invalid_lu_factors();

    if (L.rows() != L.cols())
        throw error::square_matrix_required(L.rows(), L.cols());

    if (U.rows() != U.cols())
        throw error::square_matrix_required(U.rows(), U.cols());

    if (L.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_lu");
    if (U.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_lu");

    if (L.rows() == 0)
    {
        value_code vc   = A.get_value_code();
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
        vc              = matrix_traits::unify_value_types(vc, L.get_value_code());
        vc              = matrix_traits::unify_value_types(vc, U.get_value_code());

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type() )));
    };

    bool isv    = L.all_finite() && U.all_finite() && A.all_finite();

    if (isv == false)
    {          
        value_code vc   = A.get_value_code();
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
        vc              = matrix_traits::unify_value_types(vc, L.get_value_code());
        vc              = matrix_traits::unify_value_types(vc, U.get_value_code());

        return linsolve_nan(L.rows(), vc);
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_lu_factors(A, L,U,p,q, opts)));
};

static void linsolve_lu_tridiag_impl(linsolve_obj& ret, const Matrix& D_m1, const Matrix& D_0, const Matrix& D_p1,
                                     const options& opts)
{
    if (D_0.is_vector() == false)
        throw error::error_size_tridiag_nsym(D_0.rows(),D_0.cols(), D_m1.rows(), D_m1.cols(),
                                        D_p1.rows(), D_p1.cols());

    Integer N   = D_0.length();

    if (N > 0 && D_m1.length() < N-1)
        throw error::error_size_tridiag_nsym(D_0.rows(),D_0.cols(), D_m1.rows(), D_m1.cols(),
                                        D_p1.rows(), D_p1.cols());

    if (N > 0 && D_p1.length() < N-1)
        throw error::error_size_tridiag_nsym(D_0.rows(),D_0.cols(), D_m1.rows(), D_m1.cols(),
                                        D_p1.rows(), D_p1.cols());

    value_code vc_1 = D_0.get_value_code();
    value_code vc_2 = D_m1.get_value_code();
    value_code vc_3 = D_p1.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_1, vc_2);
    vc              = matrix_traits::unify_value_types(vc, vc_3);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        ret = linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, D_0.get_type())));
        return;
    };

    Matrix dum = zeros(0,0,vc);
    return details::linsolve_lu_tridiag_vis::make<const Matrix&>(dum, ret, D_m1, D_0, D_p1, opts);
};

linsolve_obj matcl::linsolve_tridiag(const Matrix& D_m10, const Matrix& D00, const Matrix& D_p10, 
                                     const options& opts)
{
    Matrix D_m1(D_m10);
    Matrix D0(D00);
    Matrix D_p1(D_p10);

    linsolve_obj ret;
    linsolve_lu_tridiag_impl(ret, D_m1, D0, D_p1, opts);

    return ret;
};

};