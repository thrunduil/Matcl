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

#include "matcl-linalg/decompositions/chol.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/decompositions/chol_impl.h"
#include "matcl-linalg/decompositions/chol_rr_impl.h"
#include "matcl-linalg/decompositions/chol_update.h"
#include "matcl-linalg/linear_eq/linsolve_objects_decomp.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/decompositions/balancing.h"

namespace matcl { namespace details
{

//---------------------------------------------------------------
//                  ldl_tridiag_impl
//---------------------------------------------------------------
template<class V>
struct ldl_tridiag_impl
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_R = raw::Matrix<VR,struct_dense>;
    using Mat_B = raw::Matrix<V,struct_banded>;

    static void eval(mat_tup_2& ret, Mat& D0_0, Mat& D1, bool upper, Integer& k)
    {
        using VR    = typename md::real_type<V>::type;
        using Mat_R = raw::Matrix<VR,struct_dense>;

        Mat_R D0        = Mat_R(raw::converter<Mat_R, Mat>::eval(D0_0), Mat_R::copy_is_safe());

        Integer N       = D0.length();
        VR* ptr_D       = D0.ptr();
        V* ptr_E        = D1.ptr();

        Integer info;
        lapack::pttrf(N, lap(ptr_D), lap(ptr_E), info);

        if (info < 0)
            throw error::error_general("invalid argument passed to pttrf");

        if (info > 0)
        {
            if (k == 0)
                throw error::error_singular();
            else
                k = info - 1;
        }
        else
        {
            k       = N;
        };

        Matrix D(D0, false);
        D   = bdiag(D);

        using Mat_B = raw::Matrix<V,struct_banded>;

        if (upper == true)
        {
            Mat_B S0(D1.get_type(), N, N, 0, 1);

            V* ptr_S_0      = S0.rep_ptr() + S0.first_elem_diag(0);
            V* ptr_S_1      = S0.rep_ptr() + S0.first_elem_diag(1);
            Integer S_ld    = S0.ld();

            for (Integer i = 0; i < N; ++i)
            {
                ptr_S_0[0]  = V(1.0);
                ptr_S_0     += S_ld;
            }

            for (Integer i = 0; i < N-1; ++i)
            {
                ptr_S_1[0]  = conj(ptr_E[i]);
                ptr_S_0     += S_ld;
            }

            Matrix S(S0,false);
            ret = mat_tup_2(S,D);
            return;
        }
        else
        {
            Mat_B S0(D1.get_type(), N, N, -1, 0);

            V* ptr_S_0      = S0.rep_ptr() + S0.first_elem_diag(0);
            V* ptr_S_1      = S0.rep_ptr() + S0.first_elem_diag(-1);
            Integer S_ld    = S0.ld();

            for (Integer i = 0; i < N; ++i)
            {
                ptr_S_0[0]  = V(1.0);
                ptr_S_0     += S_ld;
            }

            for (Integer i = 0; i < N-1; ++i)
            {
                ptr_S_1[0]  = ptr_E[i];
                ptr_S_0     += S_ld;
            }

            Matrix S(S0,false);
            ret = mat_tup_2(S,D);
            return;
        };
    };

    static void eval_linsolve(linsolve_obj& ret, const Mat_R& D0, const Mat& D1, const options& opts)
    {                
        Integer N       = D0.length();
        const VR* ptr_D = D0.ptr();
        const V* ptr_E  = D1.ptr();        

        Mat_B A(D0.get_type(), N, N, -1, 1);

        V* ptr_A0       = A.rep_ptr() + A.first_elem_diag(0);
        V* ptr_AL       = A.rep_ptr() + A.first_elem_diag(-1);
        V* ptr_AU       = A.rep_ptr() + A.first_elem_diag(1);
        Integer A_ld    = A.ld();

        for (Integer i = 0; i < N-1; ++i)
        {
            ptr_A0[0]   = ptr_D[i];
            ptr_AL[0]   = ptr_E[i];
            ptr_AU[0]   = conj(ptr_E[i]);

            ptr_A0      += A_ld;
            ptr_AL      += A_ld;
            ptr_AU      += A_ld;
        };
        ptr_A0[0]       = ptr_D[N-1];

        Mat_R D0_f      = D0.make_explicit().make_unique();
        Mat D1_f        = D1.make_explicit().make_unique();
        VR* ptr_Df      = D0_f.ptr();
        V* ptr_Ef       = D1_f.ptr();

        Integer info;
        lapack::pttrf(N, lap(ptr_Df), lap(ptr_Ef), info);

        if (info < 0)
            throw error::error_general("invalid argument passed to pttrf");

        Real tol        = opts.get_option<Real>(opt::linsolve::tol_sing());

        if (info > 0 && (info < N || tol == 0.0))
            throw error::error_singular();

        bool isv    = D0_f.all_finite() && D1_f.all_finite();

        if (isv == false)
        {          
            value_code vc   = matrix_traits::value_code<VR>::value;
            using data_ptr  = linsolve_obj::linsolve_data_ptr;

            ret = linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, 
                        ti::convert_ti_object<VR>(D0.get_type()))));
            return;
        };

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        ret = linsolve_obj(data_ptr(new details::linsolve_obj_chol_tridiag_fac<V>(A, D0_f, D1_f, opts)));
    };

    static void eval_linsolve_fac(linsolve_obj& ret, const Mat_B& A, const Mat_R& D0, const Mat& D1, 
                                  const options& opts)
    {
        bool isv    = A.all_finite() && D0.all_finite() && D1.all_finite();

        if (isv == false)
        {
            value_code vc   = matrix_traits::value_code<VR>::value;
            using data_ptr = linsolve_obj::linsolve_data_ptr;
            ret = linsolve_obj(data_ptr(new details::linsolve_obj_nan(D0.length(), vc, 
                            ti::convert_ti_object<VR>(D0.get_type()) ))); 
            return;
        };

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        ret = linsolve_obj(data_ptr(new details::linsolve_obj_chol_tridiag_fac<V>(A, D0, D1, opts)));
    };
};

//---------------------------------------------------------------
//                  visitors
//---------------------------------------------------------------
struct unary_visitor_chol : public extract_type_switch<void,unary_visitor_chol,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, chol2_return_type& ret, bool upper, 
                     const options& opts)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        using VR    = typename real_type_int_real<V>::type;

        bool isv    = mat.all_finite();
        
        if (isv == false)
        {
            Matrix S    = md::make_nan_matrix<VR>(mat.rows(), mat.rows());
            ret = chol2_return_type(S,permvec::identity(mat.rows()), 0);
            return;
        };

        return details::chol<V,S>::eval(ret, mat, upper, opts);
    };
    
    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, chol2_return_type& ret, bool,
                            const options&)
    {
        using VR    = typename real_type_int_real<T>::type;

        bool isv = matcl::is_nan(mat) == false && matcl::is_inf(mat) == false;
        
        if (isv == false)
        {
            permvec p = details::pv_constructor::make(1);
            ret = chol2_return_type(matcl::constants::nan<VR>(), p, 0);
            return;
        };

        auto v = matcl::real(mat);

        if(v > 0)
            ret = chol2_return_type(raw::details::sqrt_helper<T>::eval(v), details::pv_constructor::make(1), 1);
        else
            ret = chol2_return_type(v, details::pv_constructor::make(1), 0);

        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, chol2_return_type&, bool, const options&)
    {
        throw error::object_value_type_not_allowed("chol");
    };

    static void eval_scalar(const Matrix&, const Object&, chol2_return_type&, bool, const options&)
    {
        throw error::object_value_type_not_allowed("chol");
    };
};

struct ldl_tridiag_vis : public extract_type_switch<void, ldl_tridiag_vis,true>
{
    template<class V>
    static void eval_impl(mat_tup_2& ret, const Matrix& D0_0, const Matrix& D1_0,
                          bool upper, Integer& k)
    {
        using VT    = typename md::unify_types<V, Float>::type;
        using VR    = typename real_type_int_real<VT>::type;
        using S     = struct_dense;        
        using Mat   = raw::Matrix<VT,struct_dense>;

        Mat D0      = convert(D0_0, Mat::matrix_code).get_impl<Mat>().make_explicit().make_unique();
        Mat D1      = convert(D1_0, Mat::matrix_code).get_impl<Mat>().make_explicit().make_unique();

        bool isv    = D0.all_finite() && D1.all_finite();
        
        Integer N   = D0.length();

        if (isv == false)
        {            
            Matrix D    = md::make_nan_matrix<VR>(N, N);
            Matrix S    = D;
                
            ret         = mat_tup_2(S,D);
            return;
        };

        return details::ldl_tridiag_impl<VT>::eval(ret, D0, D1, upper, k);
    };

    template<class T>
    static void eval(const Matrix&, const T&, mat_tup_2& ret, const Matrix& D0, 
                     const Matrix& D1, bool upper, Integer& k)
    {
        using V     = typename T::value_type;
        return eval_impl<V>(ret, D0, D1, upper, k);
    };
    
    template<class T>
    static void eval_scalar(const Matrix&, const T&, mat_tup_2& ret, const Matrix& D0, 
                            const Matrix& D1,bool upper, Integer& k)
    {
        return eval_impl<T>(ret, D0, D1, upper, k);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, mat_tup_2&, const Matrix&, 
                     const Matrix&, bool, Integer&)
    {
        throw error::object_value_type_not_allowed("ldl_tridiag");
    };

    static void eval_scalar(const Matrix&, const Object&, mat_tup_2&, const Matrix&, const Matrix&, 
                            bool, Integer&)
    {
        throw error::object_value_type_not_allowed("ldl_tridiag");
    };
};

struct linsolve_tridiag_vis : public extract_type_switch<void, linsolve_tridiag_vis,true>
{
    template<class V>
    static void eval_impl(linsolve_obj& ret, const Matrix& D0_0, const Matrix& D1_0, const options& opts)
    {
        using VT    = typename md::unify_types<V, Float>::type;
        using VR    = typename real_type_int_real<VT>::type;
        using S     = struct_dense;        
        using Mat   = raw::Matrix<VT,struct_dense>;        
        using Mat_R = raw::Matrix<VR,struct_dense>;        

        const Mat_R& D0     = convert(D0_0, Mat_R::matrix_code).get_impl<Mat_R>().make_explicit();
        const Mat& D1       = convert(D1_0, Mat::matrix_code).get_impl<Mat>().make_explicit();

        bool isv    = D0.all_finite() && D1.all_finite();

        Integer N   = D0.length();

        if (isv == false)
        {          
            value_code vc   = matrix_traits::value_code<V>::value;

            using data_ptr = linsolve_obj::linsolve_data_ptr;
            ret = linsolve_obj(data_ptr(new details::linsolve_obj_nan(N,vc, 
                        ti::convert_ti_object<VR>(D0.get_type()))));
            return;
        };

        return details::ldl_tridiag_impl<VT>::eval_linsolve(ret, D0, D1, opts);
    };

    template<class T>
    static void eval(const Matrix&, const T&, linsolve_obj& ret, const Matrix& D0, const Matrix& D1, 
                     const options& opts)
    {
        using V     = typename T::value_type;
        return eval_impl<V>(ret, D0, D1, opts);
    };
    
    template<class T>
    static void eval_scalar(const Matrix&, const T&, linsolve_obj& ret, const Matrix& D0, 
                            const Matrix& D1, const options& opts)
    {
        return eval_impl<T>(ret, D0, D1, opts);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, linsolve_obj&, const Matrix&, 
                     const Matrix&, const options&)
    {
        throw error::object_value_type_not_allowed("linsolve_chol_tridiag");
    };

    static void eval_scalar(const Matrix&, const Object&, linsolve_obj&, const Matrix&, const Matrix&, const options&)
    {
        throw error::object_value_type_not_allowed("linsolve_chol_tridiag");
    };
};

struct linsolve_tridiag_fac_vis : public extract_type_switch<void, linsolve_tridiag_fac_vis,true>
{
    template<class V>
    static void eval_impl(linsolve_obj& ret, const Matrix& A, const Matrix& D, const Matrix& S, bool upper,
                          const options& opts)
    {
        using VT    = typename md::unify_types<V, Float>::type;
        using VR    = typename real_type_int_real<VT>::type;
        using Mat   = raw::Matrix<VT,struct_dense>;
        using Mat_R = raw::Matrix<VR,struct_dense>;
        using Mat_B = raw::Matrix<VT,struct_banded>;

        Matrix D0   = D.diag(0);
        Matrix D1   = (upper == false) ? S.diag(-1) : S.diag(1);

        const Mat_R& D0m    = convert(D0, Mat_R::matrix_code).get_impl<Mat_R>().make_explicit();
        const Mat& D1m      = convert(D1, Mat::matrix_code).get_impl<Mat>().make_explicit();
        const Mat_B& Am     = convert(A, Mat_B::matrix_code).get_impl<Mat_B>();

        return details::ldl_tridiag_impl<VT>::eval_linsolve_fac(ret, Am, D0m, D1m, opts);
    };

    template<class T>
    static void eval(const Matrix&, const T&, linsolve_obj& ret, const Matrix& A, const Matrix& D0, 
                     const Matrix& D1, bool upper, const options& opts)
    {
        using V     = typename T::value_type;
        return eval_impl<V>(ret, A, D0, D1, upper, opts);
    };
    
    template<class T>
    static void eval_scalar(const Matrix&, const T&, linsolve_obj& ret, const Matrix& A, const Matrix& D0, 
                            const Matrix& D1, bool upper, const options& opts)
    {
        return eval_impl<T>(ret, A, D0, D1, upper, opts);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, linsolve_obj&, const Matrix&, const Matrix&, 
                     const Matrix&, bool, const options&)
    {
        throw error::object_value_type_not_allowed("linsolve_chol_tridiag");
    };

    static void eval_scalar(const Matrix&, const Object&, linsolve_obj&, const Matrix&, const Matrix&, 
                            const Matrix&, bool, const options&)
    {
        throw error::object_value_type_not_allowed("linsolve_chol_tridiag");
    };
};

struct unary_visitor_chol_rr : public extract_type_switch<void,unary_visitor_chol_rr,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, chol_rr_return_type& ret, bool upper, Real tol)
    {
        using V         = typename T::value_type;
        using S         = typename T::struct_type;
        using VR        = typename md::real_type_int_real<V>::type;

        bool isv        = mat.all_finite();
        
        if (isv == false)
        {
            Matrix S    = md::make_nan_matrix<VR>(mat.rows(), mat.rows());
            ret         = chol_rr_return_type(S,permvec::identity(mat.rows()), 0);
            return;
        };

        return details::chol_rr<V,S>::eval(ret, mat, upper, VR(tol));
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, chol_rr_return_type& ret, bool, Real tol)
    {
        using VR        = typename md::real_type_int_real<T>::type;

        bool isv = matcl::is_nan(mat) == false && matcl::is_inf(mat) == false;
        
        if (isv == false)
        {
            permvec p = details::pv_constructor::make(1);
            ret = chol_rr_return_type(matcl::constants::nan<VR>(), p, 0);
            return;
        };

        auto v      = matcl::real(mat);
        tol         = std::max(tol, 0.);

        if(v >= tol)
        {
            Matrix S    = raw::details::sqrt_helper<T>::eval(v);
            permvec p   = details::pv_constructor::make(1);   
            ret = chol_rr_return_type(S,p, 1);
            return;
        }
        else
        {
            permvec p = details::pv_constructor::make(1);
            ret = chol_rr_return_type(0.,p, 0);
            return;
        }
    };

    static void eval_scalar(const Matrix&, const Object&, chol_rr_return_type&, bool, Real)
    {
        throw error::object_value_type_not_allowed("chol");
    };
    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, chol_rr_return_type&, bool, Real)
    {
        throw error::object_value_type_not_allowed("chol");
    };
};

struct visitor_chol_remove : public extract_type_switch<void,visitor_chol_remove,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, matcl::Matrix& ret, bool upper, Integer k)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::chol_update<V,S>::eval_remove(ret, mat, upper, k);
    };

    template<class T>
    static void eval(const Matrix&, const T& mat, matcl::Matrix& ret, bool upper, const matcl::Matrix& K)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::chol_update<V,S>::eval_remove(ret, mat, upper, K);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix& A, const T&, matcl::Matrix& ret, bool, Integer)
    {
        ret = matcl::zeros(0,0, A.get_value_code());
    };
    template<class T>
    static void eval_scalar(const matcl::Matrix& A, const T&, matcl::Matrix& ret, bool, 
                                     const matcl::Matrix&)
    {
        //K matrix should contain valid row indices
        ret = matcl::zeros(0,0, A.get_value_code());
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, matcl::Matrix&, bool, Integer)
    {
        throw error::object_value_type_not_allowed("chol_remove");
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, matcl::Matrix&, bool, const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("chol_remove");
    };

    static void eval_scalar(const Matrix&, const Object&, matcl::Matrix&, bool, Integer)
    {
        throw error::object_value_type_not_allowed("chol_remove");
    };
    static void eval_scalar(const Matrix&, const Object&, matcl::Matrix&, 
                                     bool, const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("chol_remove");
    };
};

struct visitor_chol_update : public extract_type_switch<void,visitor_chol_update,true>
{
    template<class T>
    static void eval(const Matrix& h, const T& mat, Matrix& ret, const Matrix& w, bool upper, 
                     const Matrix& sigma)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        value_code vt1  = matrix_traits::value_code<V>::value;
        value_code vt2  = w.get_value_code();
        value_code vt0  = matrix_traits::unify_value_types(vt1, vt2);
        value_code vt   = matrix_traits::unify_value_types(vt0, value_code::v_float);

        if (vt != vt1)
        {
            //conversion is required
            struct_code st  = h.get_struct_code();
            mat_code mc     = matrix_traits::get_matrix_type(vt, st);
            Matrix Ac       = convert(h, mc);
            return visitor_chol_update::make<const Matrix&>(Ac,ret,w,upper,sigma);
        };

        return details::chol_update<V,S>::eval_update(ret, mat, w, upper, sigma);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix&, const T& A, Matrix& ret, const Matrix& w, 
                            bool, const Matrix& Sigma)
    {
        value_code vc = w.get_value_code();

        using TR0   = typename md::real_type<T>::type;

        if (vc == value_code::v_float || vc == value_code::v_float_complex)
        {
            using TR    = typename md::unify_types<TR0,Float>::type;
            TR sig      = Sigma.get_scalar<TR>();
            TR val      = abs2(A) + sig * abs2(w).get_scalar<TR>();

            if (val < 0)
                throw error::error_invalid_chol_update_nonposdef();

            ret = sqrt(val);
        }
        else
        {   
            using TR    = typename md::unify_types<TR0,Real>::type;
            TR sig      = Sigma.get_scalar<TR>();
            TR val      = abs2(A) + sig * abs2(w).get_scalar<TR>();

            if (val < 0)
                throw error::error_invalid_chol_update_nonposdef();

            ret = sqrt(val);
        };
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, 
                     const Matrix&, bool, const Matrix& )
    {
        throw error::object_value_type_not_allowed("chol_update");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, const Matrix&, bool, 
                            const Matrix& )
    {
        throw error::object_value_type_not_allowed("chol_update");
    };
};

};};

namespace matcl
{

static void chol_impl(chol2_return_type& ret, const Matrix &A, bool upper, const options& opts)
{
    if (A.rows() != A.cols())
        throw error::error_nonsymh_chol();

    if (A.structural_nnz() == 0)
    {
        permvec p = permvec::identity(A.rows());
        ret = chol2_return_type(A, p, 0);
        return;
    };

    return details::unary_visitor_chol::make<const Matrix&>(A, ret, upper, opts);
}

static void ldl_tridiag_impl(mat_tup_2& ret, const Matrix& D0, const Matrix& D1, bool upper, Integer& k)
{
    if (D0.is_vector() == false || D0.is_vector() == false)
        throw error::error_size_tridiag_sym(D0.rows(),D0.cols(),D1.rows(), D1.cols());

    Integer N   = D0.length();

    if (N > 0 && D1.length() < N-1)
        throw error::error_size_tridiag_sym(D0.rows(), D0.cols(), D1.rows(), D1.cols());

    value_code vc_1 = D0.get_value_code();
    value_code vc_2 = D1.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_1, vc_2);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (N == 0)
    {
        vc          = matrix_traits::real_value_type(vc);

        Matrix S    = zeros(0,0,vc);
        Matrix D    = S;
        ret         = mat_tup_2(S,D);
        return;
    };

    Matrix dum = zeros(0,0,vc);
    return details::ldl_tridiag_vis::make<const Matrix&>(dum, ret, D0, D1, upper, k);
};

static void linsolve_tridiag_impl(linsolve_obj& ret, const Matrix& D0, const Matrix& D1, const options& opts)
{
    if (D0.is_vector() == false)
        throw error::error_size_tridiag_sym(D0.rows(),D0.cols(),D1.rows(), D1.cols());

    Integer N   = D0.length();

    if (N > 0 && D1.length() < N-1)
        throw error::error_size_tridiag_sym(D0.rows(), D0.cols(), D1.rows(), D1.cols());

    value_code vc_1 = D0.get_value_code();
    value_code vc_2 = D1.get_value_code();
    value_code vc   = matrix_traits::unify_value_types(vc_1, vc_2);
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        ret = linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, D0.get_type())));
        return;
    };

    Matrix dum = zeros(0,0,vc);
    return details::linsolve_tridiag_vis::make<const Matrix&>(dum, ret, D0, D1, opts);
};

chol_return_type matcl::chol(const Matrix& A0, bool upper, const options& opts)
{
    //increase refcount
    Matrix A(A0);

    chol2_return_type ret;
    chol_impl(ret, A, upper, opts);

    Integer N = ret.get<3>();
    if (N < A.rows())
        throw error::error_nonposdef(false);

    return chol_return_type(ret.get<1>(), ret.get<2>());
};

chol_return_type matcl::chol(Matrix&& A0, bool upper, const options& opts)
{
    Matrix A(std::move(A0));

    chol2_return_type ret;
    chol_impl(ret, A, upper, opts);

    Integer N = ret.get<3>();
    if (N < A.rows())
        throw error::error_nonposdef(false);

    return chol_return_type(ret.get<1>(), ret.get<2>());
};

ldl_tridiag_return matcl::ldl_tridiag(const Matrix& D0_0, const Matrix& D1_0, bool upper)
{
    Matrix D0(D0_0);
    Matrix D1(D1_0);

    ldl_tridiag_return ret;
    Integer k   = 0;
    ldl_tridiag_impl(ret, D0, D1, upper, k);
    return ret;
};

ldl_tridiag2_return matcl::ldl_tridiag2(const Matrix& D0_0, const Matrix& D1_0, bool upper)
{
    Matrix D0(D0_0);
    Matrix D1(D1_0);

    ldl_tridiag_return ret;
    Integer k   = 1;
    ldl_tridiag_impl(ret, D0, D1, upper, k);

    return ldl_tridiag2_return(ret.get<1>(), ret.get<2>(), k);
};

chol2_return_type matcl::chol2(const Matrix& A0, bool upper, const options& opts)
{
    Matrix A(A0);

    chol2_return_type ret;
    chol_impl(ret, A, upper, opts);
    return ret;
};

chol2_return_type matcl::chol2(Matrix&& A0, bool upper, const options& opts)
{
    Matrix A(std::move(A0));

    chol2_return_type ret;
    chol_impl(ret, A, upper, opts);
    return ret;
};

static void chol_rr_impl(chol_rr_return_type& ret, const Matrix &A, bool upper, Real tol)
{
    if (A.rows() != A.cols())
        throw error::error_nonsymh_chol();

    if (A.structural_nnz() == 0)
    {
        permvec p = permvec::identity(A.rows());
        ret =  chol_rr_return_type(A,p,0);
        return;
    };

    return details::unary_visitor_chol_rr::make<const Matrix&>(A, ret, upper, tol);
}

chol_rr_return_type matcl::chol_rr(const Matrix& A0, bool upper, Real tol)
{
    //increase refcount
    Matrix A(A0);

    chol_rr_return_type ret;
    chol_rr_impl(ret, A,upper,tol);
    return ret;
};

chol_rr_return_type matcl::chol_rr(Matrix&& A0, bool upper, Real tol)
{
    //increase refcount
    Matrix A(std::move(A0));

    chol_rr_return_type ret;
    chol_rr_impl(ret, A, upper, tol);
    return ret;
};

//------------------------------------------------------------------------------------
//                      UPDATE
//------------------------------------------------------------------------------------
static void chol_remove_impl(Matrix& ret, const Matrix& T, bool upper, const matcl::Matrix& K)
{
    if (T.rows() != T.cols())
        throw error::error_invalid_chol_factor(upper);

    if ((upper == true && is_triu(T) == false) || (upper == false && is_tril(T) == false))
        throw error::error_invalid_chol_factor(upper);

    if (K.is_empty() == true)
    {
        ret = T;
        return;
    };
    
    //prepare vector of unique, sorted, valid positions
    matcl::Matrix K2        = sort(vec(K));
    Integer N               = K2.length();
    raw::Matrix<Integer,struct_dense> Ks(ti::ti_empty(), K2.length(), 1);

    const Integer* K2_ptr   = K2.get_array<Integer>();
    Integer* Ks_ptr         = Ks.ptr();

    Integer pos             = 0;
    Integer last_elem       = 0;
    for (Integer i = 0; i < N; ++i)
    {
        Integer k           = K2_ptr[i];
        
        if (k < 1)
            continue;
        if (k > T.rows())
            break;

        if (k > last_elem)
        {
            Ks_ptr[pos]     = k;
            last_elem       = k;
            ++pos;
        };
    };

    if (pos == 0)
    {
        ret = T;
        return;
    }

    if (pos == 1)
        return details::visitor_chol_remove::make<const Matrix&>(T, ret, upper, Ks_ptr[0]);

    matcl::Matrix K_in(Ks.resize(pos, 1), true);

    return details::visitor_chol_remove::make<const Matrix&>(T, ret, upper, K_in);
};

static void chol_update_impl(Matrix& ret, const Matrix& U, const Matrix& w, bool upper, 
                             const Matrix& sigma)
{
    if (U.rows() != U.cols())
        throw error::error_invalid_chol_factor(upper);

    if ((upper == true && is_triu(U) == false) || (upper == false && is_tril(U) == false))
        throw error::error_invalid_chol_factor(upper);

    if (w.rows() != U.rows())
        throw error::error_invalid_chol_update(U.rows(), w.rows(), w.cols());

    Integer K_vec   = w.cols();
    Integer K_r     = sigma.rows();
    Integer K_c     = sigma.cols();
    bool sigma_ok   = (K_r == 1 && (K_c == 1 || K_c == K_vec) )
                    || (K_c == 1 && (K_r == 1 || K_r == K_vec) );

    if (sigma_ok == false)
        throw error::error_invalid_update_sigma(K_vec, K_r, K_c);

    if (sigma.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("chol_update");

    if (matrix_traits::is_float_complex(sigma.get_value_code()) == true)
        throw error::error_complex_value_type_not_allowed();

    if (K_vec == 0)
    {
        ret = U;
        return;
    };

    return details::visitor_chol_update::make<const Matrix&>(U,ret,w,upper,sigma);
};

Matrix matcl::chol_update(const Matrix& U0, const Matrix& w, const Matrix& sigma, bool upper)
{
    //increase refcount
    Matrix U(U0);
    
    Matrix ret;
    chol_update_impl(ret, U, w, upper, sigma);
    return ret;
};

Matrix matcl::chol_update(Matrix&& U0, const Matrix& w, const Matrix& sigma, bool upper)
{
    Matrix U(std::move(U0));

    Matrix ret;
    chol_update_impl(ret, U, w, upper, sigma);
    return ret;
};

Matrix matcl::chol_remove(const Matrix& T0, const matcl::Matrix& K, bool upper)
{
    //increase refcount
    Matrix T(T0);

    Matrix ret;
    chol_remove_impl(ret, T,upper,K);
    return ret;
};

Matrix matcl::chol_remove(Matrix&& T0, const matcl::Matrix& K, bool upper)
{
    Matrix T(std::move(T0));

    Matrix ret;
    chol_remove_impl(ret, T,upper,K);
    return ret;
};

//------------------------------------------------------------------------------------
//                      LINSOLVE
//------------------------------------------------------------------------------------

linsolve_obj matcl::linsolve_chol(const Matrix& A, const Matrix& S, const permvec& p, bool upper,
                                  const options& opts)
{
    if (upper == false)
    {
        if (A.rows() != S.rows())
            throw error::invalid_cholecky_factors();
        if (A.cols() != S.rows())
            throw error::invalid_cholecky_factors();
    }
    else
    {
        if (A.rows() != S.cols())
            throw error::invalid_cholecky_factors();
        if (A.cols() != S.cols())
            throw error::invalid_cholecky_factors();
    }
    if (upper == true && matcl::is_triu(S) == false)
        throw error::invalid_cholecky_factors();

    if (upper == false && matcl::is_tril(S) == false)
        throw error::invalid_cholecky_factors();

    if (S.rows() != S.cols())
        throw error::square_matrix_required(S.rows(), S.cols());

    if (S.rows() != p.length())
        throw error::invalid_cholecky_factors();

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_chol");
    if (S.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_chol");

    Integer N = S.rows();

    if (N == 0)
    {
        value_code vc   = S.get_value_code();
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), value_code::v_float);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, S.get_type() )));
    };

    bool isv            = S.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = S.get_value_code();
        vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), value_code::v_float);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, S.get_type() )));
    };    

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_chol(A, S, p, upper, opts)));
};

linsolve_obj matcl::linsolve_chol_tridiag_fac(const Matrix& D0, const Matrix& D1,
                                    const Matrix& S, const Matrix& D, bool upper, const options& opts)
{
    if (D0.is_vector() == false)
        throw error::vector_required(D0.rows(), D0.cols());

    Integer N   = D0.length();

    if (N > 0 && D1.is_vector() == false)
        throw error::vector_required(D1.rows(), D1.cols());
    if (N == 0 && D1.length() != 0)
        throw error::invalid_size2(D1.rows(), D1.cols(), 0, 1);

    Matrix A = bzeros(N, N, -1, 1);

    if (N == 0)
        return linsolve_chol_tridiag_fac(A, S, D, upper, opts);

    A.diag(0) = D0;

    if (N > 1)
    {
        A.diag(-1)  = D1;
        A.diag(1)   = conj(D1);
    };

    return linsolve_chol_tridiag_fac(A, S, D, upper, opts);
};

linsolve_obj matcl::linsolve_chol_tridiag_fac(const Matrix& A, const Matrix& S, const Matrix& D, bool upper,
                                              const options& opts)
{
    if (upper == true && matcl::is_triu(S) == false)
        throw error::invalid_cholecky_factors();

    if (upper == false && matcl::is_tril(S) == false)
        throw error::invalid_cholecky_factors();

    if (D.rows() != D.cols())
        throw error::invalid_cholecky_factors();

    if (S.rows() != S.cols())
        throw error::square_matrix_required(S.rows(), S.cols());

    if (S.rows() != D.rows())
        throw error::invalid_cholecky_factors();

    if (is_diag(D) == false)
        throw error::invalid_cholecky_factors();

    if (upper == false && A.rows() != S.rows())
        throw error::invalid_cholecky_factors();
    else if (upper == true && A.rows() != S.cols())
        throw error::invalid_cholecky_factors();

    if (get_ld(A,1) > 1 || get_ud(A,1) > 1)
    {
        Integer ld  = get_ld(A,-1);
        Integer ud  = get_ud(A,-1);
        throw error::tridiagonal_matrix_required(A.rows(), A.cols(), ld, ud);
    };

    if (S.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_chol_tridiag");
    if (D.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_chol_tridiag");
    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_chol_tridiag");

    Integer N = S.rows();

    value_code vc1  = S.get_value_code();
    value_code vc2  = D.get_value_code();    
    value_code vc   = matrix_traits::unify_value_types(vc1, A.get_value_code());
    vc              = matrix_traits::unify_value_types(vc, value_code::v_float);
    vc              = matrix_traits::unify_value_types(vc1, vc2);

    if (N == 0)
    {
        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, S.get_type() )));
    };

    Matrix dum = zeros(0,0,vc);
    linsolve_obj ret;
    details::linsolve_tridiag_fac_vis::make<const Matrix&>(dum, ret, A, D, S, upper, opts);
    return ret;
};

struct linsolve_chol_impl
{
    static linsolve_obj eval(const Matrix& A, bool pivot, const options& opts)
    {
        if (A.rows() != A.cols())
            throw error::square_matrix_required(A.rows(), A.cols());

        Integer N = A.rows();

        if (N == 0)
        {
            value_code vc   = A.get_value_code();
            vc              = matrix_traits::unify_value_types(vc, value_code::v_float);

            using data_ptr = linsolve_obj::linsolve_data_ptr;
            return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
        };

        bool balance    = false;
        Real tol_sig    = 0.0;

        if (pivot == false && opts.get_option<bool>(opt::linsolve::do_balancing()))
        {
            balance     = true;
            tol_sig     = 0.0;
        }
        else if ((pivot == true) && opts.get_option<bool>(opt::linsolve::do_balancing_rr()))
        {
            balance     = true;
            tol_sig     = opts.get_option<Real>(opt::linsolve::tol_sing());
            if (tol_sig < 0.0)
            {
                tol_sig = - tol_sig;
            }
            else if (tol_sig > 0.0)
            {
                tol_sig = norm(A, -2.0) * pow(constants::eps(A.get_value_code()), tol_sig);
            };
        }

        if (pivot == true)
        {
            Matrix T;
            permvec p;
            Integer rank;
            bool upper  = false;

            Matrix B, D;
            if (balance == true)
            {
                tie(B,D)        = balance_posdef2(A, true, tol_sig);
                tie(T,p,rank)   = chol_rr(B, upper, -1.0);
            }
            else
            {
                tie(T,p,rank)   = chol_rr(A, upper, -1.0);
            };

            Real tol            = opts.get_option<Real>(opt::linsolve::tol_sing());

            if (rank < A.rows() && tol == 0.0 || rank == 0)
                throw error::error_singular();

            if (balance == true)
            {
                linsolve_obj lo_B   = linsolve_chol(B, T, p, upper, opts);
                return linsolve_balanced_sym(A, D, lo_B);
            }
            else
            {
                return linsolve_chol(A, T, p, upper, opts);
            };
        };

        Integer ld = matcl::get_ld(A, 0);
        if (ld == 0)
            return linsolve_triang(A, opts);

        Integer ud = matcl::get_ud(A, 0);
        if (ud == 0)
            return linsolve_triang(A, opts);

        if (ld <= 1 && ud <= 1)
        {
            Matrix D0 = A.diag(0);
            Matrix D1 = A.diag(-1);
            return linsolve_chol_tridiag(D0, D1, opts);
        };

        Matrix T;
        permvec p;
        bool upper  = false;        

        if (balance == true)
        {            
            Matrix B, D;

            tie(B,D) = balance_posdef2(A, true, tol_sig);
            tie(T,p) = chol(B, upper, opts);

            linsolve_obj lo_B   = linsolve_chol(B, T,p,upper, opts);
            return linsolve_balanced_sym(A, D,lo_B);
        }
        else
        {
            tie(T,p)    = chol(A, upper, opts);
            return linsolve_chol(A, T, p, upper, opts);
        };                
    }
};

linsolve_obj matcl::linsolve_chol(const Matrix& A, bool pivot, const options& opts)
{
    return linsolve_chol_impl::eval(A,pivot, opts);
}

linsolve_obj matcl::linsolve_chol(Matrix&& A, bool pivot, const options& opts)
{
    return linsolve_chol_impl::eval(std::move(A),pivot,opts);
}

linsolve_obj matcl::linsolve_chol_tridiag(const Matrix& D0_0, const Matrix& D1_0, const options& opts)
{
    Matrix D0(D0_0);
    Matrix D1(D1_0);

    linsolve_obj ret;
    linsolve_tridiag_impl(ret, D0, D1, opts);

    return ret;
};

linsolve_obj matcl::linsolve_chol_tridiag(const Matrix& D0_0, Matrix&& D1_0, const options& opts)
{
    Matrix D0(D0_0);
    Matrix D1(std::move(D1_0));

    linsolve_obj ret;
    linsolve_tridiag_impl(ret, D0, D1, opts);

    return ret;
};

linsolve_obj matcl::linsolve_chol_tridiag(Matrix&& D0_0, const Matrix& D1_0, const options& opts)
{
    Matrix D0(std::move(D0_0));
    Matrix D1(D1_0);

    linsolve_obj ret;
    linsolve_tridiag_impl(ret, D0, D1, opts);

    return ret;
};

linsolve_obj matcl::linsolve_chol_tridiag(Matrix&& D0_0, Matrix&& D1_0, const options& opts)
{
    Matrix D0(std::move(D0_0));
    Matrix D1(std::move(D1_0));

    linsolve_obj ret;
    linsolve_tridiag_impl(ret, D0, D1, opts);

    return ret;
};

};
