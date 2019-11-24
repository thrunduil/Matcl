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

#include "matcl-linalg/decompositions/ldl.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-internals/func/test_inf_nan.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/base/pv_constructor.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/linear_eq/linsolve_objects_decomp.h"
#include "matcl-linalg/decompositions/balancing.h"
#include "matcl-linalg/linear_eq/linsolve.h"

namespace matcl
{

namespace details
{

//--------------------------------------------------------------------------------
//                                  DENSE
//--------------------------------------------------------------------------------
template<class V, class S> struct ldl_str {};

template<class V>
struct ldl_str<V, struct_dense>
{
    using Mat   = raw::Matrix<V, struct_dense>;
    using Mat_b = raw::Matrix<V,struct_banded>;
    
    static void eval(const Mat& A, Matrix& ret_L, Matrix& ret_D, permvec& ret_P, 
                     bool hermitian, bool upper)
    {
        Integer N       = A.cols();

        Integer* ptr_PIV;
        matcl::Matrix PIV   = matcl::make_integer_dense_noinit(N,1,ptr_PIV);
        
        Mat Ac(A.get_type());
        Mat work(A.get_type(), 1, 1);

        if (A.get_refstr()->is_unique() == true)
            Ac.assign_to_fresh(A);
        else
            Ac.assign_to_fresh(A.copy());
        
        Ac.get_struct().reset();

        Integer info    = 0;
                
        const char uplo = upper ? 'U' : 'L';

        if(hermitian)
        {
            lapack::hetrf(&uplo, N, lap(Ac.ptr()), Ac.ld(), lap(ptr_PIV), lap(work.ptr()), 
                          -1, lap(&info));
            
            Integer lwork = (Integer) real(work.ptr()[0]);
            work.reset_unique(lwork, 1);
                        
            lapack::hetrf(&uplo, N, lap(Ac.ptr()), Ac.ld(), lap(ptr_PIV), lap(work.ptr()), 
                          lwork, lap(&info));
        }
        else
        {
            lapack::sytrf(&uplo, N, lap(Ac.ptr()), Ac.ld(), lap(ptr_PIV), lap(work.ptr()),
                          -1, lap(&info));

            Integer lwork = (Integer) real(work.ptr()[0]);
            work.reset_unique(lwork, 1);
            
            lapack::sytrf(&uplo, N, lap(Ac.ptr()), Ac.ld(), lap(ptr_PIV), lap(work.ptr()), 
                          lwork, lap(&info));
        }

        if(info < 0)
            throw error::error_ldl();

        Matrix P            = irange(1,N);
        Mat_b D             = Mat_b(A.get_type(), N, N, -1, 1);        

        if (upper == true)
            post_ldl_upper(Ac, D, P, ptr_PIV, hermitian);
        else
            post_ldl_lower(Ac, D, P, ptr_PIV, hermitian);

        ret_P           = details::pv_constructor::make(P);        
        ret_L           = Matrix(Ac,true);
        ret_D           = Matrix(D,true);

        if (hermitian)
        {
            if (is_real_matrix(ret_D))
                ret_D.add_struct(predefined_struct_type::sym);
            else
                ret_D.add_struct(predefined_struct_type::her);
        }
        else
        {
            ret_D.add_struct(predefined_struct_type::sym);
        }

        ret_D.add_struct(predefined_struct_type::qtril);
        ret_D.add_struct(predefined_struct_type::qtriu);

        if (upper == false)
            ret_L.add_struct(predefined_struct_type::tril);
        else
            ret_L.add_struct(predefined_struct_type::triu);
    }

    static void post_ldl_lower(Mat& Ac, Mat_b& D, Matrix& P, const Integer* ptr_PIV, 
                         bool hermitian)
    {
        Integer D_ld        = D.ld();
        Integer N           = Ac.cols();
        V* ptr_D_0          = D.rep_ptr() + D.first_elem_diag(0);
        V* ptr_D_l1         = D.rep_ptr() + D.first_elem_diag(-1);
        V* ptr_D_u1         = D.rep_ptr() + D.first_elem_diag(1);

        V* ptr_A            = Ac.ptr();
        Integer A_ld        = Ac.ld();
        
        Integer *ptr_P      = P.get_array_unique<Integer>();        

        for (Integer i = 0; i < N; i++)
        {                
            Integer ip      = ptr_PIV[i];
            bool step_2     = false;

            Integer p1, p2;

            if (ip > 0)
            {
                p1      = i;
                p2      = ip - 1;                    
            }
            else
            {
                p1      = i + 1;
                p2      = -ip - 1;
                step_2  = true;                    
            };

            std::swap(ptr_P[p1], ptr_P[p2]);

            //apply delayed permutations
            if (p1 != p2)
            {
                V* ptr          = Ac.ptr();
                for (Integer j = 0; j < i; ++j)
                {
                    std::swap(ptr[p1], ptr[p2]);
                    ptr         += A_ld;
                };
            };

            //clear upper-triangular part
            for (Integer j = 0; j < i; ++j)
                ptr_A[j]        = V(0.0);

            ptr_D_0[0]          = ptr_A[i];
            ptr_A[i]            = V(1.0);

            if (step_2)
            {
                V subdiag       = ptr_A[i + 1];

                ptr_D_l1[0]     = subdiag;
                if (hermitian)
                    ptr_D_u1[0] = conj(subdiag);
                else
                    ptr_D_u1[0] = subdiag;

                ptr_A[i + 1]    = V(0);

                ptr_A           += A_ld;

                ptr_D_0         += D_ld;
                ptr_D_l1        += D_ld;
                ptr_D_u1        += D_ld;

                //clear upper-triangular part
                for (Integer j = 0; j < i+1; ++j)
                    ptr_A[j]    = V(0.0);
            
                ptr_D_0[0]      = ptr_A[i+1];

                if (i < N - 2)
                {
                    ptr_D_l1[0] = V(0.0);
                    ptr_D_u1[0] = V(0.0);
                };

                ptr_A[i+1]      = V(1.0);
                ++i;                
            }
            else if (i < N - 1)
            {
                ptr_D_l1[0]     = V(0.0);
                ptr_D_u1[0]     = V(0.0);
            };
                
            ptr_A               += A_ld;
            ptr_D_0             += D_ld;
            ptr_D_l1            += D_ld;
            ptr_D_u1            += D_ld;
        }
    };
    static void post_ldl_upper(Mat& Ac, Mat_b& D, Matrix& P, const Integer* ptr_PIV, 
                         bool hermitian)
    {
        Integer D_ld        = D.ld();
        Integer N           = Ac.cols();
        V* ptr_D_0          = D.rep_ptr() + (N-1) * D_ld + D.first_elem_diag(0);
        V* ptr_D_l1         = D.rep_ptr() + (N-2) * D_ld + D.first_elem_diag(-1);
        V* ptr_D_u1         = D.rep_ptr() + (N-2) * D_ld + D.first_elem_diag(1);

        Integer A_ld        = Ac.ld();
        V* ptr_A            = Ac.ptr() + (N-1) * A_ld;
        
        Integer *ptr_P      = P.get_array_unique<Integer>();        

        for (Integer i = N-1; i >= 0; i--)
        {                
            Integer ip      = ptr_PIV[i];
            bool step_2     = false;

            Integer p1, p2;

            if (ip > 0)
            {
                p1      = i;
                p2      = ip - 1;                    
            }
            else
            {
                p1      = i - 1;
                p2      = -ip - 1;
                step_2  = true;                    
            };

            std::swap(ptr_P[p1], ptr_P[p2]);

            //apply delayed permutations
            if (p1 != p2)
            {
                V* ptr          = Ac.ptr() + (N-1)*A_ld;
                for (Integer j = i+1; j < N; ++j)
                {
                    std::swap(ptr[p1], ptr[p2]);
                    ptr         -= A_ld;
                };
            };

            //clear upper-triangular part
            for (Integer j = i+1; j < N; ++j)
                ptr_A[j]        = V(0.0);

            ptr_D_0[0]          = ptr_A[i];
            ptr_A[i]            = V(1.0);

            if (step_2)
            {
                V subdiag       = ptr_A[i - 1];

                ptr_D_u1[0]     = subdiag;
                if (hermitian)
                    ptr_D_l1[0] = conj(subdiag);
                else
                    ptr_D_l1[0] = subdiag;

                ptr_A[i - 1]    = V(0);

                ptr_A           -= A_ld;

                ptr_D_0         -= D_ld;
                ptr_D_l1        -= D_ld;
                ptr_D_u1        -= D_ld;

                //clear upper-triangular part
                for (Integer j = i; j < N; ++j)
                    ptr_A[j]    = V(0.0);
            
                ptr_D_0[0]      = ptr_A[i-1];

                if (i > 1)
                {
                    ptr_D_l1[0] = V(0.0);
                    ptr_D_u1[0] = V(0.0);
                };

                ptr_A[i-1]      = V(1.0);
                --i;                
            }
            else if (i > 0)
            {
                ptr_D_l1[0]     = V(0.0);
                ptr_D_u1[0]     = V(0.0);
            };
                
            ptr_A               -= A_ld;
            ptr_D_0             -= D_ld;
            ptr_D_l1            -= D_ld;
            ptr_D_u1            -= D_ld;
        }
    };

};

//--------------------------------------------------------------------------------
//                                  BAND
//--------------------------------------------------------------------------------
template<class V>
struct ldl_str<V, struct_banded>
{
    using Mat   = raw::Matrix<V, struct_banded>;
    using Mat_d = raw::Matrix<V, struct_dense>;
    
    static void eval(const Mat& A, Matrix& ret_L, Matrix& ret_D, permvec& ret_P, 
                     bool hermitian, bool upper)
    {
        using MC = raw::Matrix<V,struct_dense>;
        MC AC = raw::converter<MC,Mat>::eval(A);
        return ldl_str<V,struct_dense>::eval(AC,ret_L,ret_D,ret_P,hermitian,upper);
    };
};

//--------------------------------------------------------------------------------
//                                  SPARSE
//--------------------------------------------------------------------------------
template<class V>
struct ldl_str<V, struct_sparse>
{
    using Mat   = raw::Matrix<V, struct_sparse>;
    
    static void eval(const Mat& A, Matrix& ret_L, Matrix& ret_D, permvec& ret_P, 
                     bool hermitian, bool upper)
    {
        using MC = raw::Matrix<V,struct_dense>;
        MC AC = raw::converter<MC,Mat>::eval(A);
        return ldl_str<V,struct_dense>::eval(AC,ret_L,ret_D,ret_P,hermitian,upper);
    };
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct ldl_impl
{
    using M = raw::Matrix<V, S>;
    static void eval(const M& A, Matrix& ret_L, Matrix& ret_D, permvec& ret_P,
                     bool hermitian, bool upper)
    {
        return ldl_str<V,S>::eval(A, ret_L, ret_D, ret_P, hermitian, upper);
    }
};
template<class S>
struct ldl_impl<Integer, S>
{
    using M = raw::Matrix<Integer, S>;
    static void eval(const M& A, Matrix& ret_L, Matrix& ret_D, permvec& ret_P, 
                     bool hermitian, bool upper)
    {
        using MC = raw::Matrix<Real,struct_dense>;
        MC AC = raw::converter<MC,M>::eval(A);
        ldl_impl<Real, struct_dense>::eval(AC, ret_L, ret_D, ret_P, hermitian, upper);
    }
};
template<class S>
struct ldl_impl<Object, S>
{
    using M = raw::Matrix<Object, S>;
    static void eval(const M& , Matrix& , Matrix&, permvec&, bool, bool )
    {
        throw error::object_value_type_not_allowed("ldl");
    }
};

struct unary_visitor_ldl : public extract_type_switch<void, unary_visitor_ldl,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, ldl_return_type& ret, bool hermitian, bool upper)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        
        bool isv    = mat.all_finite();
        
        if (isv == false)
        {
            Matrix mat_nan  = details::make_nan_matrix<V>(mat.rows(), mat.cols());
            permvec p       = permvec::identity(mat.rows());
            ret             = ldl_return_type(mat_nan, mat_nan, p);
            return;
        };

        Matrix l,d;
        permvec p;
        details::ldl_impl<V,S>::eval(mat, l, d, p, hermitian, upper);               

        ret = ldl_return_type(l, d ,p) ;
        return;
    }

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, ldl_return_type& ret, bool, bool)
    {
        bool isv = is_finite(v);

        if (isv == false)
        {
            Matrix nan_val = details::make_nan_matrix<T>(1, 1);
            permvec p   = permvec::identity(1);
            ret         = ldl_return_type(nan_val, nan_val, p);
            return;
        }

        using VR = typename md::real_type_int_real<T>::type;

        permvec p   = permvec::identity(1);
        ret         = ldl_return_type(VR(1.0), handle, p); // (l d p)
        return;
    }

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, ldl_return_type&, bool, bool)
    {
        throw error::object_value_type_not_allowed("ldl");
    };

    static void eval_scalar(const Matrix&, const Object&, ldl_return_type&, bool, bool)
    {
        throw error::object_value_type_not_allowed("ldl");
    };
};

};

static void ldl_impl(ldl_return_type& ret, const Matrix& A, bool upper)
{
    if (!A.is_square())
        throw error::error_size_ldl(A.rows(), A.cols());

    if (A.rows() <= 1 || A.structural_nnz() == 0 || is_diag(A))
    {
        value_code v    = A.get_value_code();
        v               = matrix_traits::real_value_type(v);
        v               = matrix_traits::unify_value_types(v,value_code::v_float);
        
        Matrix I        = speye(A.rows(), A.rows(), v);
        permvec p       = permvec::identity(A.rows());

        value_code vd   = matrix_traits::unify_value_types(v, A.get_value_code());
        Matrix Ad       = bdiag(get_diag(A));
        Ad              = matcl::convert_value(Ad, vd);

        ret         = ldl_return_type(I, Ad, p);
        return;
    }

    return details::unary_visitor_ldl::make<const Matrix&>(A, ret, false, upper);
};

static void ldl_herm_impl(ldl_return_type& ret, const Matrix& A, bool upper)
{
    if (!A.is_square())
        throw error::error_size_ldl(A.rows(), A.cols());

    if (A.rows() <= 1 || A.structural_nnz() == 0 || is_diag(A))
    {
        value_code v    = A.get_value_code();
        v               = matrix_traits::real_value_type(v);
        v               = matrix_traits::unify_value_types(v,value_code::v_float);

        Matrix I    = speye(A.rows(),A.rows(),v);
        permvec p   = permvec::identity(A.rows());

        value_code vd   = matrix_traits::unify_value_types(v, A.get_value_code());
        Matrix Ad       = bdiag(get_diag(A));
        Ad              = matcl::convert_value(Ad, vd);

        ret         = ldl_return_type(I, Ad, p);
        return;
    }

    return details::unary_visitor_ldl::make<const Matrix&>(A, ret, true, upper);
};

ldl_return_type matcl::ldl(const Matrix& A0, bool upper)
{
    //increase refcount
    Matrix A(A0);

    ldl_return_type ret;
    ldl_impl(ret, A, upper);
    return ret;
};

ldl_return_type matcl::ldl(Matrix&& A0, bool upper)
{
    Matrix A(std::move(A0));

    ldl_return_type ret;
    ldl_impl(ret, A, upper);
    return ret;
};

ldl_return_type matcl::ldl_herm(const Matrix& A0, bool upper)
{
    //increase refcount
    Matrix A(A0);

    ldl_return_type ret;
    ldl_herm_impl(ret, A, upper);
    return ret;
};
ldl_return_type matcl::ldl_herm(Matrix&& A0, bool upper)
{
    Matrix A(std::move(A0));

    ldl_return_type ret;
    ldl_herm_impl(ret, A, upper);
    return ret;
};

linsolve_obj matcl::linsolve_ldl(const Matrix& A, const Matrix& L, const Matrix& D, 
                                 const permvec& p,bool sym, const options& opts)
{
    if (A.rows() != L.rows())
        throw error::invalid_ldl_factors();
    if (A.cols() != L.rows())
        throw error::invalid_ldl_factors();
    if (L.cols() != D.rows())
        throw error::invalid_ldl_factors();
    if (D.rows() != D.cols())
        throw error::invalid_ldl_factors();
    if (L.rows() != p.length())
        throw error::invalid_ldl_factors();

    if (matcl::is_triu(L) == false && matcl::is_tril(L) == false)
        throw error::invalid_ldl_factors();

    Integer ld  = matcl::get_ld(D, 1);
    Integer ud  = matcl::get_ud(D, 1);

    if (ld > 1 || ud > 1)
        throw error::invalid_ldl_factors();

    if (L.rows() != L.cols())
        throw error::square_matrix_required(L.rows(), L.cols());

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_ldl");
    if (L.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_ldl");
    if (D.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_ldl");

    Integer N   = L.rows();

    if (N == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(L.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(D.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, D.get_type())));
    };

    bool isv            = L.all_finite() && D.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(L.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(D.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, D.get_type() )));
    };

    using data_ptr  = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_ldl(A, L, D, p, sym, opts)));
};

struct linsolve_ldl_impl
{
    static linsolve_obj eval(const Matrix& A, bool sym, const options& opts)
    {
        if (A.rows() != A.cols())
            throw error::square_matrix_required(A.rows(), A.cols());

        if (sym == true && matrix_traits::is_float_real(A.get_value_code()))
            sym = false;

        bool upper  = false;

        bool balance    = opts.get_option<bool>(opt::linsolve::do_balancing());

        Matrix B, BD;

        if (balance)
        {
            Real min_diag   = 0.0;
            tie(B,BD)       = balance_sym2(A, true, min_diag);
        }
        else
        {
            B               = A;
        };

        Matrix L, D;
        permvec p;

        if (sym == true)
        {
            tie(L,D,p)  = ldl(B, upper);            
        }
        else
        {
            tie(L,D,p)  = ldl_herm(B, upper);
        };

        if (balance)
        {
            linsolve_obj lo_B = linsolve_ldl(B, L,D,p,sym, opts);
            return linsolve_balanced_sym(A, BD,lo_B);
        }
        else
        {
            return linsolve_ldl(A, L,D,p,sym, opts);
        };
    }
};

linsolve_obj matcl::linsolve_ldl(const Matrix& A, bool sym, const options& opts)
{
    return linsolve_ldl_impl::eval(A,sym, opts);
}
linsolve_obj matcl::linsolve_ldl(Matrix&& A, bool sym, const options& opts)
{
    return linsolve_ldl_impl::eval(std::move(A),sym, opts);
}

};