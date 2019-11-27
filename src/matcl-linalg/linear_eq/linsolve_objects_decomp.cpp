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

#include "linsolve_objects_decomp.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/decompositions/qr.h"
#include "linsolve_utils.inl"
#include "matcl-linalg/utils/linalg_utils.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/norms_error/norm.h"
#include "matcl-core/utils/workspace.h"

namespace matcl { namespace details
{

//--------------------------------------------------------------------------------
//                  linsolve_obj_lu_factors
//--------------------------------------------------------------------------------

linsolve_obj_lu_factors::linsolve_obj_lu_factors(const Matrix& A, const Matrix& L, const Matrix& U, 
        const permvec& p, const permvec& q, const options& opts)
    :linsolve_obj_seq_2()
{
    permvec L_id = permvec::identity(L.cols());
    permvec U_id = permvec::identity(U.rows());
    linsolve_obj_seq_2::reset(A, linsolve_triang(L,p,L_id, opts), linsolve_triang(U,U_id,q, opts));
};

linsolve_obj_lu_factors::~linsolve_obj_lu_factors()
{};

//--------------------------------------------------------------------------------
//                  linsolve_obj_chol
//--------------------------------------------------------------------------------
linsolve_obj_chol::linsolve_obj_chol(const Matrix& A, const Matrix& L, const permvec& p,
                                    bool upper, const options& opts)
    :linsolve_obj_seq_2()
{
    permvec id   = permvec::identity(L.cols());

    if (upper == false)
    {        
        linsolve_obj L1      = linsolve_triang(L,p,id, opts);
        linsolve_obj L2      = ctrans(L1);
        linsolve_obj_seq_2::reset(A, L1, L2);
    }
    else
    {
        linsolve_obj L2      = linsolve_triang(L,p,id, opts);
        linsolve_obj L1      = ctrans(L2);
        linsolve_obj_seq_2::reset(A, L1, L2);
    }
};

linsolve_obj_chol::~linsolve_obj_chol()
{};

Real linsolve_obj_chol::log_det() const
{
    return m_A1.log_det() * 2;
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_svd
//--------------------------------------------------------------------------------
linsolve_obj_svd::linsolve_obj_svd(const Matrix& A, const unitary_matrix& U, const Matrix& S, 
                                   const unitary_matrix& V, const options& opts)
    :linsolve_obj_seq_3(A, linsolve_unitary(U), linsolve_diag(S, opts), linsolve_unitary(ctrans(V)))
{};

linsolve_obj_svd::~linsolve_obj_svd()
{};

Real linsolve_obj_svd::normest_2() const
{
    return m_A2.normest(basic_vector_norm::norm_2);
};
Real linsolve_obj_svd::mat_normest_2() const
{
    return m_A2.mat_normest(basic_vector_norm::norm_2);
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_schur
//--------------------------------------------------------------------------------
linsolve_obj_schur::linsolve_obj_schur(const Matrix& A, const unitary_matrix& U, const Matrix& T,
                                       const options& opts)
    :linsolve_obj_seq_3(A, linsolve_unitary(U), linsolve_uhess(T, opts), linsolve_unitary(ctrans(U)))
{};

linsolve_obj_schur::~linsolve_obj_schur()
{};

Real linsolve_obj_schur::normest_2() const
{
    return m_A2.normest(basic_vector_norm::norm_2);
};
Real linsolve_obj_schur::mat_normest_2() const
{
    return m_A2.mat_normest(basic_vector_norm::norm_2);
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_bidiag
//--------------------------------------------------------------------------------
linsolve_obj_bidiag::linsolve_obj_bidiag(const Matrix& A, const unitary_matrix& U, const Matrix& R,
                                         const unitary_matrix& V, const options& opts)
    :linsolve_obj_seq_3(A, linsolve_unitary(U), linsolve_triang(R, opts), linsolve_unitary(ctrans(V)))
{};

linsolve_obj_bidiag::~linsolve_obj_bidiag()
{};

Real linsolve_obj_bidiag::normest_2() const
{
    return m_A2.normest(basic_vector_norm::norm_2);
};
Real linsolve_obj_bidiag::mat_normest_2() const
{
    return m_A2.mat_normest(basic_vector_norm::norm_2);
};
//--------------------------------------------------------------------------------
//                  linsolve_obj_rq
//--------------------------------------------------------------------------------
linsolve_obj_rq::linsolve_obj_rq(const Matrix& A, const Matrix& RL, const unitary_matrix& Q, 
                                 const options& opts)
    :linsolve_obj_seq_2(A, linsolve_triang(RL, opts), linsolve_unitary(Q))
{};
linsolve_obj_rq::~linsolve_obj_rq()
{};

Real linsolve_obj_rq::normest_2() const
{
    return m_A1.normest(basic_vector_norm::norm_2);
};
Real linsolve_obj_rq::mat_normest_2() const
{
    return m_A1.mat_normest(basic_vector_norm::norm_2);
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_qr
//--------------------------------------------------------------------------------
linsolve_obj_qr::linsolve_obj_qr(const Matrix& A, const unitary_matrix& Q, const Matrix& RL, 
                                 const permvec& p, const options& opts)
    :linsolve_obj_seq_2(A, linsolve_unitary(Q), linsolve_triang(RL, permvec::identity(RL.rows()), p, opts))
{};

linsolve_obj_qr::~linsolve_obj_qr()
{};

Real linsolve_obj_qr::normest_2() const
{
    return m_A2.normest(basic_vector_norm::norm_2);
};
Real linsolve_obj_qr::mat_normest_2() const
{
    return m_A2.mat_normest(basic_vector_norm::norm_2);
};
//--------------------------------------------------------------------------------
//                  linsolve_obj_hess
//--------------------------------------------------------------------------------
linsolve_obj_hess::linsolve_obj_hess(const Matrix& A, const unitary_matrix& U, const Matrix& H, 
                                     const options& opts)
    :linsolve_obj_seq_3(A, linsolve_unitary(U), linsolve_uhess(H, opts), linsolve_unitary(ctrans(U)))
{};

linsolve_obj_hess::~linsolve_obj_hess()
{}

Real linsolve_obj_hess::normest_2() const
{
    return m_A2.normest(basic_vector_norm::norm_2);
};
Real linsolve_obj_hess::mat_normest_2() const
{
    return m_A2.mat_normest(basic_vector_norm::norm_2);
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_ldl
//--------------------------------------------------------------------------------
linsolve_obj_ldl::linsolve_obj_ldl(const Matrix& A, const Matrix& L, const Matrix& D, 
                                   const permvec& p, bool sym, const options& opts)
    :m_sym(sym)
{
    bool is_real    = matrix_traits::is_float_real(this->get_value_code());
    permvec id      = permvec::identity(L.cols());

    if (sym == false || is_real == true)
    {
        linsolve_obj L1      = linsolve_triang(L,p,id, opts);
        linsolve_obj L2      = linsolve_diag_22(D, opts);
        linsolve_obj L3      = ctrans(L1);
        linsolve_obj_seq_3::reset(A, L1, L2, L3);
    }
    else
    {
        linsolve_obj L1      = linsolve_triang(L,p,id, opts);
        linsolve_obj L2      = linsolve_diag_22(D, opts);
        linsolve_obj L3      = trans(L1);
        linsolve_obj_seq_3::reset(A, L1, L2, L3);
    };
};

linsolve_obj_ldl::linsolve_obj_ldl(const Matrix& A, const linsolve_obj& L1, const linsolve_obj& D, 
                                         const linsolve_obj& L2, bool sym)
    :linsolve_obj_seq_3(A, L1, D, L2), m_sym(sym)
{};

linsolve_obj_ldl::~linsolve_obj_ldl()
{};

Real linsolve_obj_ldl::log_det() const
{
    Real ld = m_A1.log_det() * 2 + m_A2.log_det();
    return ld;
};

linsolve_obj_ldl::data_ptr 
linsolve_obj_ldl::convert(value_code vc) const
{
    Matrix Ac = details::convert_value(m_base_mat, vc);
    return data_ptr(new linsolve_obj_ldl(Ac, m_A1.convert(vc), m_A2.convert(vc), m_A3.convert(vc), m_sym));
};

//--------------------------------------------------------------------------------
//                  helpers
//--------------------------------------------------------------------------------

template<class Val, class Fac>
struct linsolve_obj_convert : public extract_type_switch<void, linsolve_obj_convert<Val, Fac>,true>
{
    using umatrix   = Fac;
    using data_ptr  = typename umatrix::data_ptr;

    template<class T>
    static void eval(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        using VT    = typename T::value_type;
        using VM    = typename md::unify_types<VT, Float>::type;

        ret = data.convert_impl<VM>();
    };
    template<class T>
    static void eval_scalar(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        using VM    = typename md::unify_types<T, Float>::type;
        ret = data.convert_impl<VM>();
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("linsolve_obj::convert");
    };
    static void eval_scalar(const Matrix&, const Object&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("linsolve_obj::convert");
    };
};

template<class Val, class Fac>
struct linsolve_obj_solve : public extract_type_switch<void, linsolve_obj_solve<Val, Fac>,true>
{
    using umatrix   = Fac;

    template<class T>
    static void eval(const Matrix&, const T&, const umatrix& data, Matrix& ret, Matrix& X, trans_type tA)
    {
        using VM    = typename T::value_type;
        using VR    = typename md::unify_types<VM, Float>::type;

        return data.solve_impl<VR>(ret, X, tA);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, const umatrix& data, Matrix& ret, Matrix& X, trans_type tA)
    {
        using VM    = T;
        using VR    = typename md::unify_types<VM, Float>::type;

        return data.solve_impl<VR>(ret, X, tA);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const umatrix&, Matrix&, Matrix&, trans_type)
    {
        throw error::object_value_type_not_allowed("linsolve_obj::solve");
    };
    static void eval_scalar(const Matrix&, const Object&, const umatrix&, Matrix&, Matrix&, trans_type)
    {
        throw error::object_value_type_not_allowed("linsolve_obj::solve");
    };
};

template<class Val, class Fac>
struct linsolve_obj_solve_rev : public extract_type_switch<void, linsolve_obj_solve_rev<Val, Fac>,true>
{
    using umatrix   = Fac;

    template<class T>
    static void eval(const Matrix&, const T&, const umatrix& data, Matrix& ret, Matrix& X, trans_type tA)
    {
        using VM    = typename T::value_type;
        using VR    = typename md::unify_types<VM, Float>::type;

        return data.solve_rev_impl<VR>(ret, X, tA);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, const umatrix& data, Matrix& ret, Matrix& X, trans_type tA)
    {
        using VM    = T;
        using VR    = typename md::unify_types<VM, Float>::type;

        return data.solve_rev_impl<VR>(ret, X, tA);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const umatrix&, Matrix&, Matrix&, trans_type)
    {
        throw error::object_value_type_not_allowed("linsolve_obj::solve_rev");
    };
    static void eval_scalar(const Matrix&, const Object&, const umatrix&, Matrix&, Matrix&, trans_type)
    {
        throw error::object_value_type_not_allowed("linsolve_obj::solve_rev");
    };
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_lu_dense
//--------------------------------------------------------------------------------
template<class V>
linsolve_obj_lu_dense<V>::linsolve_obj_lu_dense(const Mat& A, const Mat& A_dec, const Mat_I& ipiv, 
                                                const options& opts)
    : linsolve_obj_base(A.rows(), A.cols(), matrix_traits::value_code<V>::value, 
                           ti::convert_ti_object<V>(A.get_type()), false)
    , m_A(A), m_A_decomp(A_dec), m_piv(ipiv)
{
    correct_singular(m_A_decomp, opts);
};

template<class V>
linsolve_obj_lu_dense<V>::linsolve_obj_lu_dense(const Mat& A, const Mat& A_dec, const Mat_I& ipiv, 
                                                bool modif, from_inv)
    : linsolve_obj_base(A.rows(), A.cols(), matrix_traits::value_code<V>::value, 
                           ti::convert_ti_object<V>(A.get_type()), modif)
    , m_A(A), m_A_decomp(A_dec), m_piv(ipiv)
{};

template<class V>
linsolve_obj_lu_dense<V>::~linsolve_obj_lu_dense()
{};

template<class V>
void linsolve_obj_lu_dense<V>::correct_singular(Mat& A, const options& opts)
{
    using VR    = typename md::real_type<V>::type;

    Real tol    = opts.get_option<Real>(opt::linsolve::tol_sing());
    Integer N   = A.rows();
    Integer ld  = A.ld();
    V* ptr      = A.ptr();

    // calculate small pivot tolerance
    if (tol >= 0.0)
    {
        VR val_max  = VR(0.0);
        VR val_min  = VR(1.0);        

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i <= j; ++i)
            {
                VR val  = abs(ptr[i]);
                val_max = std::max(val, val_max);            
            };

            val_min = std::min(abs(ptr[j]), val_min);

            ptr     += ld;
        };

        if (val_max == VR(0.0))
            throw error::error_singular();    

        if (tol == 0.0 && val_min == VR(0.0))
            throw error::error_singular();

        if (tol > 0.0)
            tol         = pow(constants::eps<VR>(), tol) * val_max;
    }
    else
    {
        tol         = -tol;
    };
   
    if (tol == 0.0)
        return;

    VR tol_t        = (VR)tol;
    bool need_modif = false;

    Integer j;

    for (j = 0; j < N; ++j)
    {
        VR val  = abs(ptr[0]);

        if (abs(val) < tol_t)
        {
            need_modif = true;
            break;
        };

        ptr     += ld + 1;
    };

    if (need_modif == false)
        return;

    A.assign_to_fresh(A.make_unique());

    ptr         = A.ptr();
    ld          = A.ld();

    for (; j < N; ++j)
    {
        VR val  = abs(ptr[0]);

        if (val < tol_t)
            ptr[0]  = tol_t * (val == VR(0.0) ? VR(1.0) : sign(ptr[0]));

        ptr     += ld + 1;
    };

    m_modified  = true;
};

template<class V>
typename linsolve_obj_lu_dense<V>::data_ptr 
linsolve_obj_lu_dense<V>::convert(value_code new_val_code) const
{
    using Fac       = linsolve_obj_lu_dense;
    value_code vc   = matrix_traits::unify_value_types(new_val_code, value_code::v_float);
    Matrix dum      = zeros(0,0, vc);

    data_ptr ret;
    linsolve_obj_convert<V,Fac>::make<const Matrix&>(dum,*this, ret);
    return ret;
};

template<class V>
template<class T>
typename linsolve_obj_lu_dense<V>::data_ptr 
linsolve_obj_lu_dense<V>::convert_impl() const
{
    using Mat_C     = raw::Matrix<T,struct_dense>;

    Mat_C Ac        = raw::converter<Mat_C, Mat>::eval(m_A);
    Mat_C Ac_dec    = raw::converter<Mat_C, Mat>::eval(m_A_decomp);

    return data_ptr(new linsolve_obj_lu_dense<T>(Ac, Ac_dec, m_piv, m_modified, from_inv()));
};

template<class V>
matcl::Matrix linsolve_obj_lu_dense<V>::inv() const
{
    Integer N   = this->rows();
    Mat X       = m_A_decomp.copy();

    V work_query;
    Integer info;
    lapack::getri( N, lap(X.ptr()), X.ld(), lap(m_piv.ptr()), lap(&work_query), -1, info);

    Integer LWORK       = (Integer)real(work_query);

    using VTR_pod       = matcl::pod_type<V>;
    using workspace     = matcl::pod_workspace<VTR_pod>;
    workspace WORK      = workspace(LWORK);
    V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

    lapack::getri( N, lap(X.ptr()), X.ld(), lap(m_piv.ptr()), lap(ptr_WORK), -1, info);

    if (info > 0)
        throw error::error_singular();

    if (is_hermitian() == true)
    {
        X.add_struct(predefined_struct_type::her);

        if (is_posdef() == true)
            X.add_struct(posdef_flag());
    };

    return Matrix(X,true);
};

template<class V>
matcl::Matrix linsolve_obj_lu_dense<V>::base_matrix() const
{
    return Matrix(m_A,false);
};

template<class V>
matcl::Matrix linsolve_obj_lu_dense<V>::solve(const Matrix& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_lu_dense;

    Matrix X(X0);

    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_lu_dense<V>::solve(Matrix&& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_lu_dense;

    Matrix X(std::move(X0));

    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_lu_dense<V>::solve_rev(const Matrix& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_lu_dense;

    Matrix X(X0);

    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve_rev<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_lu_dense<V>::solve_rev(Matrix&& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_lu_dense;

    Matrix X(std::move(X0));

    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve_rev<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
template<class T>
void linsolve_obj_lu_dense<V>::solve_impl(Matrix& ret, Matrix& X, trans_type tA) const
{
    Integer N       = m_A_decomp.rows();
    Integer Nrhs    = X.cols();

    using M         = raw::Matrix<T,struct_dense>;

    const M& Ac     = raw::converter<M,Mat>::eval(m_A_decomp);
    M Bc            = X.impl_unique<M>();

    Integer info;

    // Solve the system A*X = B, overwriting B with X.
    lapack::getrs( get_trans_code(tA), N, Nrhs, lap(Ac.ptr()), Ac.ld(), lap(m_piv.ptr()), 
                lap(Bc.ptr()), Bc.ld(), &info);

    if (info != 0)
        throw error::error_general("invalid argument passed to getrs");

    Bc.get_struct().reset();
    ret = Matrix(Bc,true);
    return;
};

template<class V>
template<class T>
void linsolve_obj_lu_dense<V>::solve_rev_impl(Matrix& ret, Matrix& X, trans_type tA) const
{
    Integer N       = m_A_decomp.rows();
    Integer Mrhs    = X.rows();

    using M         = raw::Matrix<T,struct_dense>;

    const M& Ac     = raw::converter<M,Mat>::eval(m_A_decomp);
    M Bc            = X.impl_unique<M>();

    Integer info;

    // Solve the system X*A = B, overwriting B with X.
    lapack::getrs_rev( get_trans_code(tA), N, Mrhs, lap(Ac.ptr()), Ac.ld(), lap(m_piv.ptr()), 
                lap(Bc.ptr()), Bc.ld(), &info);

    if (info != 0)
        throw error::error_general("invalid argument passed to getrs_rev");

    Bc.get_struct().reset();
    ret = Matrix(Bc,true);
    return;
};

template<class V>
Real linsolve_obj_lu_dense<V>::mat_normest_1() const
{
    return matcl::norm(Matrix(m_A, false), 1.0);
};

template<class V>
Real linsolve_obj_lu_dense<V>::mat_normest_inf() const
{
    return norm(Matrix(m_A, false), constants::inf());
};

template<class V>
bool linsolve_obj_lu_dense<V>::is_hermitian() const
{
    bool is_real = details::is_complex<V>::value == false;
    return m_A.get_struct().is_hermitian(m_A.rows() == m_A.cols(), is_real);
};

template<class V>
bool linsolve_obj_lu_dense<V>::is_posdef() const
{
    return matcl::is_posdef(m_A.get_struct());
};

template<class V>
Real linsolve_obj_lu_dense<V>::log_det() const
{
    Integer M       = m_A_decomp.rows();
    Integer N       = m_A_decomp.cols();
    const V* ptr    = m_A_decomp.ptr();
    Integer ld      = m_A_decomp.ld();

    Real det        = 0.0;

    if (M != N)
        return -constants::inf();

    for (Integer i = 0; i < M; ++i)
    {
        det         += matcl::log(matcl::abs(ptr[0]));
        ptr         += ld + 1;
    };
    
    return det;
};

template<class V>
Matrix linsolve_obj_lu_dense<V>::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    return mmul(Matrix(m_A, false), X, t);
};

template<class V>
Matrix linsolve_obj_lu_dense<V>::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    return mmul(Matrix(m_A, false), std::move(X), t);
};

template<class V>
Matrix linsolve_obj_lu_dense<V>::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return mmul(X, Matrix(m_A,false), trans_type::no_trans,t);
};

template<class V>
Matrix linsolve_obj_lu_dense<V>::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return mmul(std::move(X), Matrix(m_A,false), trans_type::no_trans,t);
};

template class linsolve_obj_lu_dense<Real>;
template class linsolve_obj_lu_dense<Float>;
template class linsolve_obj_lu_dense<Complex>;
template class linsolve_obj_lu_dense<Float_complex>;

//--------------------------------------------------------------------------------
//                  linsolve_obj_lu_band
//--------------------------------------------------------------------------------
template<class V>
linsolve_obj_lu_band<V>::linsolve_obj_lu_band(const Mat& A, const Mat& A_dec, const Mat_I& ipiv,
                                                    Integer N, Integer ld, Integer ud, const options& opts)
    :linsolve_obj_base(N,N,matrix_traits::value_code<V>::value, ti::convert_ti_object<V>(A.get_type()), false)
    ,m_A(A), m_A_decomp(A_dec), m_piv(ipiv), m_N(N), m_ldiags(ld), m_udiags(ud)
{
    correct_singular(m_A_decomp, opts);
};
template<class V>
linsolve_obj_lu_band<V>::linsolve_obj_lu_band(const Mat& A, const Mat& A_dec, const Mat_I& ipiv,
                                                    Integer N, Integer ld, Integer ud, bool modif, from_inv)
    :linsolve_obj_base(N,N,matrix_traits::value_code<V>::value, ti::convert_ti_object<V>(A.get_type()), modif)
    ,m_A(A), m_A_decomp(A_dec), m_piv(ipiv), m_N(N), m_ldiags(ld), m_udiags(ud)
{};

template<class V>
linsolve_obj_lu_band<V>::~linsolve_obj_lu_band()
{};

template<class V>
void linsolve_obj_lu_band<V>::correct_singular(Mat& A, const options& opts)
{
    using VR        = typename md::real_type<V>::type;
    Integer ld      = A.last_diag();
    Integer A_ld    = A.ld();    

    Real tol        = opts.get_option<Real>(opt::linsolve::tol_sing());

    // calculate small pivot tolerance
    if (tol >= 0.0)
    {
        VR val_max  = VR(0.0);
        VR val_min  = VR(1.0);
        Integer s   = A.diag_length(0);
        V* ptr      = A.rep_ptr() + A.first_elem_diag(0);

        for (Integer i = 0; i < s; ++i)
        {
            VR val  = abs(ptr[0]);
            val_max = std::max(val, val_max);
            val_min = std::min(val, val_min);
            ptr     += A_ld;
        };

        for (Integer d = 1; d <= ld; ++d)
        {
            s       = A.diag_length(d);
            ptr     = A.rep_ptr() + A.first_elem_diag(d);

            for (Integer i = 0; i < s; ++i)
            {
                VR val  = abs(ptr[0]);
                val_max = std::max(val, val_max);            
                ptr     += A_ld;
            };
        };

        if (val_max == VR(0.0))
            throw error::error_singular();    

        if (tol == 0.0 && val_min == VR(0.0))
            throw error::error_singular();

        if (tol > 0.0)
            tol         = pow(constants::eps<VR>(), tol) * val_max;
    }
    else
    {
        tol         = -tol;
    }

    if (tol == 0.0)
        return;

    VR tol_t        = (VR)tol;
    bool need_modif = false;

    Integer i;
    {
        Integer s   = A.diag_length(0);
        V* ptr      = A.rep_ptr() + A.first_elem_diag(0);

        for (i = 0; i < s; ++i)
        {
            VR val  = abs(ptr[0]);
            
            if (abs(val) < tol_t)
            {
                need_modif = true;
                break;
            };

            ptr     += A_ld;
        };
    };

    if (need_modif == false)
        return;

    A.assign_to_fresh(A.make_unique());
    ld              = A.ld();

    {
        Integer s   = A.diag_length(0);
        V* ptr      = A.rep_ptr() + A.first_elem_diag(0);

        for (; i < s; ++i)
        {
            VR val  = abs(ptr[0]);
            
            if (val < tol_t)
                ptr[0]  = tol_t * (val == 0 ? VR(1.0) : sign(ptr[0]));

            ptr     += A_ld;
        };
    };

    m_modified  = true;
};

template<class V>
typename linsolve_obj_lu_band<V>::data_ptr 
linsolve_obj_lu_band<V>::convert(value_code new_val_code) const
{
    using Fac       = linsolve_obj_lu_band;
    value_code vc   = matrix_traits::unify_value_types(new_val_code, value_code::v_float);
    Matrix dum      = zeros(0,0, vc);

    data_ptr ret;
    linsolve_obj_convert<V,Fac>::make<const Matrix&>(dum,*this, ret);
    return ret;
};

template<class V>
template<class T>
typename linsolve_obj_lu_band<V>::data_ptr 
linsolve_obj_lu_band<V>::convert_impl() const
{
    using Mat_C     = raw::Matrix<T,struct_banded>;

    Mat_C Ac        = raw::converter<Mat_C, Mat>::eval(m_A);
    Mat_C Ac_dec    = raw::converter<Mat_C, Mat>::eval(m_A_decomp);

    return data_ptr(new linsolve_obj_lu_band<T>(Ac, Ac_dec, m_piv, m_N, m_ldiags, m_udiags, 
                                                m_modified, from_inv()));
};

template<class V>
bool linsolve_obj_lu_band<V>::is_hermitian() const
{
    bool is_real = details::is_complex<V>::value == false;
    return m_A.get_struct().is_hermitian(m_A.rows() == m_A.cols(), is_real);
};

template<class V>
bool linsolve_obj_lu_band<V>::is_posdef() const
{
    return matcl::is_posdef(m_A.get_struct());
};

template<class V>
matcl::Matrix linsolve_obj_lu_band<V>::inv() const
{
    Integer N   = this->rows();
    Matrix X    = eye(N, N, this->get_value_code()); 
    X           = this->solve(std::move(X), trans_type::no_trans);

    if (is_hermitian() == true)
    {
        X.add_struct(predefined_struct_type::her);

        if (is_posdef() == true)
            X.add_struct(posdef_flag());
    };

    return X;
};

template<class V>
matcl::Matrix linsolve_obj_lu_band<V>::base_matrix() const
{
    return Matrix(m_A,false);
};

template<class V>
matcl::Matrix linsolve_obj_lu_band<V>::solve(const Matrix& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_lu_band;

    Matrix X(X0);

    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_lu_band<V>::solve(Matrix&& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_lu_band;

    Matrix X(std::move(X0));

    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_lu_band<V>::solve_rev(const Matrix& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_lu_band;

    Matrix X(X0);

    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve_rev<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_lu_band<V>::solve_rev(Matrix&& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_lu_band;

    Matrix X(std::move(X0));

    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve_rev<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
template<class T>
void linsolve_obj_lu_band<V>::solve_impl(Matrix& ret, Matrix& X, trans_type tA) const
{
    Integer N       = m_N;
    Integer Nrhs    = X.cols();
    Integer ld      = m_ldiags;
    Integer ud      = m_udiags;

    using M         = raw::Matrix<T,struct_dense>;
    using MB        = raw::Matrix<T,struct_banded>;

    const MB& Ac    = raw::converter<MB,Mat>::eval(m_A_decomp);
    M Bc            = X.impl_unique<M>();

    Integer info;

    lapack::gbtrs( get_trans_code(tA), N, ld, ud, Nrhs, lap(Ac.rep_ptr()), Ac.ld(), lap(m_piv.ptr()),
                    lap(Bc.ptr()), Bc.ld(), &info );

    if (info != 0)
        throw error::error_general("invalid argument passed to gbtrs");

    Bc.get_struct().reset();
    ret = Matrix(Bc,true);
    return;
};

template<class V>
template<class T>
void linsolve_obj_lu_band<V>::solve_rev_impl(Matrix& ret, Matrix& X, trans_type tA) const
{
    Integer N       = m_N;
    Integer Mrhs    = X.rows();
    Integer ld      = m_ldiags;
    Integer ud      = m_udiags;

    using M         = raw::Matrix<T,struct_dense>;
    using MB        = raw::Matrix<T,struct_banded>;

    //solve
    const MB& Ac    = raw::converter<MB,Mat>::eval(m_A_decomp);
    M Bc            = X.impl_unique<M>();

    Integer info;

    lapack::gbtrs_rev( get_trans_code(tA), N, ld, ud, Mrhs, lap(Ac.rep_ptr()), Ac.ld(), lap(m_piv.ptr()),
                    lap(Bc.ptr()), Bc.ld(), &info );

    if (info != 0)
        throw error::error_general("invalid argument passed to gbtrs_rev");

    Bc.get_struct().reset();
    ret = Matrix(Bc,true);
    return;
};

template<class V>
Real linsolve_obj_lu_band<V>::log_det() const
{
    const V* ptr    = m_A_decomp.rep_ptr() + m_A_decomp.first_elem_diag(0);
    Integer ld      = m_A_decomp.ld();
    Integer M       = m_N;

    Real det        = 0.0;

    for (Integer i = 0; i < M; ++i)
    {
        det         += matcl::log(matcl::abs(ptr[0]));
        ptr         += ld;
    };
    
    return det;
};

template<class V>
Real linsolve_obj_lu_band<V>::mat_normest_1() const
{
    return matcl::norm(Matrix(m_A,false), 1.0);
};

template<class V>
Real linsolve_obj_lu_band<V>::mat_normest_inf() const
{
    return matcl::norm(Matrix(m_A,false), constants::inf());
};

template<class V>
Matrix linsolve_obj_lu_band<V>::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    return mmul(Matrix(m_A, false), X, t);
};

template<class V>
Matrix linsolve_obj_lu_band<V>::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    return mmul(Matrix(m_A, false), std::move(X), t);
};

template<class V>
Matrix linsolve_obj_lu_band<V>::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return mmul(X, Matrix(m_A,false), trans_type::no_trans,t);
};

template<class V>
Matrix linsolve_obj_lu_band<V>::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return mmul(std::move(X), Matrix(m_A,false), trans_type::no_trans,t);
};

template class linsolve_obj_lu_band<Real>;
template class linsolve_obj_lu_band<Float>;
template class linsolve_obj_lu_band<Complex>;
template class linsolve_obj_lu_band<Float_complex>;

//--------------------------------------------------------------------------------
//                  linsolve_obj_chol_tridiag_fac
//--------------------------------------------------------------------------------
template<class V>
linsolve_obj_chol_tridiag_fac<V>::linsolve_obj_chol_tridiag_fac(const Mat_B& A, 
                                const Mat_R& D0_i, const Mat& D1_i, bool modif, from_inv)
    :linsolve_obj_base(D0_i.length(),D0_i.length(),matrix_traits::value_code<V>::value,
                          ti::convert_ti_object<V>(D0_i.get_type()), modif)
    ,m_D0_i(D0_i), m_D1_i(D1_i), m_A(A)
{
    test_singular();
};

template<class V>
linsolve_obj_chol_tridiag_fac<V>::linsolve_obj_chol_tridiag_fac(const Mat_B& A, 
                                const Mat_R& D0_i, const Mat& D1_i, const options& opts)
    :linsolve_obj_base(D0_i.length(),D0_i.length(),matrix_traits::value_code<V>::value,
                          ti::convert_ti_object<V>(D0_i.get_type()), false)
    ,m_D0_i(D0_i.make_explicit()), m_D1_i(D1_i.make_explicit()), m_A(A)
{
    Real tol        = opts.get_option<Real>(opt::linsolve::tol_sing());

    // calculate small pivot tolerance
    if (tol > 0.0)
    {
        Real max_0  = norm_vec_all(Matrix(D0_i,false), basic_vector_norm::norm_inf);
        Real max_1  = norm_vec_all(Matrix(D1_i,false), basic_vector_norm::norm_inf);
        Real max    = std::max(max_0, max_1);
            
        if (max == 0.0)
            throw error::error_singular();

        tol         = pow(constants::eps<VR>(), tol) * max;
    }
    else if (tol < 0)
    {
        tol         = -tol;
    };

    VR tol_t        = (VR)tol;

    Integer N       = D0_i.length();
    const VR* ptr_D = D0_i.ptr();

    // perturb small pivots
    if (tol != 0.0)
    {
        bool need_modif = false;
        Integer j;

        for (j = 0; j < N; ++j)
        {
            V val   = ptr_D[j];

            if (abs(val) < tol_t)
            {
                need_modif = true;
                break;
            };                
        }

        if (need_modif == true)
        {
            m_D0_i.assign_to_fresh(m_D0_i.make_unique());
            VR* ptr_Dm  = m_D0_i.ptr();

            for (; j < N; ++j)
            {
                V val   = ptr_Dm[j];

                if (abs(val) < tol_t)
                    ptr_Dm[j]   = tol_t * (val == 0 ? VR(1.0) : sign(ptr_Dm[j]));
            }

            m_modified  = true;
        };        
    };

    test_singular();
};

template<class V>
void linsolve_obj_chol_tridiag_fac<V>::test_singular() const
{
    Integer N       = m_D0_i.rows();
    const VR* ptr_D = m_D0_i.ptr();

    for (Integer i = 0; i < N; ++i)
    {
        if (ptr_D[i] == VR(0.0))
            throw error::error_singular();
    };
};

template<class V>
linsolve_obj_chol_tridiag_fac<V>::~linsolve_obj_chol_tridiag_fac()
{};

template<class V>
matcl::Matrix linsolve_obj_chol_tridiag_fac<V>::inv() const
{
    Integer N   = this->rows();
    Matrix X    = eye(N, N, this->get_value_code()); 
    X           = this->solve(std::move(X), trans_type::no_trans);

    X.add_struct(predefined_struct_type::her);
    X.add_struct(posdef_flag());

    return X;
};

template<class V>
matcl::Matrix linsolve_obj_chol_tridiag_fac<V>::solve(const Matrix& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_chol_tridiag_fac;

    Matrix X(X0);

    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_chol_tridiag_fac<V>::solve(Matrix&& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_chol_tridiag_fac;

    Matrix X(std::move(X0));

    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_chol_tridiag_fac<V>::solve_rev(const Matrix& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_chol_tridiag_fac;

    Matrix X(X0);

    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve_rev<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_chol_tridiag_fac<V>::solve_rev(Matrix&& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_chol_tridiag_fac;

    Matrix X(std::move(X0));

    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve_rev<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
template<class T>
void linsolve_obj_chol_tridiag_fac<V>::solve_impl(Matrix& ret, Matrix& X, trans_type tA) const
{
    Integer N       = m_N;
    Integer Nrhs    = X.cols();

    using TR        = typename md::real_type<T>::type;
    using M         = raw::Matrix<T,struct_dense>;
    using MR        = raw::Matrix<TR,struct_dense>;

    const MR& D0c   = raw::converter<MR,Mat_R>::eval(m_D0_i);
    const M& D1c    = raw::converter<M,Mat>::eval(m_D1_i);
    M Bc            = X.impl_unique<M>();

    Integer info;

    const char* UPLO    = (tA == trans_type::no_trans)? "L" 
                        : ((tA == trans_type::conj_trans)? "L" : "U");
    lapack::pttrs(UPLO, N, Nrhs, lap(D0c.ptr()), lap(D1c.ptr()), lap(Bc.ptr()), Bc.ld(), info);

    if (info != 0)
        throw error::error_general("invalid argument passed to pttrs");

    Bc.get_struct().reset();
    ret = Matrix(Bc,true);
    return;
};

template<class V>
template<class T>
void linsolve_obj_chol_tridiag_fac<V>::solve_rev_impl(Matrix& ret, Matrix& X, trans_type tA) const
{
    Integer N       = m_N;
    Integer Mrhs    = X.rows();

    using TR        = typename md::real_type<T>::type;
    using M         = raw::Matrix<T,struct_dense>;
    using MR        = raw::Matrix<TR,struct_dense>;
    
    const MR& D0c   = raw::converter<MR,Mat_R>::eval(m_D0_i);
    const M& D1c    = raw::converter<M,Mat>::eval(m_D1_i);
    M Bc            = X.impl_unique<M>();

    Integer info;

    //solve
    const char* UPLO    = (tA == trans_type::no_trans)? "L" 
                        : ((tA == trans_type::conj_trans)? "L" : "U");

    lapack::pttrs_rev(UPLO, N, Mrhs, lap(D0c.ptr()), lap(D1c.ptr()), lap(Bc.ptr()), 
                      Bc.ld(), info);

    if (info != 0)
        throw error::error_general("invalid argument passed to pttrs_rev");

    Bc.get_struct().reset();
    ret = Matrix(Bc,true);
    return;
};

template<class V>
Real linsolve_obj_chol_tridiag_fac<V>::log_det() const
{
    Integer M       = m_N;

    const VR* ptr   = m_D0_i.ptr();
    Integer st      = m_D0_i.cols() == M ? m_D0_i.ld() : 1;

    Real det        = 0.0;

    for (Integer i = 0; i < M; ++i)
    {
        det         += matcl::log(matcl::abs(ptr[0]));
        ptr         += st;
    };
    
    return det;
};

template<class V>
Real linsolve_obj_chol_tridiag_fac<V>::normest_inf() const
{
    return normest_1();
};

template<class V>
Real linsolve_obj_chol_tridiag_fac<V>::normest_1() const
{
    using VL            = typename details::lapack_value_type<V>::type;
    Integer N           = this->rows();
    const VR* D         = m_D0_i.ptr();
    const V* E          = m_D1_i.ptr();

    using workspace     = matcl::pod_workspace<VR>;
    workspace WORK      = workspace(N);
    VR* ptr_WORK        = WORK.ptr();

    Integer info;
    VR norm;
    lapack::ptinorm<VL>(N, lap(D), lap(E), *lap(&norm), lap(ptr_WORK), info);

    return norm;
};

template<class V>
matcl::Matrix linsolve_obj_chol_tridiag_fac<V>::base_matrix() const
{
    return Matrix(m_A,false);
};
template<class V>
Real linsolve_obj_chol_tridiag_fac<V>::mat_normest_1() const
{
    return matcl::norm(Matrix(m_A,false), 1.0);
};

template<class V>
Real linsolve_obj_chol_tridiag_fac<V>::mat_normest_inf() const
{
    return matcl::norm(Matrix(m_A,false), constants::inf());
};

template<class V>
Matrix linsolve_obj_chol_tridiag_fac<V>::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return mmul(Matrix(m_A, false), X, t);
};

template<class V>
Matrix linsolve_obj_chol_tridiag_fac<V>::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return mmul(Matrix(m_A, false), std::move(X), t);
};

template<class V>
Matrix linsolve_obj_chol_tridiag_fac<V>::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return mmul(X, Matrix(m_A,false), trans_type::no_trans,t);
};

template<class V>
Matrix linsolve_obj_chol_tridiag_fac<V>::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return mmul(std::move(X), Matrix(m_A,false), trans_type::no_trans,t);
};

template<class V>
typename linsolve_obj_chol_tridiag_fac<V>::data_ptr 
linsolve_obj_chol_tridiag_fac<V>::convert(value_code vc0) const
{
    using Fac       = linsolve_obj_chol_tridiag_fac;
    value_code vc   = matrix_traits::unify_value_types(vc0, value_code::v_float);
    Matrix dum      = zeros(0,0, vc);

    data_ptr ret;
    linsolve_obj_convert<V,Fac>::make<const Matrix&>(dum,*this, ret);
    return ret;
};

template<class V>
template<class T>
typename linsolve_obj_chol_tridiag_fac<V>::data_ptr 
linsolve_obj_chol_tridiag_fac<V>::convert_impl() const
{
    using TR        = typename real_type<T>::type;
    using Mat_C     = raw::Matrix<T,struct_dense>;
    using Mat_BC    = raw::Matrix<T,struct_banded>;
    using Mat_RC    = raw::Matrix<TR,struct_dense>;

    Mat_BC Ac       = raw::converter<Mat_BC, Mat_B>::eval(m_A);
    Mat_RC D0_ic    = raw::converter<Mat_RC, Mat_R>::eval(m_D0_i);    
    Mat_C D1_ic     = raw::converter<Mat_C, Mat>::eval(m_D1_i);

    return data_ptr(new linsolve_obj_chol_tridiag_fac<T>(Ac, D0_ic, D1_ic, m_modified, from_inv()));
};

template class linsolve_obj_chol_tridiag_fac<Real>;
template class linsolve_obj_chol_tridiag_fac<Float>;
template class linsolve_obj_chol_tridiag_fac<Complex>;
template class linsolve_obj_chol_tridiag_fac<Float_complex>;

//--------------------------------------------------------------------------------
//                  linsolve_obj_tridiag_fac
//--------------------------------------------------------------------------------
template<class V>
linsolve_obj_tridiag_fac<V>::linsolve_obj_tridiag_fac(const Mat_B& A, const Mat& Dm1_i, const Mat& D0_i, 
            const Mat& Dp1_i, const Mat& Dp2_i, const Mat_I& piv, const options& opts)
    :linsolve_obj_base(D0_i.length(),D0_i.length(),matrix_traits::value_code<V>::value,
                          ti::convert_ti_object<V>(D0_i.get_type()), false)
    , m_A(A), m_Dm1_i(Dm1_i), m_D0_i(D0_i), m_Dp1_i(Dp1_i), m_Dp2_i(Dp2_i), m_piv(piv)
{
    Real tol        = opts.get_option<Real>(opt::linsolve::tol_sing());

    // calculate small pivot tolerance
    if (tol > 0.0)
    {
        Real max_0  = norm_vec_all(Matrix(D0_i,false), basic_vector_norm::norm_inf);
        Real max_1  = norm_vec_all(Matrix(Dp1_i,false), basic_vector_norm::norm_inf);
        Real max_2  = norm_vec_all(Matrix(Dp2_i,false), basic_vector_norm::norm_inf);
        Real max    = std::max(std::max(max_0, max_1), max_2);
            
        if (max == 0.0)
            throw error::error_singular();

        tol         = pow(constants::eps<VR>(), tol) * max;
    }
    else if (tol < 0)
    {
        tol         = -tol;
    };

    VR tol_t        = (VR)tol;

    // perturb small pivots
    if (tol != 0.0)
    {
        bool need_modif = false;
        Integer N       = m_D0_i.length();
        const V* ptr_D  = m_D0_i.ptr();
        Integer j;

        for (j = 0; j < N; ++j)
        {
            V val   = ptr_D[j];

            if (abs(val) < tol_t)
            {
                need_modif = true;
                break;
            };
        }

        if (need_modif == true)
        {
            m_D0_i.assign_to_fresh(m_D0_i.make_unique());
            V* ptr_D2   = m_D0_i.ptr();

            for (; j < N; ++j)
            {
                V val   = ptr_D2[j];

                if (abs(val) < tol_t)
                    ptr_D2[j]   = tol_t * (val == 0 ? VR(1.0) : sign(ptr_D2[j]));
            }            

            m_modified  = true;
        };
    };
};

template<class V>
linsolve_obj_tridiag_fac<V>::linsolve_obj_tridiag_fac(const Mat_B& A, const Mat& Dm1_i, const Mat& D0_i, 
            const Mat& Dp1_i, const Mat& Dp2_i, const Mat_I& piv, bool modif, from_inv)
    :linsolve_obj_base(D0_i.length(),D0_i.length(),matrix_traits::value_code<V>::value,
                          ti::convert_ti_object<V>(D0_i.get_type()), modif)
    , m_A(A), m_Dm1_i(Dm1_i), m_D0_i(D0_i), m_Dp1_i(Dp1_i), m_Dp2_i(Dp2_i), m_piv(piv)
{};

template<class V>
linsolve_obj_tridiag_fac<V>::~linsolve_obj_tridiag_fac()
{};

template<class V>
matcl::Matrix linsolve_obj_tridiag_fac<V>::inv() const
{
    Integer N   = this->rows();
    Matrix X    = eye(N, N, this->get_value_code()); 
    X           = this->solve(std::move(X), trans_type::no_trans);

    return X;
};

template<class V>
matcl::Matrix linsolve_obj_tridiag_fac<V>::solve(const Matrix& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_tridiag_fac;

    Matrix X(X0);

    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_tridiag_fac<V>::solve(Matrix&& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_tridiag_fac;

    Matrix X(std::move(X0));

    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_tridiag_fac<V>::solve_rev(const Matrix& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_tridiag_fac;

    Matrix X(X0);

    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve_rev<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
matcl::Matrix linsolve_obj_tridiag_fac<V>::solve_rev(Matrix&& X0, trans_type tA) const
{
    using Fac       = linsolve_obj_tridiag_fac;

    Matrix X(std::move(X0));

    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    value_code vc   = matrix_traits::unify_value_types(this->get_value_code(), X.get_value_code());    
    Matrix t        = zeros(0,0,vc);    

    Matrix ret;
    linsolve_obj_solve_rev<V,Fac>::make<const Matrix&>(t,*this, ret, X, tA);
    return ret;
};

template<class V>
template<class T>
void linsolve_obj_tridiag_fac<V>::solve_impl(Matrix& ret, Matrix& X, trans_type tA) const
{
    Integer N       = m_N;
    Integer Nrhs    = X.cols();

    using M         = raw::Matrix<T,struct_dense>;

    const M& Dm1c   = raw::converter<M,Mat>::eval(m_Dm1_i);
    const M& D0c    = raw::converter<M,Mat>::eval(m_D0_i);
    const M& Dp1c   = raw::converter<M,Mat>::eval(m_Dp1_i);
    const M& Dp2c   = raw::converter<M,Mat>::eval(m_Dp2_i);

    M Bc            = X.impl_unique<M>();

    Integer info;

    using VL        = typename md::lapack_value_type<T>::type;

    lapack::gttrs<VL>(get_trans_code(tA), N, Nrhs, lap(Dm1c.ptr()), lap(D0c.ptr()), lap(Dp1c.ptr()), 
                  lap(Dp2c.ptr()), m_piv.ptr(), lap(Bc.ptr()), Bc.ld(), info);

    if (info != 0)
        throw error::error_general("invalid argument passed to gttrs");

    Bc.get_struct().reset();
    ret = Matrix(Bc,true);
    return;
};

template<class V>
template<class T>
void linsolve_obj_tridiag_fac<V>::solve_rev_impl(Matrix& ret, Matrix& X, trans_type tA) const
{
    Integer N       = m_N;
    Integer Mrhs    = X.rows();

    using M         = raw::Matrix<T,struct_dense>;
    
    const M& Dm1c   = raw::converter<M,Mat>::eval(m_Dm1_i);
    const M& D0c    = raw::converter<M,Mat>::eval(m_D0_i);
    const M& Dp1c   = raw::converter<M,Mat>::eval(m_Dp1_i);
    const M& Dp2c   = raw::converter<M,Mat>::eval(m_Dp2_i);

    M Bc            = X.impl_unique<M>();

    using VL        = typename md::lapack_value_type<T>::type;

    Integer info;

    lapack::gttrs_rev<VL>(get_trans_code(tA), N, Mrhs, lap(Dm1c.ptr()), lap(D0c.ptr()), lap(Dp1c.ptr()), 
                  lap(Dp2c.ptr()), m_piv.ptr(), lap(Bc.ptr()), Bc.ld(), info);

    if (info != 0)
        throw error::error_general("invalid argument passed to gttrs_rev");

    Bc.get_struct().reset();
    ret = Matrix(Bc,true);
    return;
};

template<class V>
template<class T>
typename linsolve_obj_tridiag_fac<V>::data_ptr 
linsolve_obj_tridiag_fac<V>::convert_impl() const
{
    using Mat_C     = raw::Matrix<T,struct_dense>;
    using Mat_BC    = raw::Matrix<T,struct_banded>;

    Mat_BC Ac       = raw::converter<Mat_BC, Mat_B>::eval(m_A);
    Mat_C Dm1_ic    = raw::converter<Mat_C, Mat>::eval(m_Dm1_i);    
    Mat_C D0_ic     = raw::converter<Mat_C, Mat>::eval(m_D0_i);    
    Mat_C Dp1_ic    = raw::converter<Mat_C, Mat>::eval(m_Dp1_i);    
    Mat_C Dp2_ic    = raw::converter<Mat_C, Mat>::eval(m_Dp2_i);    

    return data_ptr(new linsolve_obj_tridiag_fac<T>(Ac, Dm1_ic, D0_ic, Dp1_ic, Dp2_ic, m_piv, 
                                                    m_modified, from_inv()));
};

template<class V>
Real linsolve_obj_tridiag_fac<V>::log_det() const
{
    Integer M       = m_N;

    const V* ptr    = m_D0_i.ptr();
    Integer st      = m_D0_i.cols() == M ? m_D0_i.ld() : 1;

    Real det        = 0.0;

    for (Integer i = 0; i < M; ++i)
    {
        det         += matcl::log(matcl::abs(ptr[0]));
        ptr         += st;
    };
    
    return det;
};

template<class V>
matcl::Matrix linsolve_obj_tridiag_fac<V>::base_matrix() const
{
    return Matrix(m_A,false);
};
template<class V>
Real linsolve_obj_tridiag_fac<V>::mat_normest_1() const
{
    return matcl::norm(Matrix(m_A,false), 1.0);
};

template<class V>
Real linsolve_obj_tridiag_fac<V>::mat_normest_inf() const
{
    return matcl::norm(Matrix(m_A,false), constants::inf());
};

template<class V>
Matrix linsolve_obj_tridiag_fac<V>::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return mmul(Matrix(m_A, false), X, t);
};

template<class V>
Matrix linsolve_obj_tridiag_fac<V>::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return mmul(Matrix(m_A, false), std::move(X), t);
};

template<class V>
Matrix linsolve_obj_tridiag_fac<V>::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return mmul(X, Matrix(m_A,false), trans_type::no_trans,t);
};

template<class V>
Matrix linsolve_obj_tridiag_fac<V>::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return mmul(std::move(X), Matrix(m_A,false), trans_type::no_trans,t);
};

template<class V>
typename linsolve_obj_tridiag_fac<V>::data_ptr 
linsolve_obj_tridiag_fac<V>::convert(value_code vc0) const
{
    using Fac       = linsolve_obj_tridiag_fac;
    value_code vc   = matrix_traits::unify_value_types(vc0, value_code::v_float);
    Matrix dum      = zeros(0,0, vc);

    data_ptr ret;
    linsolve_obj_convert<V,Fac>::make<const Matrix&>(dum,*this, ret);
    return ret;
};

template class linsolve_obj_tridiag_fac<Real>;
template class linsolve_obj_tridiag_fac<Float>;
template class linsolve_obj_tridiag_fac<Complex>;
template class linsolve_obj_tridiag_fac<Float_complex>;

//--------------------------------------------------------------------------------
//                  linsolve_forwarding
//--------------------------------------------------------------------------------
linsolve_forwarding::linsolve_forwarding(const linsolve_obj& base)
    :linsolve_obj_base(base.rows(), base.cols(), base.get_value_code(), base.get_type(), base.is_modified())
    ,m_base(base)
{};

linsolve_forwarding::data_ptr linsolve_forwarding::convert(value_code new_val_code) const
{
    return data_ptr(new linsolve_forwarding(m_base.convert(new_val_code)));
};
matcl::Matrix linsolve_forwarding::solve(const Matrix& X, trans_type tA) const
{
    return m_base.solve(X, tA);
}
matcl::Matrix linsolve_forwarding::solve(Matrix&& X, trans_type tA) const
{
    return m_base.solve(std::move(X), tA);
}
matcl::Matrix linsolve_forwarding::solve_rev(const Matrix& X, trans_type tA) const
{
    return m_base.solve_rev(X, tA);
}
matcl::Matrix linsolve_forwarding::solve_rev(Matrix&& X, trans_type tA) const
{
    return m_base.solve_rev(std::move(X), tA);
}
Matrix linsolve_forwarding::mmul_right(const Matrix& X, trans_type t) const
{
    return m_base.mmul_right(X, t);
}
Matrix linsolve_forwarding::mmul_right(Matrix&& X, trans_type t) const
{
    return m_base.mmul_right(std::move(X), t);
}
Matrix linsolve_forwarding::mmul_left(const Matrix& X, trans_type t) const
{
    return m_base.mmul_left(X, t);
}
Matrix linsolve_forwarding::mmul_left(Matrix&& X, trans_type t) const
{
    return m_base.mmul_left(std::move(X), t);
}

//--------------------------------------------------------------------------------
//                  linsolve_perm
//--------------------------------------------------------------------------------

linsolve_perm::linsolve_perm(const Matrix& A, const linsolve_obj& Apq, const permvec& p, const permvec& q,
                bool symperm)
    :linsolve_forwarding(Apq), m_A(A), m_has_A(true), m_p(p), m_q(q), m_symperm(symperm)
{};

linsolve_perm::linsolve_perm(const linsolve_obj& Apq, const permvec& p, const permvec& q, bool symperm)
    :linsolve_forwarding(Apq), m_has_A(false), m_p(p), m_q(q), m_symperm(symperm)
{};

linsolve_perm::~linsolve_perm()
{};

bool linsolve_perm::is_hermitian() const
{
    if (m_symperm == false)
        return false;
    else
        return base().is_hermitian();
};

bool linsolve_perm::is_posdef() const
{
    if (m_symperm == false)
        return false;
    else
        return base().is_posdef();
};

linsolve_perm::data_ptr linsolve_perm::convert(value_code nvc) const
{
    if (m_has_A)
    {
        Matrix Ac   = details::convert_value(m_A, nvc);
        return data_ptr(new linsolve_perm(Ac, base().convert(nvc), m_p, m_q, m_symperm) );
    }
    else
    {
        return data_ptr(new linsolve_perm(base().convert(nvc), m_p, m_q, m_symperm) );
    }
};
Matrix linsolve_perm::base_matrix() const
{
    if (m_has_A)
        return m_A;

    Matrix Ab   = base().base_matrix();
    return Ab(m_p.invperm(), m_q.invperm());
};

Matrix linsolve_perm::inv() const
{
    Matrix Aib  = base().inv();
    Matrix Ai   = Aib(m_q.invperm(), m_p.invperm());
    return Ai;
};

Matrix linsolve_perm::mmul_right(const Matrix& X, trans_type t) const
{
    if (m_has_A)
        return mmul(m_A, X, t);

    if (t == trans_type::no_trans)
    {
        Matrix Apx  = base().mmul_right(X(m_q, colon()),t);
        Matrix Ax   = std::move(Apx)(m_p.invperm(), colon());
        return Ax;
    }
    else
    {
        Matrix Apx  = base().mmul_right(X(m_p, colon()),t);
        Matrix Ax   = std::move(Apx)(m_q.invperm(), colon());
        return Ax;
    };
};
Matrix linsolve_perm::mmul_right(Matrix&& X, trans_type t) const
{
    if (m_has_A)
        return mmul(m_A, std::move(X), t);

    if (t == trans_type::no_trans)
    {
        Matrix Apx  = base().mmul_right(std::move(X)(m_q, colon()),t);
        Matrix Ax   = std::move(Apx)(m_p.invperm(), colon());
        return Ax;
    }
    else
    {
        Matrix Apx  = base().mmul_right(std::move(X)(m_p, colon()),t);
        Matrix Ax   = std::move(Apx)(m_q.invperm(), colon());
        return Ax;
    };
};
Matrix linsolve_perm::mmul_left(const Matrix& X, trans_type t) const
{
    if (m_has_A)
        return mmul(X, m_A, trans_type::no_trans, t);

    if (t == trans_type::no_trans)
    {
        Matrix Apx  = base().mmul_left(X(colon(), m_p),t);
        Matrix Ax   = std::move(Apx)(colon(), m_q.invperm());
        return Ax;
    }
    else
    {
        Matrix Apx  = base().mmul_left(X(colon(), m_q),t);
        Matrix Ax   = std::move(Apx)(colon(), m_p.invperm());
        return Ax;
    };
};
Matrix linsolve_perm::mmul_left(Matrix&& X, trans_type t) const
{
    if (m_has_A)
        return mmul(std::move(X), m_A, trans_type::no_trans, t);
    
    if (t == trans_type::no_trans)
    {
        Matrix Apx  = base().mmul_left(std::move(X)(colon(), m_p),t);
        Matrix Ax   = std::move(Apx)(colon(), m_q.invperm());
        return Ax;
    }
    else
    {
        Matrix Apx  = base().mmul_left(std::move(X)(colon(), m_q),t);
        Matrix Ax   = std::move(Apx)(colon(), m_p.invperm());
        return Ax;
    };
};

Matrix linsolve_perm::solve(const Matrix& X, trans_type t) const
{
    if (t == trans_type::no_trans)
    {
        Matrix Apx  = base().solve(X(m_p, colon()),t);
        Matrix Ax   = std::move(Apx)(m_q.invperm(), colon());
        return Ax;
    }
    else
    {
        Matrix Apx  = base().solve(X(m_q, colon()),t);
        Matrix Ax   = std::move(Apx)(m_p.invperm(), colon());
        return Ax;
    };
};

Matrix linsolve_perm::solve(Matrix&& X, trans_type t) const
{
    if (t == trans_type::no_trans)
    {
        Matrix Apx  = base().solve(std::move(X)(m_p, colon()),t);
        Matrix Ax   = std::move(Apx)(m_q.invperm(), colon());
        return Ax;
    }
    else
    {
        Matrix Apx  = base().solve(std::move(X)(m_q, colon()),t);
        Matrix Ax   = std::move(Apx)(m_p.invperm(), colon());
        return Ax;
    };
};
Matrix linsolve_perm::solve_rev(const Matrix& X, trans_type t) const
{
    if (t == trans_type::no_trans)
    {
        Matrix Apx  = base().solve_rev(X(m_q, colon()),t);
        Matrix Ax   = std::move(Apx)(m_p.invperm(), colon());
        return Ax;
    }
    else
    {
        Matrix Apx  = base().solve_rev(X(m_p, colon()),t);
        Matrix Ax   = std::move(Apx)(m_q.invperm(), colon());
        return Ax;
    };
};

Matrix linsolve_perm::solve_rev(Matrix&& X, trans_type t) const
{
    if (t == trans_type::no_trans)
    {
        Matrix Apx  = base().solve_rev(std::move(X)(m_q, colon()),t);
        Matrix Ax   = std::move(Apx)(m_p.invperm(), colon());
        return Ax;
    }
    else
    {
        Matrix Apx  = base().solve_rev(std::move(X)(m_p, colon()),t);
        Matrix Ax   = std::move(Apx)(m_q.invperm(), colon());
        return Ax;
    };
};


};};