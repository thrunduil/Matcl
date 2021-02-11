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

#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type2_switch.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"

#include "matcl-linalg/decompositions/lu.h"
#include "matcl-linalg/decompositions/lu/lu_impl.h"

namespace matcl { namespace details
{

template<class S, class V>
struct linsolve_general_lu_struct
{};

//--------------------------------------------------------------------------------
//                  GENERAL SPARSE
//--------------------------------------------------------------------------------

struct lu_solve
{
    static void eval(Matrix& ret, const Matrix& L, const Matrix& U, const permvec& P,
                     const permvec& Q, const permvec& p, const permvec& q,
                     Matrix&& B, trans_type trans)
    {
        permvec id = permvec::identity(L.rows());

        if (trans == trans_type::no_trans)
        {
            Matrix X;

            if (p.is_id() == true)
                X = linsolve(L,P,id, std::move(B), trans_type::no_trans);
            else
                X = linsolve(L,p(P),id, std::move(B), trans_type::no_trans);

            if (q.is_id() == true)
                X = linsolve(U,id,Q, std::move(X), trans_type::no_trans);
            else
                X = linsolve(U,id,q(Q), std::move(X), trans_type::no_trans);

            ret = X;
            return;
        }
        else
        {
            Matrix X;

            if (q.is_id() == true)
                X = linsolve(U,id,Q, std::move(B), trans);            
            else
                X = linsolve(U,id,q(Q), std::move(B), trans);            

            if (p.is_id() == true)
                X = linsolve(L,P,id, std::move(X), trans);
            else
                X = linsolve(L,p(P),id, std::move(X), trans);

            ret = X;
            return;
        };
    };

    static void eval_rev(Matrix& ret, const Matrix& L, const Matrix& U, const permvec& P,
                     const permvec& Q, const permvec& p, const permvec& q,
                     Matrix&& B, trans_type trans)
    {
        permvec id = permvec::identity(L.rows());

        if (trans == trans_type::no_trans)
        {
            Matrix X;

            if (q.is_id() == true)
                X = linsolve_rev(U,id,Q, std::move(B));
            else
                X = linsolve_rev(U,id,q(Q), std::move(B));

            if (p.is_id() == true)
                X = linsolve_rev(L,P,id, std::move(X));
            else
                X = linsolve_rev(L,p(P),id, std::move(X));

            ret = X;
            return;
        }
        else
        {
            Matrix X;

            if (p.is_id() == true)
                X = linsolve_rev2(L,P,id, std::move(B), trans);
            else
                X = linsolve_rev2(L,p(P),id, std::move(B), trans);

            if (q.is_id() == true)
                X = linsolve_rev2(U,id,Q, std::move(X), trans);
            else
                X = linsolve_rev2(U,id,q(Q), std::move(X), trans);

            ret = X;
            return;
        };        
    };
};

template<class struct_type, class V1, class V2>
struct linsolve_general_sparse_impl
{
    using M1    = raw::Matrix<V1,struct_sparse>;
    using M2    = raw::Matrix<V2,struct_type>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, trans_type trans,
                     const options& opts)
    {        
        matcl::Matrix L, U;
        permvec P, Q;

        {
            lu_return_type ret_loc; 
            lu_sparse<V1>::eval(ret_loc, A, opts);
            std::tie(L,U,P,Q) = ret_loc;
        };

        if (U.rows() != U.cols())
            throw error::square_matrix_required(U.rows(), U.cols());

        return lu_solve::eval(ret,L,U,P,Q,p,q, Matrix(B,false), trans);
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, trans_type trans,
                         const options& opts)
    {
        matcl::Matrix L, U;
        permvec P, Q;

        {
            lu_return_type ret_loc;
            lu_sparse<V1>::eval(ret_loc, A,opts);
            std::tie(L,U,P,Q) = ret_loc;
        }

        if (U.rows() != U.cols())
            throw error::square_matrix_required(U.rows(), U.cols());

        return lu_solve::eval_rev(ret,L,U,P,Q,p,q, Matrix(B,false), trans);
    };
};

template<class V>
struct linsolve_general_lu_struct<struct_sparse,V>
{
    using M1    = raw::Matrix<V,struct_sparse>;

    static void eval(linsolve_obj& ret, const M1& A, const options& opts)
    {   
        matcl::Matrix L, U;
        permvec P, Q;

        matcl::Matrix As(A,false);

        {
            lu_return_type ret_loc; 
            lu_sparse<V>::eval(ret_loc, A, opts);
            std::tie(L,U,P,Q) = ret_loc;
        };

        if (U.rows() != U.cols())
            throw error::square_matrix_required(U.rows(), U.cols());

        bool isv = L.all_finite() && U.all_finite();

        if (isv == false)
        {
            ret = linsolve_nan(As.rows(), matrix_traits::value_code<V>::value);
            return;
        };

        using linsolve_ptr  = linsolve_obj::linsolve_data_ptr;
        ret = linsolve_obj(linsolve_ptr(new details::linsolve_obj_lu_factors(As, L, U, P, Q, opts)));
        return;
    };
};

template<class struct_type, class V1, class V2>
struct linsolve_general_sparse{};

//--------------------------------------------------------------------------------
//                  GENERAL SPARSE-DENSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_sparse<struct_dense,V1,V2> 
    : public linsolve_general_sparse_impl<struct_dense,V1,V2>
{};

//--------------------------------------------------------------------------------
//                  GENERAL SPARSE-BAND
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_sparse<struct_banded,V1,V2>
    : public linsolve_general_sparse_impl<struct_banded,V1,V2>
{};

//--------------------------------------------------------------------------------
//                  GENERAL SPARSE-SPARSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_sparse<struct_sparse,V1,V2>
    : public linsolve_general_sparse_impl<struct_sparse,V1,V2>
{};

//--------------------------------------------------------------------------------
//                  STRUCT TYPES
//--------------------------------------------------------------------------------
template<class S1, class S2, class V1, class V2>
struct linsolve_general_struct
{};

//--------------------------------------------------------------------------------
//                  DENSE-DENSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_dense,struct_dense,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_dense>;
    using M2    = raw::Matrix<V2,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                     matcl::permvec q, const M2& B, trans_type trans,
                     const options& opts)
    {   
        using pivot_type    = opt::linsolve::pivot_type;
        Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
        pivot_type piv      = (pivot_type)piv_int;
        bool isp            = (piv == pivot_type::partial);

        if (isp == false)
        {
            matcl::Matrix L, U;
            permvec P, Q;

            {
                lu_return_type ret_loc; 
                lu_dense<V1>::eval(ret_loc, A, opts);
                std::tie(L,U,P,Q) = ret_loc;
            };

            if (U.rows() != U.cols())
                throw error::square_matrix_required(U.rows(), U.cols());

            return lu_solve::eval(ret,L,U,P,Q,p,q, Matrix(B,false), trans);
        };

        Integer N                   = A.rows();
        Integer Nrhs                = B.cols();

        matcl::lapack::i_type info  = 0;

        raw::integer_dense ipiv(ti::ti_empty(),1,N);

        M1 A2 = A.make_unique();

        // Compute the LU factorization of A.
        lapack::getrf_rec( N, N, lap(A2.ptr()), A2.ld(), lap(ipiv.ptr()), &info );
        
        if (info != 0)
            throw error::error_singular();

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_dense>;

        const M& Ac = raw::converter<M,M1>::eval(A2);
        M Bc        = M(raw::converter<M,M2>::eval(B.make_unique()), M::copy_is_safe());

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute rows of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswp(Bc.cols(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), 1);
        };

        // Solve the system A*X = B, overwriting B with X.
        lapack::getrs( get_trans_code(trans), N, Nrhs, lap(Ac.ptr()), Ac.ld(), lap(ipiv.ptr()), 
                    lap(Bc.ptr()), Bc.ld(), &info);

        if (info != 0)
            throw error::error_singular();

        //permute rows of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswp(Bc.cols(), lap(Bc.ptr()), Bc.ld(), 1, q.length(), lap(q_ptr), -1);
        };

        Bc.get_struct().reset();
        ret = Matrix(Bc,true);
        return;
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                         matcl::permvec q, const M2& B, trans_type trans,
                         const options& opts)
    {
        using pivot_type    = opt::linsolve::pivot_type;
        Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
        pivot_type piv      = (pivot_type)piv_int;
        bool isp            = (piv == pivot_type::partial);

        if (isp == false)
        {
            matcl::Matrix L, U;
            permvec P, Q;

            {
                lu_return_type ret_loc; 
                lu_dense<V1>::eval(ret_loc, A, opts);
                std::tie(L,U,P,Q) = ret_loc;
            };

            if (U.rows() != U.cols())
                throw error::square_matrix_required(U.rows(), U.cols());

            return lu_solve::eval_rev(ret,L,U,P,Q,p,q, Matrix(B,false), trans);
        };

        Integer N                   = A.rows();
        Integer Mrhs                = B.rows();
        matcl::lapack::i_type info  = 0;

        raw::integer_dense ipiv(ti::ti_empty(),1,N);

        M1 A2 = A.make_unique();

        // Compute the LU factorization of A.
        lapack::getrf_rec( N, N, lap(A2.ptr()), A2.ld(), lap(ipiv.ptr()), &info );
        
        if (info != 0)
            throw error::error_singular();

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_dense>;

        const M& Ac = raw::converter<M,M1>::eval(A2);
        M Bc        = M(raw::converter<M,M2>::eval(B.make_unique()), M::copy_is_safe());

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute columns of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswpc(Bc.rows(), lap(Bc.ptr()), Bc.ld(), 1, q.length(), lap(q_ptr), 1);
        };

        // Solve the system X*A = B, overwriting B with X.
        lapack::getrs_rev( get_trans_code(trans), N, Mrhs, lap(Ac.ptr()), Ac.ld(), lap(ipiv.ptr()), 
                    lap(Bc.ptr()), Bc.ld(), &info);

        if (info != 0)
            throw error::error_singular();

        //permute columns of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswpc(Bc.rows(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), -1);
        };

        Bc.get_struct().reset();
        ret = Matrix(Bc,true);
        return;
    };
};

template<class V>
struct linsolve_general_lu_struct<struct_dense,V>
{
    using M1    = raw::Matrix<V,struct_dense>;

    static void eval(linsolve_obj& ret, const M1& A, const options& opts)
    {   
        using pivot_type    = opt::linsolve::pivot_type;
        Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
        pivot_type piv      = (pivot_type)piv_int;
        bool isp            = (piv == pivot_type::partial);        

        if (isp == false)
        {
            matcl::Matrix L, U;
            permvec P, Q;

            //increase refcount
            matcl::Matrix As(A,false);

            {
                lu_return_type ret_loc; 
                lu_dense<V>::eval(ret_loc, A, opts);
                std::tie(L,U,P,Q) = ret_loc;
            };

            if (U.rows() != U.cols())
                throw error::square_matrix_required(U.rows(), U.cols());

            bool isv = L.all_finite() && U.all_finite();

            if (isv == false)
            {
                ret = linsolve_nan(A.rows(), matrix_traits::value_code<V>::value);
                return;
            };

            using linsolve_ptr  = linsolve_obj::linsolve_data_ptr;
            ret = linsolve_obj(linsolve_ptr(new details::linsolve_obj_lu_factors(As, L, U, P, Q, opts)));
            return;
        };

        Integer N       = A.rows();
        Integer info    = 0;

        raw::integer_dense ipiv(ti::ti_empty(), 1, N);

        M1 A2           = A.copy();

        // Compute the LU factorization of A.
        lapack::getrf_rec( N, N, lap(A2.ptr()), A2.ld(), lap(ipiv.ptr()), &info );
       
        bool isv = A2.all_finite();

        if (isv == false)
        {
            ret = linsolve_nan(A.rows(), matrix_traits::value_code<V>::value);
            return;
        };

        using linsolve_ptr  = linsolve_obj::linsolve_data_ptr;

        ret = linsolve_obj(linsolve_ptr(new details::linsolve_obj_lu_dense<V>(A,A2, ipiv, opts)));
        return;
    };
};

//--------------------------------------------------------------------------------
//                  DENSE-SPARSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_dense,struct_sparse,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_dense>;
    using M2    = raw::Matrix<V2,struct_sparse>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, trans_type trans,
                     const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_general_struct<struct_dense,struct_dense,V1,V2>::eval(ret, A,p,q, Bc,trans, opts);
    };
    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, trans_type trans,
                         const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_general_struct<struct_dense,struct_dense,V1,V2>::eval_rev(ret,A,p,q, Bc,trans,opts);
    };
};

//--------------------------------------------------------------------------------
//                  DENSE-BAND
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_dense,struct_banded,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_dense>;
    using M2    = raw::Matrix<V2,struct_banded>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, trans_type trans,
                     const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_general_struct<struct_dense,struct_dense,V1,V2>::eval(ret, A,p,q, Bc,trans, opts);
    };
    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, trans_type trans,
                         const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_general_struct<struct_dense,struct_dense,V1,V2>::eval_rev(ret,A,p,q, Bc,trans,opts);
    };
};

//--------------------------------------------------------------------------------
//                  SPARSE-DENSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_sparse,struct_dense,V1,V2>
    : public linsolve_general_sparse<struct_dense, V1, V2>
{};

//--------------------------------------------------------------------------------
//                  BAND-DENSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_banded,struct_dense,V1,V2>
    : public linsolve_general_sparse<struct_dense, V1, V2>
{
    using M1    = raw::Matrix<V1,struct_banded>;
    using M2    = raw::Matrix<V2,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                     matcl::permvec q, const M2& B, trans_type trans,
                     const options& opts)
    {             
        using pivot_type    = opt::linsolve::pivot_type;
        Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
        pivot_type piv      = (pivot_type)piv_int;
        bool isp            = (piv == pivot_type::partial);

        if (isp == false)
        {
            matcl::Matrix L, U;
            permvec P, Q;

            {
                lu_return_type ret_loc; 
                lu_band<V1>::eval(ret_loc, A, opts);
                std::tie(L,U,P,Q) = ret_loc;
            };

            if (U.rows() != U.cols())
                throw error::square_matrix_required(U.rows(), U.cols());

            return lu_solve::eval(ret,L,U,P,Q,p,q, Matrix(B,false), trans);
        };

        using V     = typename unify_types<V1,V2>::type;
        using M     = raw::Matrix<V,struct_dense>;
        using MB    = raw::Matrix<V,struct_banded>;

        if (A.has_diag(0) == false)
            throw error::error_singular();

        Integer ld      = A.number_subdiagonals();
        Integer ud      = A.number_superdiagonals();
        Integer N       = A.cols();
        Integer Nrhs    = B.cols();

        using Mat_B     = raw::Matrix<V1,struct_banded>;

        Mat_B AFM(A.get_type(),N,max(N, ld + ud + 1), -ld,ud+ld);

        using V1T       = typename lapack_value_type<V1>::type;

        const V1T* ptr_A    = lap(A.rep_ptr());
        V1T* ptr_AF         = lap(AFM.rep_ptr());

        for (Integer J = 1; J <= N; ++J)
        {
            Integer J1 = max( J-ud, 1 );
            Integer J2 = min( J+ld, N );

            matcl::lapack::copy(J2-J1+1,ptr_A+ud-J+J1,1,ptr_AF+ld+ud-J+J1,1);

            ptr_A   += A.ld();
            ptr_AF  += AFM.ld();
        };
        
        raw::integer_dense ipiv(ti::ti_empty(),1,N);

        matcl::lapack::i_type info = 0;

        //factor
        lapack::gbtrf( N, N, ld, ud, lap(AFM.rep_ptr()), AFM.ld(), lap(ipiv.ptr()), &info );

        if (info != 0)
            throw error::error_singular();

        //solve
        const MB& Ac    = raw::converter<MB,M1>::eval(AFM);
        const M& Bc0    = raw::converter<M,M2>::eval(B.make_unique());
        M Bc            = M(Bc0, M::copy_is_safe());

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute rows of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswp(Bc.cols(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), 1);
        };

        lapack::gbtrs( get_trans_code(trans), N, ld, ud, Nrhs, lap(Ac.rep_ptr()), Ac.ld(), lap(ipiv.ptr()),
                     lap(Bc.ptr()), Bc.ld(), &info );

        if (info != 0)
            throw error::error_singular();

        //permute rows of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswp(Bc.cols(), lap(Bc.ptr()), Bc.ld(), 1, q.length(), lap(q_ptr), -1);
        };	

        Bc.get_struct().reset();
        ret = Matrix(Bc,true);
        return;
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                         matcl::permvec q, const M2& B, trans_type trans,
                         const options& opts)
    {
        using pivot_type    = opt::linsolve::pivot_type;
        Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
        pivot_type piv      = (pivot_type)piv_int;
        bool isp            = (piv == pivot_type::partial);

        if (isp == false)
        {
            matcl::Matrix L, U;
            permvec P, Q;

            {
                lu_return_type ret_loc; 
                lu_band<V1>::eval(ret_loc, A, opts);
                std::tie(L,U,P,Q) = ret_loc;
            };

            if (U.rows() != U.cols())
                throw error::square_matrix_required(U.rows(), U.cols());

            return lu_solve::eval_rev(ret,L,U,P,Q,p,q, Matrix(B,false), trans);
        };

        using V     = typename unify_types<V1,V2>::type;
        using M     = raw::Matrix<V,struct_dense>;
        using MB    = raw::Matrix<V,struct_banded>;

        if (A.has_diag(0) == false)
            throw error::error_singular();

        Integer ld      = A.number_subdiagonals();
        Integer ud      = A.number_superdiagonals();
        Integer N       = A.cols();
        Integer Mrhs    = B.rows();

        using Mat_B     = raw::Matrix<V1,struct_banded>;

        Mat_B AFM(A.get_type(),N,max(N, ld + ud + 1),-ld,ud+ld);

        using V1T       = typename lapack_value_type<V1>::type;

        const V1T* ptr_A    = lap(A.rep_ptr());
        V1T* ptr_AF         = lap(AFM.rep_ptr());

        for (Integer J = 1; J <= N; ++J)
        {
            Integer J1 = max( J-ud, 1 );
            Integer J2 = min( J+ld, N );
            matcl::lapack::copy(J2-J1+1,ptr_A+ud-J+J1,1,ptr_AF+ld+ud-J+J1,1);

            ptr_A   += A.ld();
            ptr_AF  += AFM.ld();
        };
        
        raw::integer_dense ipiv(ti::ti_empty(),1,N);

        matcl::lapack::i_type info = 0;

        //factor
        lapack::gbtrf( N, N, ld, ud, lap(AFM.rep_ptr()), AFM.ld(), lap(ipiv.ptr()), &info );

        if (info != 0)
            throw error::error_singular();

        //solve
        const MB& Ac    = raw::converter<MB,M1>::eval(AFM);
        const M& Bc0    = raw::converter<M,M2>::eval(B.make_unique());
        M Bc            = M(Bc0, M::copy_is_safe());

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute columns of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswpc(Bc.rows(), lap(Bc.ptr()), Bc.ld(), 1, q.length(), lap(q_ptr), 1);
        };

        lapack::gbtrs_rev( get_trans_code(trans), N, ld, ud, Mrhs, lap(Ac.rep_ptr()), Ac.ld(), lap(ipiv.ptr()),
                     lap(Bc.ptr()), Bc.ld(), &info );

        if (info != 0)
            throw error::error_singular();

        //permute columns of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswpc(Bc.rows(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), -1);
        };

        Bc.get_struct().reset();
        ret = Matrix(Bc,true);
        return;
    };
};

template<class V>
struct linsolve_general_lu_struct<struct_banded,V>
{
    using M1    = raw::Matrix<V,struct_banded>;

    static void eval(linsolve_obj& ret, const M1& A, const options& opts)
    {   
        using pivot_type    = opt::linsolve::pivot_type;
        Integer piv_int     = opts.get_option<Integer>(opt::linsolve::pivot());
        pivot_type piv      = (pivot_type)piv_int;
        bool isp            = (piv == pivot_type::partial);

        if (isp == false || A.has_diag(0) == false)
        {
            Matrix As(A,false);

            matcl::Matrix L, U;
            permvec P, Q;

            {
                lu_return_type ret_loc; 
                lu_band<V>::eval(ret_loc, A, opts);
                std::tie(L,U,P,Q) = ret_loc;
            };

            bool isv = L.all_finite() && U.all_finite();

            if (isv == false)
            {
                ret = linsolve_nan(As.rows(), matrix_traits::value_code<V>::value);
                return;
            };

            using linsolve_ptr  = linsolve_obj::linsolve_data_ptr;
            ret = linsolve_obj(linsolve_ptr(new details::linsolve_obj_lu_factors(As, L, U, P, Q, opts)));
            return;
        };

        Integer ld      = A.number_subdiagonals();
        Integer ud      = A.number_superdiagonals();
        Integer N       = A.cols();

        using Mat_B     = raw::Matrix<V,struct_banded>;
        Mat_B AFM(A.get_type(),N, max(N, ld + ud + 1), -ld, ud+ld);

        using V1T       = typename lapack_value_type<V>::type;

        const V1T* ptr_A    = lap(A.rep_ptr());
        V1T* ptr_AF         = lap(AFM.rep_ptr());

        for (Integer J = 1; J <= N; ++J)
        {
            Integer J1 = max( J-ud, 1 );
            Integer J2 = min( J+ld, N );

            matcl::lapack::copy(J2-J1+1,ptr_A+ud-J+J1,1,ptr_AF+ld+ud-J+J1,1);

            ptr_A   += A.ld();
            ptr_AF  += AFM.ld();
        };
        
        raw::integer_dense ipiv(ti::ti_empty(),1,N);

        matcl::lapack::i_type info = 0;

        //factor
        lapack::gbtrf( N, N, ld, ud, lap(AFM.rep_ptr()), AFM.ld(), lap(ipiv.ptr()), &info );

        if (info < 0)
            throw error::error_general("invalid argument passed to gbtrf");

        bool isv = AFM.all_finite();

        if (isv == false)
        {
            ret = linsolve_nan(A.rows(), matrix_traits::value_code<V>::value);
            return;
        };

        using linsolve_ptr  = linsolve_obj::linsolve_data_ptr;
        ret = linsolve_obj(linsolve_ptr(new details::linsolve_obj_lu_band<V>(A, AFM, ipiv, N, ld, ud, opts)));
        return;
    };
};

//--------------------------------------------------------------------------------
//                  BAND-SPARSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_banded,struct_sparse,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_banded>;
    using M2    = raw::Matrix<V2,struct_sparse>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, trans_type trans,
                     const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_general_struct<struct_banded,struct_dense,V1,V2>::eval(ret, A,p,q, Bc,trans, opts);
    };

    static void eval_rev(Matrix& ret, const M1& A, const matcl::permvec& p, 
                                  const matcl::permvec& q, const M2& B, trans_type trans,
                                  const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_general_struct<struct_banded,struct_dense,V1,V2>::eval_rev(ret, A,p,q,Bc,trans, opts);
    };
};

//--------------------------------------------------------------------------------
//                  BAND-BAND
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_banded,struct_banded,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_banded>;
    using M2    = raw::Matrix<V2,struct_banded>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, trans_type trans,
                     const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_general_struct<struct_banded,struct_dense,V1,V2>::eval(ret, A,p,q, Bc,trans, opts);
    };
    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, trans_type trans,
                         const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_general_struct<struct_banded,struct_dense,V1,V2>::eval_rev(ret, A,p,q, Bc,trans, opts);
    };
};

//--------------------------------------------------------------------------------
//                  SPARSE-BAND
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_sparse,struct_banded,V1,V2>
    : public linsolve_general_sparse<struct_banded, V1, V2>
{};

//--------------------------------------------------------------------------------
//                  SPARSE-SPARSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_general_struct<struct_sparse,struct_sparse,V1,V2>
    : public linsolve_general_sparse<struct_sparse, V1, V2>
{};

}};