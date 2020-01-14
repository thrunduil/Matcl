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

#include "matcl-linalg/decompositions/svd.h"
#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/utils/linalg_utils.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-linalg/decompositions/householder_q.h"
#include "matcl-linalg/linear_eq/linsolve_objects_decomp.h"
#include "matcl-linalg/decompositions/eig/svd_range.h"
#include "matcl-linalg/linear_eq/linsolve.h"

namespace matcl { namespace details
{

template<class V, class S> struct svd_str {};
template<class V, class S> struct make_bidiagonal_str{};

//--------------------------------------------------------------------------------
//                                  DIAG
//--------------------------------------------------------------------------------
struct svd_diag
{
    static void eval(mat_tup_3& ret, const Matrix& A, Integer M, Integer N, bool economy)
    {
        Integer K       = min(M,N);

        Matrix S,I;
        Matrix Ab       = abs(A);
        tie(S,I)        = sort2(Ab,1,false);

        value_code v    = A.get_value_code();
        v               = matrix_traits::real_value_type(v);
        v               = matrix_traits::unify_value_types(v,value_code::v_float);

        Matrix U        = economy == false ? speye(M,M,v) : speye(M,K,v);

        if (economy == false)
        {
            Matrix II   = (mat_col(), I, trans(irange(I.rows() + 1,1,M)));
            U           = U(colon(),colon(II));
        }
        else
        {
            U = U(colon(),colon(I));
        };

        Matrix J        = find(A != 0);

        Matrix V        = ones(N,1,v);
        V(colon(J),1)   = conj(div(A(colon(J)),Ab(colon(J))));

        V               = spdiag(V);

        if (economy == false)
        {
            Matrix II   = (mat_col(),I, trans(irange(I.rows() + 1,1,N)));
            V           = V(colon(),colon(II));
        }
        else
        {
            V           = V(colon(),colon(I));
        };

        if (economy == false)
            S   = bdiags(S,0,M,N);
        else
            S   = bdiags(S,0,K,K);

        if (M == N || economy == false)
        {
            struct_flag sf_u;
            sf_u.set_user(unitary_flag());

            U.add_struct(sf_u);
            V.add_struct(sf_u);
        };

        ret = mat_tup_3(U,S,V);
        return;
    };

    static void eval1(Matrix& ret, const Matrix& A)
    {
        Matrix S    = sort(abs(A),1,false);
        ret         = S;
    };
};

//--------------------------------------------------------------------------------
//                                  HERMITIAN
//--------------------------------------------------------------------------------
struct svd_her
{
    static void eval(mat_tup_3& ret, const Matrix& A, svd_algorithm alg)
    {
        Matrix U, T;
        schur_sym_alg sa= (alg == svd_algorithm::dc)? schur_sym_alg::dc : schur_sym_alg::qr;
        tie(U, T)       = schur(A, sa);
        Matrix S        = T.diag(0);

        Matrix S2,I;
        tie(S2,I)       = sort2(abs(S),1,false);
        U               = U(colon(),colon(I));

        Matrix J        = find(S(I) < 0);
        Matrix V        = U;
        V(colon(), J)   = -V(colon(), J);

        Integer N       = A.rows();
        S               = bdiags(S2,0,N, N);

        struct_flag sf_u;
        sf_u.set_user(unitary_flag());

        U.add_struct(sf_u);
        V.add_struct(sf_u);

        ret = mat_tup_3(U,S,V);
        return;
    };

    static void eval1(Matrix& ret, const Matrix& A)
    {
        Matrix S    = sort(abs(A),1,false);
        ret         = S;
    };
};

template<class V> 
struct svd_str<V,struct_dense> 
{
    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_I = raw::Matrix<Integer,struct_dense>;

    static void eval(mat_tup_3& ret, const Mat& A, bool economy, svd_algorithm alg)
    {
        char jobu   = economy == false ? 'A' : 'S';
        char jobvt  = economy == false ? 'A' : 'S';

        Integer M   = A.rows();
        Integer N   = A.cols();
        Integer K   = min(M, N);
        Integer MK  = economy ? K : M;
        Integer NK  = economy ? K : N;

        using VR    = typename real_type<V>::type;
        using MR    = raw::Matrix<VR,struct_dense>;
        
        if (K == 0)
        {
            value_code v    = matrix_traits::value_code<VR>::value;

            if (economy == false)
            {
                Matrix U    = speye(M,M,v);
                Matrix Vv   = speye(N,N,v);
                Matrix S    = spzeros(M,N,0,v);
                ret = mat_tup_3(U,S,Vv);
            }
            else
            {
                Matrix U    = spzeros(M,K,0,v);
                Matrix Vv   = spzeros(N,K,0,v);
                Matrix S    = spzeros(K,K,0,v);
                ret = mat_tup_3(U,S,Vv);
            };
            return;
        };        

        MR SV(A.get_type(),K,1);
        Mat work(A.get_type(),1,1);
        Mat U(A.get_type(),M,MK);
        Mat VT(A.get_type(),NK,N);

        Mat Ac(A.get_type());

        if (A.get_refstr()->is_unique() == true)
            Ac.assign_to_fresh(A);
        else
            Ac.assign_to_fresh(A.copy());

        using VL = typename lusol_value_type<V>::type;
        using md::lap;

        matcl::Integer info = 0;

        switch (alg)
        {
            case svd_algorithm::qr:
            {                
                lapack::gesvd<VL>(&jobu, &jobvt, M, N, lap(Ac.ptr()), Ac.ld(), lap(SV.ptr()),lap(U.ptr()), U.ld(),
                            lap(VT.ptr()), VT.ld(), lap(work.ptr()), -1, &info);

                Integer lwork = (Integer) real(work.ptr()[0]);
                work.reset_unique(lwork,1);
                lapack::gesvd<VL>(&jobu, &jobvt, M, N, lap(Ac.ptr()), Ac.ld(), SV.ptr(), lap(U.ptr()), U.ld(),
                            lap(VT.ptr()), VT.ld(), lap(work.ptr()), lwork, &info);
                break;
            }
            case svd_algorithm::dc:
            {
                Integer KI = 8*std::min(M,N) + 1;
                Mat_I iwork(A.get_type(),KI,1);

                lapack::gesdd<VL>(&jobu, M, N, lap(Ac.ptr()), Ac.ld(), lap(SV.ptr()),lap(U.ptr()), U.ld(),
                            lap(VT.ptr()), VT.ld(), lap(work.ptr()), -1, iwork.ptr(), &info);

                Integer lwork = (Integer) real(work.ptr()[0]);
                work.reset_unique(lwork,1);
                lapack::gesdd<VL>(&jobu, M, N, lap(Ac.ptr()), Ac.ld(), SV.ptr(), lap(U.ptr()), U.ld(),
                            lap(VT.ptr()), VT.ld(), lap(work.ptr()), lwork, iwork.ptr(), &info);
                break;
            }
        }

        if (info)
            throw error::error_svd();

        Matrix S    = bdiags(Matrix(SV,false), 0, MK, NK);
        Matrix Vv   = ctrans(Matrix(VT,true));

        if (M == N || economy == false)
        {
            struct_flag str_unit;
            str_unit.set_user(unitary_flag());

            U.add_struct(str_unit);
            Vv.add_struct(str_unit);
        };

        ret = mat_tup_3(Matrix(U,true),S,Vv);
        return;
    };

    static void eval1(Matrix& ret, const Mat& A, svd_algorithm alg)
    {
        char jobu   = 'N';
        char jobvt  = 'N';

        Integer M   = A.rows();
        Integer N   = A.cols();
        Integer K   = min(M, N);

        using VR    = typename real_type<V>::type;
        using MR    = raw::Matrix<VR,struct_dense>;
        
        if (K == 0)
        {
            value_code v    = matrix_traits::value_code<VR>::value;
            ret = spzeros(K,1,0,v);
            return;
        };        

        MR SV(A.get_type(),K,1);
        Mat work(A.get_type(),1,1);

        Mat Ac(A.get_type());

        if (A.get_refstr()->is_unique() == true)
            Ac.assign_to_fresh(A);
        else
            Ac.assign_to_fresh(A.copy());

        using VL = typename lusol_value_type<V>::type;
        using md::lap;

        Integer info = 0;

        switch(alg)
        {
            case svd_algorithm::qr:
            {                
                lapack::gesvd<VL>(&jobu, &jobvt, M, N, lap(Ac.ptr()), Ac.ld(), lap(SV.ptr()),0, M, 0, N, 
                            lap(work.ptr()), -1, &info);

                Integer lwork = (Integer) real(work.ptr()[0]);
                work.reset_unique(lwork,1);
                lapack::gesvd<VL>(&jobu, &jobvt, M, N, lap(Ac.ptr()), Ac.ld(), SV.ptr(), 0, 1, 0, 1, 
                            lap(work.ptr()), lwork, &info);
                break;
            }
            case svd_algorithm::dc:
            {
                Integer KI = 8*std::min(M,N) + 1;
                Mat_I iwork(A.get_type(),KI,1);

                lapack::gesdd<VL>(&jobu, M, N, lap(Ac.ptr()), Ac.ld(), lap(SV.ptr()),0, M, 0, N, 
                            lap(work.ptr()), -1, iwork.ptr(), &info);

                Integer lwork = (Integer) real(work.ptr()[0]);
                work.reset_unique(lwork,1);
                lapack::gesdd<VL>(&jobu, M, N, lap(Ac.ptr()), Ac.ld(), SV.ptr(), 0, 1, 0, 1, 
                            lap(work.ptr()), lwork, iwork.ptr(), &info);
                break;
            }
        };
        if (info)
            throw error::error_svd();

        ret = Matrix(SV,true);
        return;
    };
};
template<class V> 
struct make_bidiagonal_str<V,struct_dense> 
{
    using Mat = raw::Matrix<V,struct_dense>;
    static void eval(make_bidiagonal_ret& ret, const Mat& A)
    {
        using VR    = typename md::real_type<V>::type;
        using Mat_R = raw::Matrix<VR,struct_dense>;

        Integer M   = A.rows();
        Integer N   = A.cols();
        Integer K   = std::min(M,N);

        Mat Ac(A.get_type());

        if (A.get_refstr()->is_unique() == true)
            Ac.assign_to_fresh(A);
        else
            Ac.assign_to_fresh(A.copy());

        Mat_R DE(A.get_type(), K, 2);
        Mat TAUQ(A.get_type(), K, 1);
        Mat TAUP(A.get_type(), K, 1);

        Integer info = 0;
        V work_query = V(0);

        VR* ptr_D    = DE.ptr();
        VR* ptr_E    = DE.ptr() + DE.ld();

        using md::lap;

        lapack::gebrd(M, N, lap(Ac.ptr()), Ac.ld(), lap(ptr_D), lap(ptr_E), lap(TAUQ.ptr()), 
                    lap(TAUP.ptr()), lap(&work_query), -1, info );

        if (info)
            throw error::error_general("invalid argument is passed to gebrd");

        Integer lwork = (Integer) real(work_query);

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_w            = reinterpret_cast<V*>(WORK.ptr());

        lapack::gebrd(M, N, lap(Ac.ptr()), Ac.ld(), lap(ptr_D), lap(ptr_E), lap(TAUQ.ptr()), 
                    lap(TAUP.ptr()), lap(ptr_w), lwork, info );

        if (info)
            throw error::error_general("invalid argument is passed to gebrd");

        Integer pos_E   = (M >= N)? 1 : -1;
        Integer off_U   = (M >= N)? 0 : 1;
        Integer off_V   = (M >= N)? 1 : 0;
        Matrix d        = (mat_row(), 0, pos_E);
        Matrix S        = bdiags(matcl::Matrix(DE,false), d, M, N);

        unitary_matrix U, VV;

        using householder_impl = std::shared_ptr<householder_q<V>>;

        bool isv        = Ac.all_finite() && TAUQ.all_finite();
        value_code vc   = matrix_traits::value_code<V>::value;

        if (isv == false)
        {
            U = unitary_matrix::from_nan(Ac.rows(), M, vc);
        }
        else
        {            
            householder_impl U_impl(new details::householder_q<V>(M, Ac, TAUQ, M, off_U));

            U = unitary_matrix(U_impl);
        };

        Matrix V2       = ctrans(Matrix(Ac,false));
        V2              = convert(V2, Mat::matrix_code);
        Mat Ac2         = V2.get_impl_unique<Mat>();

        isv             = Ac2.all_finite() && TAUP.all_finite();

        if (isv == false)
        {
            VV = unitary_matrix::from_nan(Ac2.rows(), N, vc);
        }
        else
        {
            householder_impl V_impl(new details::householder_q<V>(N, Ac2, TAUP, N, off_V));
            VV = unitary_matrix(V_impl);
        };

        ret = make_bidiagonal_ret(U,S,VV);
        return;
    };
};

//--------------------------------------------------------------------------------
//                                  BAND
//--------------------------------------------------------------------------------
template<class V> 
struct svd_str<V,struct_banded> 
{
    using Mat = raw::Matrix<V,struct_banded>;
    static void eval(mat_tup_3& ret, const Mat& A, bool economy, svd_algorithm alg)
    {
        using MC    = raw::Matrix<V,struct_dense>;
        MC AC       = raw::converter<MC,Mat>::eval(A);
        return svd_str<V,struct_dense>::eval(ret, AC, economy,alg);
    };

    static void eval1(Matrix& ret, const Mat& A, svd_algorithm alg)
    {
        //TODO: check bidiagonal
        using MC    = raw::Matrix<V,struct_dense>;
        MC AC       = raw::converter<MC,Mat>::eval(A);
        return svd_str<V,struct_dense>::eval1(ret, AC, alg);
    };
};
template<class V> 
struct make_bidiagonal_str<V,struct_banded> 
{
    using Mat = raw::Matrix<V,struct_banded>;
    static void eval(make_bidiagonal_ret& ret, const Mat& A)
    {
        using MC    = raw::Matrix<V,struct_dense>;
        MC AC       = raw::converter<MC,Mat>::eval(A);
        return make_bidiagonal_str<V,struct_dense>::eval(ret, AC);
    };
};

//--------------------------------------------------------------------------------
//                                  SPARSE
//--------------------------------------------------------------------------------
template<class V> 
struct svd_str<V,struct_sparse> 
{
    using Mat = raw::Matrix<V,struct_sparse>;
    static void eval(mat_tup_3& ret, const Mat& A, bool economy, svd_algorithm alg)
    {
        using MC    = raw::Matrix<V,struct_dense>;
        MC AC       = raw::converter<MC,Mat>::eval(A);
        return svd_str<V,struct_dense>::eval(ret, AC, economy, alg);
    };

    static void eval1(Matrix& ret, const Mat& A, svd_algorithm alg)
    {
        using MC    = raw::Matrix<V,struct_dense>;
        MC AC       = raw::converter<MC,Mat>::eval(A);
        return svd_str<V,struct_dense>::eval1(ret, AC, alg);
    };
};
template<class V> 
struct make_bidiagonal_str<V,struct_sparse> 
{
    using Mat = raw::Matrix<V,struct_sparse>;
    static void eval(make_bidiagonal_ret& ret, const Mat& A)
    {
        using MC    = raw::Matrix<V,struct_dense>;
        MC AC       = raw::converter<MC,Mat>::eval(A);
        return make_bidiagonal_str<V,struct_dense>::eval(ret, AC);
    };
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct svd_impl
{
    using M = raw::Matrix<V,S>;

    static void eval(mat_tup_3& ret, const M& A, bool economy, svd_algorithm alg)
    {
        return svd_str<V,S>::eval(ret,A, economy, alg);
    };
    static void eval1(Matrix& ret, const M& A, svd_algorithm alg)
    {
        return svd_str<V,S>::eval1(ret,A, alg);
    };
};

template<class S>
struct svd_impl<Integer,S>
{
    using M = raw::Matrix<Integer,S>;
    static void eval(mat_tup_3& ret, const M& A, bool economy, svd_algorithm alg)
    {
        using MC    = raw::Matrix<Real,struct_dense>;
        MC AC       = raw::converter<MC,M>::eval(A);
        return svd_impl<Real,struct_dense>::eval(ret,AC,economy,alg);
    };
    static void eval1(Matrix& ret, const M& A, svd_algorithm alg)
    {
        using MC    = raw::Matrix<Real,struct_dense>;
        MC AC       = raw::converter<MC,M>::eval(A);
        return svd_impl<Real,struct_dense>::eval1(ret,AC,alg);
    };
};
template<class S>
struct svd_impl<Object,S>
{
    using M = raw::Matrix<Object,S>;
    static void eval(mat_tup_3&, const M&, bool, svd_algorithm)
    {
        throw error::object_value_type_not_allowed("svd");
    };
    static void eval1(Matrix&, const M&, svd_algorithm)
    {
        throw error::object_value_type_not_allowed("svd");
    };
};

template<class V>
struct gsvd_impl
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_R = raw::Matrix<VR,struct_dense>;

    static void eval(gsvd_ret& ret, const Matrix& A, const Matrix& B)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();
        Integer P           = B.rows();

        if (A.all_finite() == false || B.all_finite() == false)
        {
            Matrix mat_U    = make_nan_matrix<V>(M,M);
            Matrix mat_V    = make_nan_matrix<V>(P,P);
            Matrix mat_Q    = make_nan_matrix<V>(N,N);
            Integer K       = 0;
            Integer L       = 0;
            Matrix R        = make_nan_matrix<V>(0,0);
            Matrix D1       = make_nan_matrix<VR>(M,0);
            Matrix D2       = make_nan_matrix<VR>(P,0);
            ret = gsvd_ret(mat_U, mat_V, mat_Q, D1, D2, R, K, L);
            return;
        };

        bool is_complex     = md::is_complex<V>::value;

        const char* JOBU    = "U";
        const char* JOBV    = "V";
        const char* JOBQ    = "Q";
        Integer K, L;

        Mat mat_A           = convert(A, Mat::matrix_code).get_impl<Mat>().make_unique();
        Mat mat_B           = convert(B, Mat::matrix_code).get_impl<Mat>().make_unique();

        V* ptr_A            = mat_A.ptr();
        V* ptr_B            = mat_B.ptr();
        Integer LDA         = mat_A.ld();
        Integer LDB         = mat_B.ld();

        Mat_R alpha(ti::ti_empty(), N, 1);
        Mat_R beta(ti::ti_empty(), N, 1);
        Mat U(ti::ti_empty(), M, M);
        Mat VV(ti::ti_empty(), P, P);
        Mat Q(ti::ti_empty(), N, N);

        VR* ptr_alpha       = alpha.ptr();
        VR* ptr_beta        = beta.ptr();
        V* ptr_U            = U.ptr();
        V* ptr_V            = VV.ptr();
        V* ptr_Q            = Q.ptr();
        Integer LDU         = U.ld();
        Integer LDV         = VV.ld();
        Integer LDQ         = Q.ld();

        Integer liwork      = N;
        Integer lrwork      = is_complex ? 2 * N : 0;
        using iworkspace    = matcl::pod_workspace<Integer>;
        using rworkspace    = matcl::pod_workspace<VR>;
        iworkspace iwork    = iworkspace(liwork);
        rworkspace rwork    = rworkspace(lrwork);
        Integer* ptr_iwork  = iwork.ptr();
        VR* ptr_rwork       = rwork.ptr();

        Integer info;

        V work_query;
        lapack::ggsvd3( JOBU, JOBV, JOBQ, M, N, P, K, L, lap(ptr_A), LDA, lap(ptr_B), LDB, lap(ptr_alpha), 
                lap(ptr_beta), lap(ptr_U), LDU, lap(ptr_V), LDV, lap(ptr_Q), LDQ, lap(&work_query), -1, 
                ptr_rwork, ptr_iwork, info);

        Integer lwork       = (Integer)real(work_query);

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace work      = workspace(lwork);
        V* ptr_work         = (V*)work.ptr();

        lapack::ggsvd3( JOBU, JOBV, JOBQ, M, N, P, K, L, lap(ptr_A), LDA, lap(ptr_B), LDB, lap(ptr_alpha), 
                lap(ptr_beta), lap(ptr_U), LDU, lap(ptr_V), LDV, lap(ptr_Q), LDQ, lap(ptr_work), lwork, 
                ptr_rwork, ptr_iwork, info);

        if (info < 0)
            throw error::error_general("invalid argument passed to ggsvd3");
        if (info != 0)
            throw error::gsvd_failed();

        Matrix mat_U(U, true);
        Matrix mat_V(VV, true);
        Matrix mat_Q(Q, true);

        Matrix D1, D2, R;

        if (M - K - L >= 0)
        {
            D1  = bdiags(Matrix(alpha, false)(colon(1, K + L)), 0, M, K + L);
            D2  = bdiags(Matrix(beta, false)(colon(K + 1, K + L)), K, P, K + L);
            R   = triu(Matrix(mat_A, false), N - K - L);
        }
        else
        {
            D1  = bdiags(Matrix(alpha, false)(colon(1, M)), 0, M, K + L);
            D2  = bdiags(Matrix(beta, false)(colon(K + 1, K+L)), K, P, K + L);
            Matrix AR   = Matrix(mat_A, false);
            Matrix BR   = Matrix(mat_B, false)(colon(M - K + 1, L), colon());
            R   = mat_col().add(AR).add(BR);
            R   = triu(R, N - K - L);
        };

        struct_flag sf_u;
        sf_u.set_user(unitary_flag());

        mat_U.add_struct(sf_u);
        mat_V.add_struct(sf_u);
        mat_Q.add_struct(sf_u);

        ret = gsvd_ret(mat_U, mat_V, mat_Q, D1, D2, R, K, L);
        return;
    }

    static void eval1(Matrix& ret, const Matrix& A, const Matrix& B)
    {
        Integer M           = A.rows();
        Integer N           = A.cols();
        Integer P           = B.rows();

        if (A.all_finite() == false || B.all_finite() == false)
        {
            Matrix E        = make_nan_matrix<V>(N,1);
            ret = E;
            return;
        };

        bool is_complex     = md::is_complex<V>::value;

        const char* JOBU    = "N";
        const char* JOBV    = "N";
        const char* JOBQ    = "N";
        Integer K, L;

        Mat mat_A           = convert(A, Mat::matrix_code).get_impl<Mat>().make_unique();
        Mat mat_B           = convert(B, Mat::matrix_code).get_impl<Mat>().make_unique();

        V* ptr_A            = mat_A.ptr();
        V* ptr_B            = mat_B.ptr();
        Integer LDA         = mat_A.ld();
        Integer LDB         = mat_B.ld();

        Mat_R alpha(ti::ti_empty(), N, 1);
        Mat_R beta(ti::ti_empty(), N, 1);

        VR* ptr_alpha       = alpha.ptr();
        VR* ptr_beta        = beta.ptr();
        V* ptr_U            = nullptr;
        V* ptr_V            = nullptr;
        V* ptr_Q            = nullptr;
        Integer LDU         = 1;
        Integer LDV         = 1;
        Integer LDQ         = 1;

        Integer liwork      = N;
        Integer lrwork      = is_complex ? 2 * N : 0;
        using iworkspace    = matcl::pod_workspace<Integer>;
        using rworkspace    = matcl::pod_workspace<VR>;
        iworkspace iwork    = iworkspace(liwork);
        rworkspace rwork    = rworkspace(lrwork);
        Integer* ptr_iwork  = iwork.ptr();
        VR* ptr_rwork       = rwork.ptr();

        Integer info;

        V work_query;
        lapack::ggsvd3( JOBU, JOBV, JOBQ, M, N, P, K, L, lap(ptr_A), LDA, lap(ptr_B), LDB, lap(ptr_alpha), 
                lap(ptr_beta), lap(ptr_U), LDU, lap(ptr_V), LDV, lap(ptr_Q), LDQ, lap(&work_query), -1, 
                ptr_rwork, ptr_iwork, info);

        Integer lwork       = (Integer)real(work_query);

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace work      = workspace(lwork);
        V* ptr_work         = (V*)work.ptr();

        lapack::ggsvd3( JOBU, JOBV, JOBQ, M, N, P, K, L, lap(ptr_A), LDA, lap(ptr_B), LDB, lap(ptr_alpha), 
                lap(ptr_beta), lap(ptr_U), LDU, lap(ptr_V), LDV, lap(ptr_Q), LDQ, lap(ptr_work), lwork, 
                ptr_rwork, ptr_iwork, info);

        if (info < 0)
            throw error::error_general("invalid argument passed to ggsvd3");
        if (info != 0)
            throw error::gsvd_failed();

        for (Integer i = 0; i < N; ++i)
            ptr_alpha[i]    = ptr_alpha[i] / ptr_beta[i];


        ret = Matrix(alpha, false);
        return;
    }

};

template<class V, class S>
struct make_bidiagonal_impl
{
    using M = raw::Matrix<V,S>;

    static void eval(make_bidiagonal_ret& ret, const M& A)
    {
        return make_bidiagonal_str<V,S>::eval(ret,A);
    };
};

template<class S>
struct make_bidiagonal_impl<Integer,S>
{
    using M = raw::Matrix<Integer,S>;
    static void eval(make_bidiagonal_ret& ret, const M& A)
    {
        using MC    = raw::Matrix<Real,struct_dense>;
        MC AC       = raw::converter<MC,M>::eval(A);
        return make_bidiagonal_impl<Real,struct_dense>::eval(ret,AC);
    };
};


template<class Compl>
static void svd_scalar(mat_tup_3& ret, const Compl& mat)
{
    using constants::nan;
    using VR    = typename md::real_type_int_real<Complex>::type;

    bool isv = is_finite(mat);

    if (isv == false)
    {
        ret = mat_tup_3(nan<VR>(), nan<VR>(), nan<VR>());
        return;
    };

    Compl S     = abs(mat);
    ret         = mat_tup_3(mat/S,S,VR(1.));
    return;
};

struct unary_visitor_svd : public extract_type_switch<void,unary_visitor_svd,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, mat_tup_3& ret, bool economy, svd_algorithm alg)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;        
        using VR    = typename md::real_type_int_real<V>::type;
     
        bool isv    = mat.all_finite();
        
        if (isv == false)
        {
            Integer K   = min(mat.rows(),mat.cols());

            if (economy == false)
            {
                ret = mat_tup_3(md::make_nan_matrix<VR>(mat.rows(), mat.rows()),
                                md::make_nan_matrix<VR>(mat.rows(), mat.cols()),
                                md::make_nan_matrix<VR>(mat.cols(), mat.cols()));
            }
            else
            {
                ret = mat_tup_3(md::make_nan_matrix<VR>(mat.rows(), K),
                                md::make_nan_matrix<VR>(K, K),
                                md::make_nan_matrix<VR>(mat.cols(), K));
            }

            return;
        };

        return details::svd_impl<V,S>::eval(ret, mat, economy, alg);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, mat_tup_3&, bool, svd_algorithm)
    {
        throw error::object_value_type_not_allowed("svd");
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, mat_tup_3& ret, bool, svd_algorithm)
    {
        using constants::nan;

        using VR    = typename md::real_type_int_real<T>::type;

        bool isv = is_finite(mat);

        if (isv == false)
        {
            ret = mat_tup_3(nan<VR>(), nan<VR>(), nan<VR>());
            return;
        }

        if (mat < 0)
            ret = mat_tup_3(VR(1.),abs(mat), VR(-1.));
        else
            ret = mat_tup_3(VR(1.),abs(mat),VR(1.));

        return;
    };

    static void eval_scalar(const Matrix&, const Complex& mat, mat_tup_3& ret, bool, svd_algorithm)
    {
        return svd_scalar<Complex>(ret, mat);
    };
    static void eval_scalar(const Matrix&, const Float_complex& mat, mat_tup_3& ret, bool, svd_algorithm)
    {
        return svd_scalar<Float_complex>(ret, mat);
    };
    static void eval_scalar(const Matrix&, const Object&, mat_tup_3&, bool, svd_algorithm)
    {
        throw error::object_value_type_not_allowed("svd");
    };
};

struct unary_visitor_svd1 : public extract_type_switch<void,unary_visitor_svd1,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, svd_algorithm alg)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;        
        using VR    = typename md::real_type_int_real<V>::type;

        bool isv    = mat.all_finite();
        
        if (isv == false)
        {
            ret = details::make_nan_matrix<VR>(min(mat.rows(), mat.cols()), 1);
            return;
        };

        return details::svd_impl<V,S>::eval1(ret, mat, alg);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, svd_algorithm)
    {
        throw error::object_value_type_not_allowed("svd");
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, Matrix& ret, svd_algorithm)
    {
        using constants::nan;

        using VR    = typename md::real_type_int_real<T>::type;

        bool isv = is_finite(mat);
    
        if (isv == false)
        {
            ret = nan<VR>();
            return;
        }

        ret = abs(mat);
        return;
    };

    static void eval_scalar(const Matrix&, const Object&, Matrix&, svd_algorithm)
    {
        throw error::object_value_type_not_allowed("svd");
    };
};

struct unary_visitor_gsvd : public extract_type_switch<void,unary_visitor_gsvd,true>
{
    template<class T>
    static void eval(const Matrix&, const T&, gsvd_ret& ret, const Matrix& A, const Matrix& B)
    {
        using V0    = typename T::value_type;
        using V     = typename md::unify_types<V0, Float>::type;
        return details::gsvd_impl<V>::eval(ret, A, B);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, gsvd_ret&, const Matrix&, const Matrix&)
    {
        throw error::object_value_type_not_allowed("gsvd");
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, gsvd_ret& ret, const Matrix& A, const Matrix& B)
    {
        using V     = typename md::unify_types<T, Float>::type;
        return gsvd_impl<V>::eval(ret, A, B);
    };

    static void eval_scalar(const Matrix&, const Object&, gsvd_ret&, const Matrix&, const Matrix&)
    {
        throw error::object_value_type_not_allowed("gsvd");
    };
};

struct unary_visitor_gsvd1 : public extract_type_switch<void,unary_visitor_gsvd1,true>
{
    template<class T>
    static void eval(const Matrix&, const T&, Matrix& ret, const Matrix& A, const Matrix& B)
    {
        using V0    = typename T::value_type;
        using V     = typename md::unify_types<V0, Float>::type;
        return details::gsvd_impl<V>::eval1(ret, A, B);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const Matrix&, const Matrix&)
    {
        throw error::object_value_type_not_allowed("gsvd");
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T&, Matrix& ret, const Matrix& A, const Matrix& B)
    {
        using V     = typename md::unify_types<T, Float>::type;
        return gsvd_impl<V>::eval1(ret, A, B);
    };

    static void eval_scalar(const Matrix&, const Object&, Matrix&, const Matrix&, const Matrix&)
    {
        throw error::object_value_type_not_allowed("gsvd1");
    };
};

struct make_bidiagonal_vis : public extract_type_switch<void,make_bidiagonal_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, make_bidiagonal_ret& ret)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;        
        using VR    = typename md::real_type_int_real<V>::type;
     
        bool isv    = mat.all_finite();
        
        if (isv == false)
        {
            Integer M           = mat.rows();
            Integer N           = mat.cols();
            value_code vc       = matrix_traits::value_code<VR>::value;
            unitary_matrix U    = unitary_matrix::from_nan(M, M, vc);
            unitary_matrix V    = unitary_matrix::from_nan(N, N, vc);

            ret = make_bidiagonal_ret(U,md::make_nan_matrix<VR>(M, N),V);
            return;
        };

        return details::make_bidiagonal_impl<V,S>::eval(ret, mat);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, make_bidiagonal_ret&)
    {
        throw error::object_value_type_not_allowed("make_bidiagonal");
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, make_bidiagonal_ret& ret)
    {
        using VR    = typename md::real_type_int_real<T>::type;

        unitary_matrix U(VR(1.0),false);
        ret         = make_bidiagonal_ret(U, mat, U);
        return;
    };

    static void eval_scalar(const Matrix&, const Object&, make_bidiagonal_ret&)
    {
        throw error::object_value_type_not_allowed("make_bidiagonal");
    };
};

};};

namespace matcl
{

static void svd_impl(mat_tup_3& ret, const Matrix& A, bool economy, svd_algorithm alg)
{
    Integer M   = A.rows();
    Integer N   = A.cols();
    Integer K   = min(M,N);

    value_code v    = A.get_value_code();
    bool is_real    = matrix_traits::is_float_real(v);
    v               = matrix_traits::real_value_type(v);
    v               = matrix_traits::unify_value_types(v,value_code::v_float);

    if (A.structural_nnz() == 0)
    {
        Matrix U    = economy == false ? speye(M,M,v) : speye(M,K,v);
        Matrix S    = economy == false ? spzeros(M,N,0,v) : spzeros(K,K,0,v);
        Matrix V    = economy == false ? speye(N,N,v) : speye(N,K,v);

        ret = mat_tup_3(U,S,V);
        return;
    };

    if (A.get_struct().is_id())
    {
        Matrix U    = economy == false ? speye(M,M,v) : speye(M,K,v);
        Matrix S    = economy == false ? speye(M,N,v) : speye(K,K,v);
        Matrix V    = economy == false ? speye(N,N,v) : speye(N,K,v);

        ret = mat_tup_3(U,S,V);
        return;
    }

    if (is_diag(A))
        return details::svd_diag::eval(ret, full(get_diag(A)), A.rows(), A.cols(), economy);

    if (is_unitary(A.get_struct()) && M == N)
    {
        Matrix S    = speye(M,N,v);
        Matrix V    = speye(N,N,v);

        ret = mat_tup_3(A,S,V);
        return;
    }

    if (A.get_struct().is_hermitian(A.is_square(), is_real))
        return details::svd_her::eval(ret, A, alg);

    return details::unary_visitor_svd::make<const Matrix&>(A,ret, economy,alg);
};

static void gsvd_impl(gsvd_ret& ret, const Matrix& A, const Matrix& B)
{
    Integer NA      = A.cols();
    Integer MB      = B.rows();
    Integer NB      = B.cols();

    if (NA != NB)
        throw error::invalid_size2(MB, NB, MB, NA);

    value_code vA   = A.get_value_code();
    value_code vB   = B.get_value_code();
    value_code v    = matrix_traits::unify_value_types(vA,vB);
    v               = matrix_traits::unify_value_types(v,value_code::v_float);

    Matrix tmp      = zeros(0,0, v);
    return details::unary_visitor_gsvd::make<const Matrix&>(tmp,ret, A, B);
};

static void gsvd1_impl(Matrix& ret, const Matrix& A, const Matrix& B)
{
    Integer NA      = A.cols();
    Integer MB      = B.rows();
    Integer NB      = B.cols();

    if (NA != NB)
        throw error::invalid_size2(MB, NB, MB, NA);

    value_code vA   = A.get_value_code();
    value_code vB   = B.get_value_code();
    value_code v    = matrix_traits::unify_value_types(vA,vB);
    v               = matrix_traits::unify_value_types(v,value_code::v_float);

    Matrix tmp      = zeros(0,0, v);
    return details::unary_visitor_gsvd1::make<const Matrix&>(tmp,ret, A, B);
};

static void svd1_impl(Matrix& ret, const Matrix& A, svd_algorithm alg)
{
    value_code v    = A.get_value_code();
    bool is_real    = matrix_traits::is_float_real(v);
    v               = matrix_traits::real_value_type(v);
    v               = matrix_traits::unify_value_types(v,value_code::v_float);

    if (A.structural_nnz() == 0)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();
        Integer K   = min(M,N);

        Matrix S = spzeros(K,1,0,v);
        ret = S;
        return;
    };

    if (A.get_struct().is_id())
    {
        Integer M   = A.rows();
        Integer N   = A.cols();
        Integer K   = min(M,N);

        Matrix S    = spones(K,1,v);
        ret = S;
        return;
    }

    if (is_diag(A))
        return details::svd_diag::eval1(ret, full(get_diag(A)));

    if (A.get_struct().is_hermitian(A.is_square(), is_real))
    {
        Matrix E    = abs(eig(A));
        ret         = sort(E, 1, false);
        return;
    };

    return details::unary_visitor_svd1::make<const Matrix&>(A,ret,alg);
};

static void make_bidiagonal_impl(make_bidiagonal_ret& ret, const Matrix& A)
{
    Integer M       = A.rows();
    Integer N       = A.cols();

    value_code v    = A.get_value_code();
    v               = matrix_traits::real_value_type(v);
    v               = matrix_traits::unify_value_types(v,value_code::v_float);

    Integer ld      = get_ld(A,1);
    Integer ud      = get_ud(A,1);

    if (A.structural_nnz() == 0 || ld == 0 && ud == 0
        || (M >= N && ld == 0 && ud <= 1) || (N < N && ld <= 1 && ud == 0))
    {
        Matrix U    = speye(M,M,v);
        Matrix S    = A;
        Matrix V    = speye(N,N,v);

        ret = make_bidiagonal_ret(unitary_matrix(U,false),S,unitary_matrix(V,false));
        return;
    };

    return details::make_bidiagonal_vis::make<const Matrix&>(A,ret);
};

mat_tup_3 matcl::svd(const Matrix& A0, bool economy, svd_algorithm alg)
{
    //increase refcount
    Matrix A(A0);

    mat_tup_3 ret;
    svd_impl(ret, A, economy, alg);
    return ret;
};
mat_tup_3 matcl::svd(Matrix&& A0, bool economy, svd_algorithm alg)
{
    Matrix A(std::move(A0));

    mat_tup_3 ret;
    svd_impl(ret, A, economy, alg);
    return ret;
};
Matrix matcl::svd1(const Matrix& A0, svd_algorithm alg)
{
    //increase refcount
    Matrix A(A0);

    Matrix ret;
    svd1_impl(ret, A, alg);
    return ret;
};
Matrix matcl::svd1(Matrix&& A0, svd_algorithm alg)
{
    Matrix A(std::move(A0));

    Matrix ret;
    svd1_impl(ret, A, alg);
    return ret;
};

gsvd_ret matcl::gsvd(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    gsvd_ret ret;
    gsvd_impl(ret, A, B);
    return ret;
};
gsvd_ret matcl::gsvd(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    gsvd_ret ret;
    gsvd_impl(ret, A, B);
    return ret;
};
gsvd_ret matcl::gsvd(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    gsvd_ret ret;
    gsvd_impl(ret, A, B);
    return ret;
};
gsvd_ret matcl::gsvd(Matrix&& A0, Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    gsvd_ret ret;
    gsvd_impl(ret, A, B);
    return ret;
};

Matrix matcl::gsvd1(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    Matrix ret;
    gsvd1_impl(ret, A, B);
    return ret;
};
Matrix matcl::gsvd1(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    Matrix ret;
    gsvd1_impl(ret, A, B);
    return ret;
};
Matrix matcl::gsvd1(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    Matrix ret;
    gsvd1_impl(ret, A, B);
    return ret;
};
Matrix matcl::gsvd1(Matrix&& A0, Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    Matrix ret;
    gsvd1_impl(ret, A, B);
    return ret;
};

make_bidiagonal_ret matcl::make_bidiagonal(const Matrix& A0)
{
    Matrix A(A0);

    make_bidiagonal_ret ret;
    make_bidiagonal_impl(ret, A);
    return ret;
};
make_bidiagonal_ret matcl::make_bidiagonal(Matrix&& A0)
{
    Matrix A(std::move(A0));

    make_bidiagonal_ret ret;
    make_bidiagonal_impl(ret, A);
    return ret;
};

linsolve_obj matcl::linsolve_svd(const Matrix& A, const Matrix& U, const Matrix& S, const Matrix& V,
                                 const options& opts)
{
    return linsolve_svd(A, unitary_matrix(U,true), S, unitary_matrix(V,true), opts);
};
linsolve_obj matcl::linsolve_svd(const Matrix& A, const unitary_matrix& U, const Matrix& S, 
                                const unitary_matrix& V, const options& opts)
{
    Integer N   = U.rows();

    if (U.cols() != S.rows())
        throw error::invalid_svd_factors();
    if (V.cols() != S.cols())
        throw error::invalid_svd_factors();
    if (U.rows() < U.cols())
        throw error::invalid_svd_factors();
    if (V.rows() < V.cols())
        throw error::invalid_svd_factors();
    if (A.rows() != U.rows())
        throw error::invalid_svd_factors();
    if (A.cols() != V.rows())
        throw error::invalid_svd_factors();

    if (matcl::is_diag(S) == false)
        throw error::invalid_svd_factors();

    if (S.rows() != S.cols())
        throw error::square_matrix_required(S.rows(), S.cols());
    if (U.rows() != U.cols())
        throw error::square_matrix_required(U.rows(), U.cols());
    if (V.rows() != V.cols())
        throw error::square_matrix_required(V.rows(), V.cols());

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_svd");
    if (U.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_svd");
    if (S.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_svd");
    if (V.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_svd");

    if (U.rows() == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(S.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(V.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, U.get_type())));
    };

    bool isv            = S.all_finite() && U.all_finite() == true && V.all_finite() == true
                        && A.all_finite() == true;

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(S.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(V.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, S.get_type() )));
    };

    using data_ptr  = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_svd(A, U, S, V, opts)));
};

linsolve_obj matcl::linsolve_bidiagonal(const Matrix& A, const Matrix& U, const Matrix& R, const Matrix& V,
                                        const options& opts)
{
    return linsolve_bidiagonal(A, unitary_matrix(U,true), R, unitary_matrix(V,true), opts);
};
linsolve_obj matcl::linsolve_bidiagonal(const Matrix& A, const unitary_matrix& U, const Matrix& R, 
                                        const unitary_matrix& V, const options& opts)
{
    Integer N   = U.rows();

    if (A.rows() != U.rows())
        throw error::invalid_bidiag_factors();
    if (U.cols() != R.rows())
        throw error::invalid_bidiag_factors();
    if (V.cols() != R.cols())
        throw error::invalid_bidiag_factors();
    if (U.rows() < U.cols())
        throw error::invalid_bidiag_factors();
    if (V.rows() < V.cols())
        throw error::invalid_bidiag_factors();
    if (A.cols() != V.rows())
        throw error::invalid_bidiag_factors();

    Integer ld  = matcl::get_ld(R, 1);
    Integer ud  = matcl::get_ud(R, 1);

    //test only if R is bidiagonal
    if ((ld <= 1  && ud == 0) == false && (ld == 0  && ud <= 1) == false)
        throw error::invalid_bidiag_factors();

    if (R.rows() != R.cols())
        throw error::square_matrix_required(R.rows(), R.cols());
    if (U.rows() != U.cols())
        throw error::square_matrix_required(U.rows(), U.cols());
    if (V.rows() != V.cols())
        throw error::square_matrix_required(V.rows(), V.cols());

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_bigiag");
    if (U.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_bigiag");
    if (R.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_bigiag");
    if (V.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_bigiag");

    if (U.rows() == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(R.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(V.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, U.get_type())));
    };

    bool isv            = R.all_finite() && U.all_finite() == true && V.all_finite() == true
                        && A.all_finite() == true;

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(R.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(V.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, R.get_type() )));
    };

    using data_ptr  = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_bidiag(A, U, R, V, opts)));
};

struct linsolve_svd_impl
{
    static linsolve_obj eval(const Matrix& A, const options& opts)
    {
        if (A.rows() != A.cols())
            throw error::square_matrix_required(A.rows(), A.cols());        

        if (A.rows() == 0)
        {
            value_code vc   = matrix_traits::unify_value_types(A.get_value_code(), value_code::v_float);

            using data_ptr = linsolve_obj::linsolve_data_ptr;
            return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, A.get_type())));
        };

        Matrix U, S, V;     
        tie(U,S,V)  = svd(A, true);

        return linsolve_svd(A,U,S,V, opts);
    }
};

linsolve_obj matcl::linsolve_svd(const Matrix& A, const options& opts)
{
    return linsolve_svd_impl::eval(A, opts);
};
linsolve_obj matcl::linsolve_svd(Matrix&& A, const options& opts)
{
    return linsolve_svd_impl::eval(std::move(A), opts);
};

Matrix matcl::null(const Matrix& A, bool right, Real tol)
{
    operator_spaces os(A);
    if (right == true)
        return os.null_right(tol);
    else
        return os.null_left(tol);
};

Matrix matcl::null(Matrix&& A, bool right, Real tol)
{
    operator_spaces os(std::move(A));
    if (right == true)
        return os.null_right(tol);
    else
        return os.null_left(tol);
};

Matrix matcl::orth(const Matrix& A, bool right, Real tol)
{
    operator_spaces os(A);
    if (right == true)
        return os.orth_right(tol);
    else
        return os.orth_left(tol);
};

Matrix matcl::orth(Matrix&& A, bool right, Real tol)
{
    operator_spaces os(std::move(A));
    if (right == true)
        return os.orth_right(tol);
    else
        return os.orth_left(tol);
};

Integer matcl::rank(const Matrix& A, Real tol)
{
    operator_spaces os(A);
    return os.rank(tol);
};

Integer matcl::rank(Matrix&& A, Real tol)
{
    operator_spaces os(std::move(A));
    return os.rank(tol);
};

#if 0
Matrix matcl::singsel_range(const Matrix& A0, Real VL, Real VU)
{
    Matrix A(A0);
    Matrix E;
    svd_range::eval_range(E, A, VL, VU);
    return E;
}
Matrix matcl::singsel_range(Matrix&& A0, Real VL, Real VU)
{
    Matrix A(std::move(A0));

    Matrix E;
    svd_range::eval_range(E, A, VL, VU);
    return E;
}

Matrix matcl::singsel_index(const Matrix& A0, Integer IF, Integer IL)
{
    Matrix A(A0);

    Matrix E;
    svd_range::eval_index(E, A, IF, IL);
    return E;
}
Matrix matcl::singsel_index(Matrix&& A0, Integer IF, Integer IL)
{
    Matrix A(std::move(A0));

    Matrix E;
    svd_range::eval_index(E, A, IF, IL);
    return E;
}

mat_tup_3 matcl::singsel_range2(const Matrix& A0, Real VL, Real VU)
{
    Matrix A(A0);

    mat_tup_3 ret;
    svd_range::eval2_range(ret, A, VL, VU);
    return ret;
}
mat_tup_3 matcl::singsel_range2(Matrix&& A0, Real VL, Real VU)
{
    Matrix A(std::move(A0));

    mat_tup_3 ret;
    svd_range::eval2_range(ret, A, VL, VU);
    return ret;
}

mat_tup_3 matcl::singsel_index2(const Matrix& A0, Integer IF, Integer IL)
{
    Matrix A(A0);

    mat_tup_3 ret;
    svd_range::eval2_index(ret, A, IF, IL);
    return ret;
}
mat_tup_3 matcl::singsel_index2(Matrix&& A0, Integer IF, Integer IL)
{
    Matrix A(std::move(A0));

    mat_tup_3 ret;
    svd_range::eval2_index(ret, A, IF, IL);
    return ret;
}
#endif

//--------------------------------------------------------------------------------
//                                  operator_spaces
//--------------------------------------------------------------------------------

operator_spaces::operator_spaces()
    :operator_spaces(Matrix(1.0))
{};
operator_spaces::operator_spaces(const Matrix& A)
{
    tie(m_U, m_S, m_V)  = svd(A, false);
    m_S                 = m_S.diag(0);
};

operator_spaces::operator_spaces(Matrix&& A)
{
    tie(m_U, m_S, m_V)  = svd(std::move(A), false);
    m_S                 = m_S.diag(0);
};

operator_spaces::~operator_spaces()
{};

Real operator_spaces::default_tol() const
{
    if (m_S.length() == 0)
        return 0.0;

    Integer K   = std::max(m_U.rows(), m_U.cols());
    Real tol    = m_S(1).get_scalar<Real>() * constants::eps(m_S.get_value_code()) * K;
    return tol;
};

Matrix operator_spaces::null_right(Real tol) const
{
    if (tol < 0)
        tol = default_tol();

    Matrix I    = find(m_S <= tol);
    Matrix N    = m_V(colon(), I);
    return N;
};
Matrix operator_spaces::null_left(Real tol) const
{
    if (tol < 0)
        tol = default_tol();

    Matrix I    = find(m_S <= tol);
    Matrix N    = m_U(colon(), I);
    return ctrans(N);
};
Matrix operator_spaces::orth_right(Real tol) const
{
    if (tol < 0)
        tol = default_tol();

    Matrix I    = find(m_S > tol);
    Matrix N    = m_V(colon(), I);
    return N;
};
Matrix operator_spaces::orth_left(Real tol) const
{
    if (tol < 0)
        tol = default_tol();

    Matrix I    = find(m_S > tol);
    Matrix N    = m_U(colon(), I);
    return ctrans(N);
};

Integer operator_spaces::rank(Real tol) const
{
    if (tol < 0)
        tol = default_tol();

    Matrix I    = find(m_S > tol);
    Integer ret = I.length();
    return ret;
};

};
