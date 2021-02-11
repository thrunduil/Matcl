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

#include "matcl-linalg/decompositions/schur.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/decompositions/eig/schur_utils.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-matrep/details/struct_flag_predefined.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/linear_eq/linsolve_objects_decomp.h"
#include "matcl-linalg/linear_eq/linsolve.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details 
{

template<class V, class S, bool Is_compl = md::is_float_real_scalar<V>::value>
struct schur_is_sym_her {};

template<class V, class S>
struct schur_is_sym_her<V,S,true>
{
    using Mat = raw::Matrix<V, S>;
    static bool eval(const Mat &A)
    {
        return is_sym(Matrix(A, false), 0.0);
    }
};
template<class V, class S>
struct schur_is_sym_her<V,S,false>
{
    using Mat = raw::Matrix<V, S>;
    static bool eval(const Mat &A)
    {
        return is_her(Matrix(A, false), 0.0);
    }
};

template<class V, bool Is_real = md::is_float_real_scalar<V>::value>
struct eval_syevd_heevd
{
    static void eval(const char* job_U, const char* uplo, Integer N, V* TA, Integer TA_ld,
                     V* eig, Integer& info)
    {
        V       work_query;
        Integer iwork_query;

        lapack::syevd(job_U, uplo, N, lap(TA), TA_ld, lap(eig), lap(&work_query), -1, 
                      lap(&iwork_query), -1, &info);
    
        Integer lwork   = (Integer) real(work_query);
        Integer liwork  = iwork_query;

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

        using iworkspace    = matcl::pod_workspace<Integer>;
        iworkspace IWORK    = iworkspace(liwork);
        Integer* ptr_IWORK  = IWORK.ptr();

        lapack::syevd(job_U, uplo, N, lap(TA), TA_ld, lap(eig), lap(ptr_WORK), lwork, ptr_IWORK,
                      liwork, &info);
    }
};
template<class V>
struct eval_syevd_heevd<V,false>
{
    using VR = typename md::real_type<V>::type;
    static void eval(const char* job_U, const char* uplo, Integer N, V* TA, Integer TA_ld,
                     VR* eig, Integer& info)
    {
        V       work_query;
        Integer iwork_query;
        VR      rwork_query;

        lapack::heevd(job_U, uplo, N, lap(TA), TA_ld, lap(eig), lap(&work_query), -1,
                      lap(&rwork_query), -1, lap(&iwork_query), -1, &info);
    
        Integer lwork   = (Integer) real(work_query);
        Integer liwork  = iwork_query;
        Integer lrwork  = (Integer)rwork_query;

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

        using iworkspace    = matcl::pod_workspace<Integer>;
        iworkspace IWORK    = iworkspace(liwork);
        Integer* ptr_IWORK  = IWORK.ptr();

        using rworkspace    = matcl::pod_workspace<VR>;
        rworkspace RWORK    = rworkspace(lrwork);
        VR* ptr_RWORK       = RWORK.ptr();

        lapack::heevd(job_U, uplo, N, lap(TA), TA_ld, lap(eig), lap(ptr_WORK), lwork, 
                      ptr_RWORK, lrwork, ptr_IWORK, liwork, &info);
    }
};

template<class V, bool Is_real = md::is_float_real_scalar<V>::value>
struct eval_sbev_hbev
{
    static void eval(const char* job_U, const char* uplo, Integer N, Integer LD, V* TA, Integer TA_ld,
                     V* eig, V* Z, Integer Z_ld, Integer& info)
    {
        Integer lwork       = std::max(1,3*N);

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

        lapack::sbev(job_U, uplo, N, LD, lap(TA), TA_ld,lap(eig), lap(Z), Z_ld, lap(ptr_WORK), &info);
    };
};

template<class V>
struct eval_sbev_hbev<V,false>
{
    using VR = typename md::real_type<V>::type;

    static void eval(const char* job_U, const char* uplo, Integer N, Integer LD, V* TA, Integer TA_ld,
                     VR* eig, V* Z, Integer Z_ld, Integer& info)
    {
        Integer lwork       = N;
        Integer lrwork      = std::max(1,3*N);

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

        using rworkspace    = matcl::pod_workspace<VR>;
        rworkspace RWORK    = rworkspace(lrwork);
        VR* ptr_RWORK       = RWORK.ptr();

        lapack::hbev(job_U, uplo, N, LD, lap(TA), TA_ld,lap(eig), lap(Z), Z_ld, 
                     lap(ptr_WORK), lap(ptr_RWORK), &info);
    };
};

template<class V, bool Is_real = md::is_float_real_scalar<V>::value>
struct eval_sbevd_hbevd
{
    static void eval(const char* job_U, const char* uplo, Integer N, Integer LD, V* TA, Integer TA_ld,
                     V* eig, V* Z, Integer Z_ld, Integer& info)
    {
        V       work_query;
        Integer iwork_query;

        lapack::sbevd(job_U, uplo, N, LD, lap(TA), TA_ld,lap(eig), lap(Z), Z_ld, lap(&work_query), 
                      -1, &iwork_query, -1, &info);

        Integer lwork       = (Integer)real(work_query);
        Integer liwork      = iwork_query;

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

        using iworkspace    = matcl::pod_workspace<Integer>;
        iworkspace IWORK    = iworkspace(liwork);
        Integer* ptr_IWORK  = IWORK.ptr();

        lapack::sbevd(job_U, uplo, N, LD, lap(TA), TA_ld,lap(eig), lap(Z), Z_ld, lap(ptr_WORK), lwork, 
                      ptr_IWORK, liwork, &info);

    };
};
template<class V>
struct eval_sbevd_hbevd<V,false>
{
    using VR = typename md::real_type<V>::type;

    static void eval(const char* job_U, const char* uplo, Integer N, Integer LD, V* TA, Integer TA_ld,
                     VR* eig, V* Z, Integer Z_ld, Integer& info)
    {
        V       work_query;
        VR      rwork_query;
        Integer iwork_query;

        lapack::hbevd(job_U, uplo, N, LD, lap(TA), TA_ld,lap(eig), lap(Z), Z_ld, lap(&work_query), 
                      -1, lap(&rwork_query), -1, &iwork_query, -1, &info);

        Integer lwork       = (Integer)real(work_query);
        Integer liwork      = iwork_query;
        Integer lrwork      = (Integer)rwork_query;

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

        using rworkspace    = matcl::pod_workspace<VR>;
        rworkspace RWORK    = rworkspace(lrwork);
        VR* ptr_RWORK       = RWORK.ptr();

        using iworkspace    = matcl::pod_workspace<Integer>;
        iworkspace IWORK    = iworkspace(liwork);
        Integer* ptr_IWORK  = IWORK.ptr();

        lapack::hbevd(job_U, uplo, N, LD, lap(TA), TA_ld,lap(eig), lap(Z), Z_ld, lap(ptr_WORK), lwork, 
                      lap(ptr_RWORK), lrwork, ptr_IWORK, liwork, &info);

    };
};

struct reorder_diag
{
    static void eval(const matcl::Matrix& ind, schur_decomposition& sd)
    {
        Matrix J1   = matcl::find(ind == 1);
        Matrix J2   = matcl::find(ind != 1);

        Matrix J    = (mat_col(), vec(J1), vec(J2));

        Matrix TA2  = sd.m_TA(J, J);
        Matrix eig  = TA2;

        TA2.add_struct(predefined_struct_type::diag);

        sd.m_TA     = TA2;
        sd.m_eig    = eig;

        if (sd.m_has_U)
        {
            Matrix U        = sd.m_U_factor;
            bool has_unit   = is_unitary(U.get_struct());

            U               = U(colon(), J);

            if (has_unit)
            {
                struct_flag sf_u;
                sf_u.set_user(unitary_flag());
                U.add_struct(sf_u);
            };

            sd.m_U_factor   = U;
        };
    };
};

template<class V, class S>
struct schur_str
{
    using Mat   = raw::Matrix<V,S>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval(const Mat& A, schur_decomposition& sd, schur_sym_alg alg, bool with_U)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval(Ac,sd,alg,with_U);
    };
    static void eval_reorder(const Mat& TA, const matcl::Matrix& ind, schur_decomposition& sd)
    {
        bool is_diag = matcl::is_diag(Matrix(TA,false));

        if (is_diag == false)
        {
            Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(TA);
            return schur_str<V,struct_dense>::eval_reorder(Ac,ind,sd);
        };

        return reorder_diag::eval(ind, sd);
    };
    static void eval_eigenvec(const Mat& TA, matcl::Matrix& ind, const schur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(TA);
        return schur_str<V,struct_dense>::eval_eigenvec(Ac,ind,sd,XL,XR,comp_left,comp_right);
    };

    static void eval_cond_eig(Matrix& ret, const Mat& A, const schur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& VL, const Matrix& VR)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval_cond_eig(ret, Ac, sd, ind, VL, VR);
    };
    static void eval_cond_vec(Matrix& ret, const Mat& A, const schur_decomposition& sd, const matcl::Matrix& ind)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval_cond_vec(ret, Ac, sd, ind);
    };

    static Real eval_cond_eig_cluster(const Mat& A, Integer M)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval_cond_eig_cluster(Ac, M);
    };
    static Real eval_cond_subspace(const Mat& A, Integer M)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval_cond_eig_cluster(Ac, M);
    };
};

template<class V>
struct schur_str<V,struct_banded>
{
    using Mat   = raw::Matrix<V,struct_banded>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval(const Mat& A, schur_decomposition& sd, schur_sym_alg alg, bool with_U)
    {
        if (schur_is_sym_her<V,struct_banded>::eval(A))
            return compute_sym_her(A, sd, alg, with_U);

        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval(Ac,sd,alg,with_U);
    };
    
    static void eval_reorder(const Mat& TA, const matcl::Matrix& ind, schur_decomposition& sd)
    {
        bool is_diag = matcl::is_diag(Matrix(TA,false));

        if (is_diag == false)
        {
            Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(TA);
            return schur_str<V,struct_dense>::eval_reorder(Ac,ind,sd);
        };

        return reorder_diag::eval(ind, sd);
    };

    static void eval_eigenvec(const Mat& TA, matcl::Matrix& ind, const schur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(TA);
        return schur_str<V,struct_dense>::eval_eigenvec(Ac,ind,sd,XL,XR,comp_left,comp_right);
    };
    static void eval_cond_eig(Matrix& ret, const Mat& A, const schur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& VL, const Matrix& VR)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval_cond_eig(ret, Ac, sd, ind, VL, VR);
    };
    static void eval_cond_vec(Matrix& ret, const Mat& A, const schur_decomposition& sd, const matcl::Matrix& ind)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval_cond_vec(ret, Ac, sd, ind);
    };
    static Real eval_cond_eig_cluster(const Mat& A, Integer M)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval_cond_eig_cluster(Ac, M);
    };
    static Real eval_cond_subspace(const Mat& A, Integer M)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_str<V,struct_dense>::eval_cond_eig_cluster(Ac, M);
    };

    static void compute_sym_her(const Mat& A, schur_decomposition& sd, schur_sym_alg alg, bool with_U)
    {
        using RC        =  typename details::real_type<V>::type;
        using mat_real  = raw::Matrix<RC,struct_dense>;
        
        using matcl::details::lap;

        Integer N           = A.rows();

        Integer U_N         = with_U? N : 1;
        const char* job_U   = with_U ? "V" : "N";

        if (A.has_diag(0) == false)
            throw error::band_matrix_with_main_diag_required(A.first_diag(), A.last_diag());

        Mat TA              = A.make_unique();
        Integer LD          = TA.number_subdiagonals();

        TA.set_struct(struct_flag());        

        mat_real    eig(A.get_type(), N, 1);
        Mat_D       Q(A.get_type(),U_N,U_N);    

        matcl::lapack::i_type info = 0;

        V* TA_ptr           = TA.rep_ptr() + TA.first_elem_diag(0);

        // rrr algorithm is not available for band matrices
        if (alg == schur_sym_alg::rrr)
            alg = schur_sym_alg::dc;

        if (alg == schur_sym_alg::qr)
        {            
            eval_sbev_hbev<V>::eval(job_U, "L", N, LD, TA_ptr, TA.ld(), eig.ptr(), Q.ptr(), Q.ld(), info);

            if (info)
            {
                sd.clear();
                throw error::error_schur();
            }

            if (with_U)
            {
                sd.m_U_factor   = Matrix(Q,true);
                sd.m_has_U      = true;

                struct_flag sf_u;
                sf_u.set_user(unitary_flag());

                sd.m_U_factor.add_struct(sf_u);
            }
            else
            {
                sd.m_has_U  = false;
            };
        }        
        else if (alg == schur_sym_alg::dc)
        {
            eval_sbevd_hbevd<V>::eval(job_U, "L", N, LD, TA_ptr, TA.ld(), eig.ptr(), Q.ptr(), Q.ld(), info);

            if (info)
            {
                sd.clear();
                throw error::error_schur();
            }

            if (with_U)
            {
                sd.m_U_factor   = Matrix(Q,true);
                sd.m_has_U      = true;

                struct_flag sf_u;
                sf_u.set_user(unitary_flag());

                sd.m_U_factor.add_struct(sf_u);
            }
            else
            {
                sd.m_has_U  = false;
            };
        }
        else
        {
            matcl_assert(0,"unknown schur_sym_alg");
            throw error::error_general("unknown schur_sym_alg");
        };

        sd.m_eig    = Matrix(eig, true);
        sd.m_TA     = bdiag(sd.m_eig);
    
        sd.m_TA.add_struct(predefined_struct_type::diag);

        sd.test_factors();
    };
};

template<class V>
struct schur_str<V,struct_dense>
{
    using Mat   = raw::Matrix<V,struct_dense>;

    static void eval(const Mat& A, schur_decomposition& sd, schur_sym_alg alg, bool with_U)
    {    
        if (schur_is_sym_her<V,struct_dense>::eval(A))
            return compute_sym_her(A, sd, alg, with_U);
        else
            return compute_nsym(A,sd, with_U);
    };
    
    static void compute_nsym(const Mat& A, schur_decomposition& sd, bool with_U)
    {
        using matcl::details::lap;

        Integer N       = A.rows();

        using VC        = typename complex_type<V>::type;
        using mat_compl = raw::Matrix<VC,struct_dense> ;

        Mat TA          = A.make_unique();
        TA.set_struct(struct_flag());
    
        Integer U_N         = with_U? N : 1;
        const char* job_U   = with_U ? "V" : "N";

        Integer     Sdim;
        Mat         Q(A.get_type(),U_N,U_N);    
        mat_compl   eig(A.get_type(), N, 1);
        V           work_query;        

        matcl::lapack::i_type info = 0;  

        lapack::gees(job_U, "N", 0, N, lap(TA.ptr()), TA.ld(), lap(&Sdim),
                          lap(eig.ptr()), lap(Q.ptr()), Q.ld(), lap(&work_query), -1, 0, &info);

        Integer lwork       = (Integer) real(work_query);

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

        lapack::gees(job_U, "N", 0, N, lap(TA.ptr()), TA.ld(), lap(&Sdim),
                          lap(eig.ptr()), lap(Q.ptr()), Q.ld(), lap(ptr_WORK), lwork, 0, &info);

        if (info)
        {
            sd.clear();
            throw error::error_schur();
        };
        
        sd.m_TA     = Matrix(TA,true);
        sd.m_eig    = Matrix(eig, true);

        if (with_U)
        {
            sd.m_U_factor   = Matrix(Q,true);
            sd.m_has_U      = true;

            struct_flag sf_u;
            sf_u.set_user(unitary_flag());

            sd.m_U_factor.add_struct(sf_u);
        }
        else
        {
            sd.m_has_U  = false;
        };
    
        if (md::is_complex<V>::value == false)
        {
            Integer ld  = get_ld(sd.m_TA, 1);

            if (ld == 1)
                sd.m_TA.add_struct(md::predefined_struct::get_qtriu(sd.m_TA.get_struct()));
            else
                sd.m_TA.add_struct(md::predefined_struct::get_triu(sd.m_TA.get_struct(),0,
                                                                   is_real_matrix(sd.m_TA)));
        }
        else
        {
            sd.m_TA.add_struct(md::predefined_struct::get_triu(sd.m_TA.get_struct(),0,
                                                               is_real_matrix(sd.m_TA)));
        };

        sd.test_factors();
    };

    static void eval_tridiagonal(const Mat& A, schur_decomposition& sd, schur_sym_alg alg, bool with_U)
    {
        Matrix mat_A(A, false);
        mat_A   = matcl::select_band(mat_A, -1, 1);        

        using Mat_B = raw::Matrix<V,struct_banded>;
        mat_A       = convert(mat_A, Mat_B::matrix_code);
        Mat_B& AB   = mat_A.get_impl_unique<Mat_B>(); 

        mat_A = Matrix();
        return schur_str<V,struct_banded>::eval(AB, sd, alg, with_U);
    };

    static void compute_sym_her(const Mat& A, schur_decomposition& sd, schur_sym_alg alg, bool with_U)
    {
        if (get_ld(Matrix(A,false), 1) == 1)
            return eval_tridiagonal(A, sd, alg, with_U);

        using RC        =  typename details::real_type<V>::type;
        using mat_real  = raw::Matrix<RC,struct_dense>;
        
        using matcl::details::lap;

        Integer N           = A.rows();

        Integer U_N         = with_U? N : 1;
        const char* job_U   = with_U ? "V" : "N";

        Mat TA              = A.make_unique();
        TA.set_struct(struct_flag());

        mat_real eig(A.get_type(), N, 1);

        matcl::lapack::i_type info = 0;

        if (alg == schur_sym_alg::qr)
        {
            V work_query;
            lapack::heev(job_U, "L", N, lap(TA.ptr()), TA.ld(),lap(eig.ptr()),
                          lap(&work_query), -1, &info);
    
            Integer lwork = (Integer) real(work_query);

            using VTR_pod       = matcl::pod_type<V>;
            using workspace     = matcl::pod_workspace<VTR_pod>;
            workspace WORK      = workspace(lwork);
            V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

            lapack::heev(job_U, "L", N, lap(TA.ptr()), TA.ld(), lap(eig.ptr()),
                          lap(ptr_WORK), lwork, &info);

            if (info)
            {
                sd.clear();
                throw error::error_schur();
            }

            if (with_U)
            {
                sd.m_U_factor   = Matrix(TA,true);
                sd.m_has_U      = true;

                struct_flag sf_u;
                sf_u.set_user(unitary_flag());

                sd.m_U_factor.add_struct(sf_u);
            }
            else
            {
                sd.m_has_U  = false;
            };
        }
        else if (alg == schur_sym_alg::dc)
        {
            eval_syevd_heevd<V>::eval(job_U, "L", N, TA.ptr(), TA.ld(), eig.ptr(), info);

            if (info)
            {
                sd.clear();
                throw error::error_schur();
            }

            if (with_U)
            {
                sd.m_U_factor   = Matrix(TA,true);
                sd.m_has_U      = true;

                struct_flag sf_u;
                sf_u.set_user(unitary_flag());

                sd.m_U_factor.add_struct(sf_u);
            }
            else
            {
                sd.m_has_U  = false;
            };
        }
        else if (alg == schur_sym_alg::rrr)
        {
            using VC        = typename complex_type<V>::type;
            using mat_compl = raw::Matrix<VC,struct_dense>;

            using VTR_pod       = matcl::pod_type<V>;
            using workspace     = matcl::pod_workspace<VTR_pod>;
            using r_workspace   = matcl::pod_workspace<RC>;
            using i_workspace   = matcl::pod_workspace<Integer>;

            Mat         Q(A.get_type(),U_N,U_N);
            V           work_query;
            RC          rwork_query;
            Integer     iwork_query;
            Integer     M;

            i_workspace isuppz = i_workspace(2 * max(1, N));

            lapack::heevr(job_U, "A", "L", N, lap(TA.ptr()), TA.ld(), 0, 0, 0, 0, 0, &M, lap(eig.ptr()),
                          lap(Q.ptr()), Q.ld(), lap(isuppz.ptr()),
                          lap(&work_query), -1, lap(&rwork_query), -1, lap(&iwork_query), -1, &info);
    
            Integer lwork   = (Integer) real(work_query);
            Integer lrwork  = (Integer) rwork_query;
            Integer liwork  = iwork_query;

            workspace WORK      = workspace(lwork);
            r_workspace rwork   = r_workspace(lrwork);
            i_workspace iwork   = i_workspace(liwork);

            V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

            lapack::heevr(job_U, "A", "L", N, lap(TA.ptr()), TA.ld(), 0, 0, 0, 0, 0, &M, lap(eig.ptr()),
                          lap(Q.ptr()), Q.ld(), lap(isuppz.ptr()),
                          lap(ptr_WORK), lwork, lap(rwork.ptr()), lrwork, lap(iwork.ptr()), liwork, &info);

            if (info)
            {
                sd.clear();
                throw error::error_schur();
            }

            if (with_U)
            {
                // NOTE: don't set unitary flag for Q, RRR algorithm may yield a numerically non-unitary matrix
                sd.m_U_factor   = Matrix(Q,true);
                sd.m_has_U      = true;
            }
            else
            {
                sd.m_has_U  = false;
            };            
        }
        else
        {
            matcl_assert(0,"unknown schur_sym_alg");
            throw error::error_general("unknown schur_sym_alg");
        };

        sd.m_eig    = Matrix(eig, true);
        sd.m_TA     = bdiag(sd.m_eig);
    
        sd.m_TA.add_struct(predefined_struct_type::diag);

        sd.test_factors();
    };

    static void eval_reorder(const Mat& TA0, const matcl::Matrix& ind, schur_decomposition& sd)
    {    
        Matrix mat_TA0  = Matrix(TA0,false);
        Integer ld      = get_ld(mat_TA0, 1);
        Integer ud      = get_ud(mat_TA0, 1);

        bool TA_was_diagonal = (ld == 0 && ud == 0);

        bool with_U     = sd.m_has_U;        
        bool has_unit   = is_unitary(sd.m_U_factor.get_struct());;

        Mat TA          = TA0.make_unique();
        TA.get_struct().reset();

        if (with_U == true)
            sd.m_U_factor   = convert(sd.m_U_factor, Mat::matrix_code);

        Mat Q           = with_U? Mat(sd.m_U_factor.get_impl_unique<Mat>(), Mat::copy_is_safe())
                                : Mat(TA.get_type());
        Q.get_struct().reset();

        V* r_TA         = TA.ptr();
        Integer TA_ld   = TA.ld();
        Integer N       = TA.rows();
        V* r_Q          = Q.ptr();
        Integer Q_ld    = std::max(1,Q.ld());
        Integer Q_r     = Q.rows();
        const Integer* I= ind.get_array<Integer>();

        using VC        = typename md::complex_type<V>::type;
        using VR        = typename md::real_type<V>::type;
        using DM_compl  = raw::Matrix<VC,struct_dense>;

        DM_compl w(ti::get_ti(TA),N,1);

        //transform quasi-triangular matrix to schur form
        lapack::qtriu2schur(TA.rows(), lap(r_TA), TA_ld, with_U, lap(r_Q), Q_r, Q_ld);
        
        matcl::lapack::i_type info = 0, M;    

        char compq  = with_U? 'V' : 'N';    // update Schur vectors

        using matcl::details::lap;

        V work_query;

        lapack::trsen3(&compq, lap(I), N, Q_r, lap(r_TA), TA_ld, lap(r_Q),
                           Q_ld, lap(w.ptr()), &M, lap(&work_query), -1, &info);

        Integer lwork = (Integer) real(work_query);

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(lwork);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());
        
        lapack::trsen3(&compq, lap(I), N, Q_r, lap(r_TA), TA_ld, lap(r_Q),
                           Q_ld, lap(w.ptr()), &M, lap(ptr_WORK), lwork, &info);        

        if (info)
        {
            sd.clear();
            throw error::error_schur_sel_comp();
        };

        sd.m_eig = Matrix(w,false);

        Matrix mat_TA = Matrix(TA,false);

        if (TA_was_diagonal)
        {
            mat_TA.add_struct(predefined_struct_type::diag);
        }
        else if (md::is_complex<V>::value == false)
        {            
            if (ld == 1)
                mat_TA.add_struct(predefined_struct_type::qtriu);
            else
                mat_TA.add_struct(predefined_struct_type::triu);
        }
        else
        {
            mat_TA.add_struct(predefined_struct_type::triu);
        }

        sd.m_TA = mat_TA;

        if (with_U)
        {
            struct_flag sf_u;

            if (has_unit)
                sf_u.set_user(unitary_flag());

            sd.m_U_factor = Matrix(Q,false);
            sd.m_U_factor.set_struct(sf_u);
        };
    };

    static void eval_eigenvec(const Mat& TA0, matcl::Matrix& ind, const schur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        static_assert(std::is_same<V,typename std::decay<V>::type>::value,"");
        // matrix TA is modified by trevc but restored on exit;
        // this should be thread safe after making TA unique (instances of Matrix cannot be shared
        // between threads)
        Mat TA(TA0.get_type());

        if (TA0.is_unique() == true)
        {
            TA.assign_to_fresh(TA0);
        }
        else
        {
            TA.assign_to_fresh(TA0.copy());
            // we reset TA in schur_decomposition to avoid creating new copy if this function
            // is called again
            sd.m_TA             = Matrix(TA,false);
        };

        const char* SIDE        = (comp_left && comp_right)? "B" : (comp_left? "L" : "R");
        Integer* SELECT         = ind.get_array_unique<Integer>();
        Integer N               = TA.rows();
        V* ptr_T                = TA.ptr();
        Integer LDT             = TA.ld();

        Integer M               = make_complex_eigenvectors<V>::calc_sel_size(TA, SELECT, N);
        Integer ML              = comp_left? M : 1;
        Integer MR              = comp_right? M : 1;

        Mat VL                  = Mat(TA.get_type(), N, ML);
        Mat VR                  = Mat(TA.get_type(), N, MR);

        Integer info;

        using VTR_pod           = matcl::pod_type<V>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        Integer lwork           = 3*N;
        workspace WORK          = workspace(lwork);
        V* ptr_WORK             = reinterpret_cast<V*>(WORK.ptr());

        lapack::trevc(SIDE, "S", SELECT, N, lap(ptr_T), LDT, lap(VL.ptr()), VL.ld(), lap(VR.ptr()), VR.ld(), 
                      M, M, lap(ptr_WORK), &info);

        if (info != 0)
            throw error::error_general("illegal argument passed to trevc");

        Matrix mat_VL           = Matrix(VL,false);
        Matrix mat_VR           = Matrix(VR,false);

        if (sd.m_has_U)
        {
            // form U*VL, U*VR
            if (comp_left)
                mat_VL          = mmul(sd.m_U_factor, mat_VL, trans_type::no_trans);

            if (comp_right)
                mat_VR          = mmul(sd.m_U_factor, mat_VR, trans_type::no_trans);
        };

        make_complex_eigenvectors<V>::eval(TA, ind, mat_VL, mat_VR, XL, XR, comp_left, comp_right);
        return;
    };
    static void eval_cond_eig(Matrix& ret, const Mat& TA, const schur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& WL, const Matrix& WR)
    {
        (void)sd;

        using VR                = typename md::real_type<V>::type;
        using Mat_R             = raw::Matrix<VR,struct_dense>;

        const char* JOB         = "E";
        const char* HOWMNY      = "S";
        const Integer* SELECT   = ind.get_array<Integer>();
        Integer N               = TA.rows();
        Integer M               = make_complex_eigenvectors<V>::calc_sel_size(TA, SELECT, N);
        const V* ptr_TA         = TA.ptr();
        Integer LDT             = TA.ld();

        Mat_R S(TA.get_type(), M, 1);
        VR* ptr_S               = S.ptr();

        Mat WL_R(TA.get_type());
        Mat WR_R(TA.get_type());

        make_complex_eigenvectors<V>::complex_vectors_to_real(TA, WL, WR, WL_R, WR_R, ind);

        V* ptr_VL               = WL_R.ptr();
        V* ptr_VR               = WR_R.ptr();
        Integer ld_VL           = WL_R.ld();
        Integer ld_VR           = WR_R.ld();

        V* ptr_null             = nullptr;
        VR* ptr_null_r          = nullptr;
        Integer* ptr_null_i     = nullptr;

        Integer info;

        lapack::trsna( JOB, HOWMNY, SELECT, N, lap(ptr_TA), LDT, lap(ptr_VL), ld_VL, lap(ptr_VR), ld_VR, 
                      lap(ptr_S), lap(ptr_null_r), M, M, lap(ptr_null), N, ptr_null_i, ptr_null_r, &info);

        if (info != 0)
            throw error::error_general("illegal argument passed to trsna");

        ret = matcl::Matrix(S,false);
        return;
    };
    static void eval_cond_vec(Matrix& ret, const Mat& TA, const schur_decomposition& sd, const matcl::Matrix& ind)
    {
        (void)sd;

        using VR                = typename md::real_type<V>::type;
        using Mat_R             = raw::Matrix<VR,struct_dense>;

        bool is_real_type       = md::is_float_real_scalar<V>::value;

        const char* JOB         = "V";
        const char* HOWMNY      = "S";
        const Integer* SELECT   = ind.get_array<Integer>();
        Integer N               = TA.rows();
        Integer M               = make_complex_eigenvectors<V>::calc_sel_size(TA, SELECT, N);
        const V* ptr_TA         = TA.ptr();
        Integer LDT             = TA.ld();

        Mat_R SEP(TA.get_type(), M, 1);
        VR* ptr_SEP             = SEP.ptr();

        using VTR_pod           = matcl::pod_type<V>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        Integer lwork           = N*(N+6);
        workspace WORK          = workspace(lwork);
        V* ptr_WORK             = reinterpret_cast<V*>(WORK.ptr());

        using iworkspace        = matcl::pod_workspace<Integer>;
        Integer liwork          = is_real_type? 2*(N-1) : 0;
        iworkspace IWORK        = iworkspace(liwork);
        Integer* ptr_IWORK      = IWORK.ptr();

        using rworkspace        = matcl::pod_workspace<VR>;
        Integer lrwork          = is_real_type? 0 : N;
        rworkspace RWORK        = rworkspace(lrwork);
        VR* ptr_RWORK           = RWORK.ptr();

        V* ptr_null             = nullptr;
        VR* ptr_null_r          = nullptr;

        Integer info;

        lapack::trsna( JOB, HOWMNY, SELECT, N, lap(ptr_TA), LDT, lap(ptr_null), 1, lap(ptr_null), 1, 
                      lap(ptr_null_r), lap(ptr_SEP), M, M, lap(ptr_WORK), N, ptr_IWORK, ptr_RWORK, &info);

        if (info != 0)
            throw error::error_general("illegal argument passed to trsna");

        ret = matcl::Matrix(SEP,false);
        return;
    };

    static Real eval_cond_eig_cluster(const Mat& TA, Integer M)
    {
        using VR        = typename md::real_type<V>::type;

        const char* JOB = "E";
        Integer N       = TA.rows();
        const V* ptr_TA = TA.ptr();
        Integer ld_TA   = TA.ld();

        //check if M-th eigenvalue is complex
        if (M > 0 && M < N)
        {
            if (ptr_TA[M + (M-1)*ld_TA] != V(0.0))
                M       = M + 1;
        };

        VR S, SEP;

        V       work_query;
        Integer iwork_query;
        Integer info;

        lapack::trsen_cond(JOB, N, M, lap(ptr_TA), ld_TA, S, SEP, lap(&work_query), -1, 
                           &iwork_query, -1, info);

        Integer lwork   = (Integer)real(work_query);
        Integer liwork  = iwork_query;

        using VTR_pod           = matcl::pod_type<V>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        workspace WORK          = workspace(lwork);
        V* ptr_WORK             = reinterpret_cast<V*>(WORK.ptr());

        using iworkspace        = matcl::pod_workspace<Integer>;
        iworkspace IWORK        = iworkspace(liwork);
        Integer* ptr_IWORK      = IWORK.ptr();

        lapack::trsen_cond(JOB, N, M, lap(ptr_TA), ld_TA, S, SEP, lap(ptr_WORK), lwork, 
                           ptr_IWORK, liwork, info);

        if (info != 0)
            throw error::error_general("illegal argument passed to trsen_cond");

        return S;
    };

    static Real eval_cond_subspace(const Mat& TA, Integer M)
    {
        using VR        = typename md::real_type<V>::type;

        const char* JOB = "V";
        Integer N       = TA.rows();
        const V* ptr_TA = TA.ptr();
        Integer ld_TA   = TA.ld();

        //check if M-th eigenvalue is complex
        if (M > 0 && M < N)
        {
            if (ptr_TA[M + (M-1)*ld_TA] != V(0.0))
                M       = M + 1;
        };

        VR S, SEP;

        V       work_query;
        Integer iwork_query;
        Integer info;

        lapack::trsen_cond(JOB, N, M, lap(ptr_TA), ld_TA, S, SEP, lap(&work_query), -1, 
                           &iwork_query, -1, info);

        Integer lwork   = (Integer)real(work_query);
        Integer liwork  = iwork_query;

        using VTR_pod           = matcl::pod_type<V>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        workspace WORK          = workspace(lwork);
        V* ptr_WORK             = reinterpret_cast<V*>(WORK.ptr());

        using iworkspace        = matcl::pod_workspace<Integer>;
        iworkspace IWORK        = iworkspace(liwork);
        Integer* ptr_IWORK      = IWORK.ptr();

        lapack::trsen_cond(JOB, N, M, lap(ptr_TA), ld_TA, S, SEP, lap(ptr_WORK), lwork, 
                           ptr_IWORK, liwork, info);

        if (info != 0)
            throw error::error_general("illegal argument passed to trsen_cond");

        return SEP;
    };
};

template<class V>
struct estim_rcond_vec_diag
{
    using VR    = typename md::real_type_int_real<V>::type;
    using Mat_R = raw::Matrix<VR,struct_dense>;
    using Mat   = raw::Matrix<V,struct_dense>;

    static void eval(Matrix& ret, const Mat& D0, const Matrix& ind)
    {
        //  Then the reciprocal condition number is

        //          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
        //  where sigma-min denotes the smallest singular value. We approximate
        //  the smallest singular value by the reciprocal of an estimate of the
        //  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is
        //  defined to be abs(T(1,1)).

        // for diagonal matrix this is min_{k != j}(|D_k - D_j|), where D_i is i-th eigenvalue        

        const Mat& D = D0.make_explicit();

        Integer N = D.size();

        if (N == 0)
        {
            Mat_R res = Mat_R(D.get_type(), 0, 1);
            ret = Matrix(res,false);
            return;
        };        

        const Integer* ptr_I    = ind.get_array<Integer>();
        const V* ptr_D          = D.ptr();

        if (N == 1)
        {
            ret = VR(abs(ptr_D[0]));
            return;
        };        

        Integer M   = 0;
        for (Integer i = 0; i < N; ++i)
        {
            if (ptr_I[i] != 0)
                ++M;
        };

        Mat_R res   = Mat_R(D.get_type(), M, 1);

        VR* ptr_res = res.ptr();

        // for diagonal matrix this is min_{k != j}(|D_k - D_j|), where D_i is i-th eigenvalue        
        for (Integer i = 0, pos = 0; i < N; ++i)
        {
            bool sel    = ptr_I[i] != 0;
            
            if (sel == false)
                continue;

            V D_i       = ptr_D[i];
            VR ma       = constants::inf<VR>();

            for (Integer k = 0; k < N; ++k)
            {
                if (k == i)
                    continue;

                V D_k   = ptr_D[k];
                VR dif  = abs(D_k - D_i);

                if (dif < ma)
                    ma  = dif;
            };

            ptr_res[pos] = ma;
            ++pos;
        };

        ret = Matrix(res,false);
    };

    static Real eval_subspace(const Mat& D0, Integer M)
    {
        //  The reciprocal condition number of the right invariant subspace
        //  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
        //  SEP is defined as the separation of T11 and T22:
        //
        //                     sep( T11, T22 ) = sigma-min( C )
        //
        //  where sigma-min(C) is the smallest singular value of the
        //  n1*n2-by-n1*n2 matrix
        //
        //     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
        //
        //  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
        //  product. We estimate sigma-min(C) by the reciprocal of an estimate of
        //  the 1-norm of inverse(C).

        // for diagonal matrix this is min_{k in I1, j in I2}(|D_k - D_j|), where D_i is i-th eigenvalue        

        const Mat& D = D0.make_explicit();

        Integer N = D.size();

        if (N == 0)
            return 0.0;

        const V* ptr_D          = D.ptr();

        if (N == 1)
            return abs(ptr_D[0]);

        // for diagonal matrix this is min_{k in I1, j in I2}(|D_k - D_j|), where D_i is i-th eigenvalue 
        VR ma       = constants::inf<VR>();

        if (M == 0 || M == N)
        {
            for (Integer k = 0; k < M; ++k)
            {
                V D_k   = ptr_D[k];            

                VR dif  = abs(D_k);

                if (dif < ma)
                    ma  = dif;
            };

            return ma;
        };

        for (Integer k = 0; k < M; ++k)
        {
            V D_k       = ptr_D[k];            

            for (Integer j = M; j < N; ++j)
            {
                V D_j   = ptr_D[j];
                VR dif  = abs(D_k - D_j);

                if (dif < ma)
                    ma  = dif;
            };
        };

        return ma;
    };
};

template<class V, class S>
struct schur_impl
{
    using Mat = raw::Matrix<V,S>;

    static void eval(const Mat& A, schur_decomposition& sd, schur_sym_alg alg, bool with_U)
    {
        return schur_str<V,S>::eval(A,sd, alg, with_U);
    };

    static void eval_reorder(const Mat& A, const matcl::Matrix& ind, schur_decomposition& sd)
    {
        return schur_str<V,S>::eval_reorder(A,ind,sd);
    };
    static void eval_eigenvec(const Mat& A, matcl::Matrix& ind, const schur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        return schur_str<V,S>::eval_eigenvec(A,ind,sd,XL,XR,comp_left,comp_right);
    };
    static void eval_cond_eig(Matrix& ret, const Mat& A, const schur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& VL, const Matrix& VR)
    {
        return schur_str<V,S>::eval_cond_eig(ret,A,sd,ind,VL,VR);
    };
    static void eval_cond_vec(Matrix& ret, const Mat& A, const schur_decomposition& sd, const matcl::Matrix& ind)
    {
        return schur_str<V,S>::eval_cond_vec(ret,A,sd,ind);
    };

    static Real eval_cond_eig_cluster(const Mat& TA, Integer M)
    {
        return schur_str<V,S>::eval_cond_eig_cluster(TA, M);
    };
    static Real eval_cond_subspace(const Mat& TA, Integer M)
    {
        return schur_str<V,S>::eval_cond_subspace(TA, M);
    };
};

template<class S>
struct schur_impl<Integer, S>
{
    using Mat   = raw::Matrix<Integer,S>;
    using Mat_D = raw::Matrix<Real,struct_dense>;

    static void eval(const Mat& A, schur_decomposition& sd, schur_sym_alg alg, bool with_U)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_impl<Real,struct_dense>::eval(Ac,sd, alg, with_U);
    };

    static void eval_reorder(const Mat& A, const matcl::Matrix& ind, schur_decomposition& sd)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_impl<Real,struct_dense>::eval_reorder(Ac,ind,sd);
    };

    static void eval_eigenvec(const Mat& A, matcl::Matrix& ind, const schur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_impl<Real,struct_dense>::eval_eigenvec(Ac,ind,sd,XL,XR,comp_left,comp_right);
    };

    static void eval_cond_eig(Matrix& ret, const Mat& A, const schur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& VL, const Matrix& VR)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_impl<Real,struct_dense>::eval_cond_eig(ret, Ac,sd, ind, VL, VR);
    };
    static void eval_cond_vec(Matrix& ret, const Mat& A, const schur_decomposition& sd, const matcl::Matrix& ind)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return schur_impl<Real,struct_dense>::eval_cond_vec(ret, Ac,sd, ind);
    };

    static Real eval_cond_eig_cluster(const Mat& TA, Integer M)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(TA);
        return schur_impl<Real,struct_dense>::eval_cond_eig_cluster(Ac, M);
    };
    static Real eval_cond_subspace(const Mat& TA, Integer M)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(TA);
        return schur_impl<Real,struct_dense>::eval_cond_subspace(Ac, M);
    };
};

struct schur_visitor : public extract_type_switch<void, schur_visitor,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, schur_decomposition& sd, schur_sym_alg alg,
                     bool with_U)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::schur_impl<V,S>::eval(mat, sd, alg, with_U);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, schur_decomposition& sd, schur_sym_alg alg,
                            bool with_U)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, sd, alg, with_U);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, schur_decomposition&, schur_sym_alg,bool)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, schur_decomposition&, schur_sym_alg,bool)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};

struct schur_visitor_reorder : public extract_type_switch<void, schur_visitor_reorder,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, const matcl::Matrix& ind, schur_decomposition& sd)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::schur_impl<V,S>::eval_reorder(mat, ind, sd);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, const matcl::Matrix& ind, schur_decomposition& sd)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, ind, sd);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const matcl::Matrix&, schur_decomposition&)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, const matcl::Matrix&, schur_decomposition&)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};

struct schur_visitor_eigvec : public extract_type_switch<void, schur_visitor_eigvec,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, matcl::Matrix& ind, const schur_decomposition& sd,
                    Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::schur_impl<V,S>::eval_eigenvec(mat, ind, sd, XL, XR, comp_left, comp_right);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, matcl::Matrix& ind, const schur_decomposition& sd,
                    Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, ind, sd,XL,XR,comp_left,comp_right);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, matcl::Matrix&, const schur_decomposition&,
                     Matrix&, Matrix&, bool, bool)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, matcl::Matrix&, const schur_decomposition&,
                            Matrix&, Matrix&, bool, bool)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};

struct schur_visitor_cond_eig : public extract_type_switch<void, schur_visitor_cond_eig,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, const schur_decomposition& sd,
                     const Matrix& VL, const Matrix& VR, const matcl::Matrix& ind)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::schur_impl<V,S>::eval_cond_eig(ret, mat, sd, ind, VL, VR);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, Matrix& ret, const schur_decomposition& sd,
                     const Matrix& VL, const Matrix& VR, const matcl::Matrix& ind)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, ret, sd, VL, VR, ind);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const schur_decomposition&,
                     const Matrix&, const Matrix&, const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, Matrix&, const schur_decomposition&,
                     const Matrix&, const Matrix&, const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};
struct schur_visitor_cond_vec : public extract_type_switch<void, schur_visitor_cond_vec,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, const schur_decomposition& sd,
                     const matcl::Matrix& ind)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::schur_impl<V,S>::eval_cond_vec(ret, mat, sd, ind);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, Matrix& ret, const schur_decomposition& sd,
                     const matcl::Matrix& ind)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, ret, sd, ind);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const schur_decomposition&,
                     const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, Matrix&, const schur_decomposition&,
                     const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};
struct estim_rcond_vec_diag_vis : public extract_type_switch<void, estim_rcond_vec_diag_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& D, Matrix& ret, const Matrix& ind)
    {
        using V         = typename T::value_type;
        using Mat       = raw::Matrix<V,struct_dense>;
        const Mat& DD   = raw::converter<Mat,T>::eval(D);
        estim_rcond_vec_diag<V>::eval(ret, DD, ind);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& scal, Matrix& ret, const Matrix&)
    {
        ret = abs(scal);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const Matrix&)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};
struct estim_rcond_subspace_diag_vis : public extract_type_switch<Real, estim_rcond_subspace_diag_vis,true>
{
    template<class T>
    static Real eval(const Matrix&, const T& D, Integer M)
    {
        using V         = typename T::value_type;
        using Mat       = raw::Matrix<V,struct_dense>;
        const Mat& DD   = raw::converter<Mat,T>::eval(D);
        return estim_rcond_vec_diag<V>::eval_subspace(DD, M);
    };

    template<class T>
    static Real eval_scalar(const Matrix&, const T& scal, Integer)
    {
        return abs(scal);
    };

    static Real eval_scalar(const Matrix&, const Object&, Integer)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    };
    template<class S>
    static Real eval(const Matrix&, const raw::Matrix<Object,S>&, Integer)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};

struct schur_visitor_rcond_eig_cluster : public extract_type_switch<Real, schur_visitor_rcond_eig_cluster,true>
{
    template<class T>
    static Real eval(const Matrix&, const T& TA, Integer M)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::schur_impl<V,S>::eval_cond_eig_cluster(TA, M);
    };

    template<class T>
    static Real eval_scalar(const Matrix& h, const T& scal, Integer M)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, M);
    };

    static Real eval_scalar(const Matrix&, const Object&, Integer)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    };
    template<class S>
    static Real eval(const Matrix&, const raw::Matrix<Object,S>&, Integer)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};

struct schur_visitor_rcond_subspace : public extract_type_switch<Real, schur_visitor_rcond_subspace,true>
{
    template<class T>
    static Real eval(const Matrix&, const T& TA, Integer M)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::schur_impl<V,S>::eval_cond_subspace(TA, M);
    };

    template<class T>
    static Real eval_scalar(const Matrix& h, const T& scal, Integer M)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, M);
    };

    static Real eval_scalar(const Matrix&, const Object&, Integer)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    };
    template<class S>
    static Real eval(const Matrix&, const raw::Matrix<Object,S>&, Integer)
    {
        throw error::object_value_type_not_allowed("schur_decomposition");
    }
};

}; // end namespace details

schur_decomposition::schur_decomposition()
{
    clear();
};
void schur_decomposition::clear()
{
    m_has_U     = false;
    m_is_nan    = false;
    m_U_factor  = zeros(0,0);
    m_TA        = m_U_factor;
    m_eig       = zeros(0,1);
};
schur_decomposition::schur_decomposition(const Matrix &A0, schur_sym_alg alg, bool with_U)
{
    Matrix A(A0);

    clear();
    compute(A, alg, with_U);
};
schur_decomposition::schur_decomposition(Matrix && A0, schur_sym_alg alg, bool with_U)
{
    Matrix A(std::move(A0));

    clear();
    compute(A, alg, with_U);
};

schur_decomposition& schur_decomposition::operator()(const Matrix &A0, schur_sym_alg alg, bool with_U)
{
    Matrix A(A0);

    clear();
    compute(A, alg, with_U);
    return *this;
};
schur_decomposition& schur_decomposition::operator()(Matrix &&A0, schur_sym_alg alg, bool with_U)
{
    Matrix A(std::move(A0));

    clear();
    compute(A, alg, with_U);
    return *this;
};
schur_decomposition::schur_decomposition(const Matrix &U, const Matrix &T)
{
    clear();
    set(U, T, true);
};
schur_decomposition& schur_decomposition::set_factors(const Matrix &U, const Matrix &T)
{
    clear();
    set(U, T, true);
    return *this;
};
schur_decomposition& schur_decomposition::set_factors(const Matrix &T)
{
    clear();
    set(1.0, T, false);
    return *this;
};
void schur_decomposition::set(const Matrix &U, const Matrix &T, bool with_U)
{
    if (!T.is_square() )
        throw error::error_schur_incorrect_factors();

    if (with_U == true)
    {
        if (U.cols() != T.rows() || U.rows() < T.rows() )
            throw error::error_schur_incorrect_factors();
    }

    Integer ld  = get_ld(T,1);

    if (ld > 1)
    {
        //S is not quasi upper triangular
        throw error::error_schur_incorrect_factors();
    };

    bool compl  =  with_U && matrix_traits::is_float_complex(U.get_value_code())
                || matrix_traits::is_float_complex(T.get_value_code());

    bool obj    =  (with_U && U.get_value_code() == value_code::v_object)
                || (T.get_value_code() == value_code::v_object);
    
    if (compl == true && ld > 0)
        throw error::error_schur_incorrect_factors();

    if (obj == true)
        throw error::object_value_type_not_allowed("schur_decomposition");

    m_U_factor  = U;
    m_TA        = T;
    m_has_U     = with_U;

    if (test_factors() == false)
        return;

    set_eig();
}

bool schur_decomposition::test_factors()
{    
    bool isv    = m_TA.all_finite();

    if (m_has_U == true)
        isv     = isv && m_U_factor.all_finite();

    Integer N   = m_TA.rows();

    if (isv == false)
    {
        m_TA        = details::make_nan_matrix(N, N, m_TA.get_value_code());        
        m_eig       = details::make_nan_matrix(N, 1, m_TA.get_value_code());

        if (m_has_U)
            m_U_factor  = details::make_nan_matrix(m_U_factor.rows(), m_U_factor.cols(), 
                                                   m_U_factor.get_value_code());

        m_is_nan    = true;
        return false;
    }
    else
    {
        m_is_nan    = false;
    };
    
    return true;
};

schur_decomposition::~schur_decomposition()
{};

Matrix schur_decomposition::U() const
{
    if (m_has_U)
        return m_U_factor;
    else
        throw error::error_schur_U_not_computed();
};
Matrix schur_decomposition::TA() const
{
    return m_TA;
};

Matrix schur_decomposition::eig() const
{
    return m_eig;
};

void schur_decomposition::compute(const Matrix &A, schur_sym_alg alg, bool with_U)
{
    if (!A.is_square())
        throw error::error_size_eig(A.rows(),A.cols());

    m_has_U             = with_U;

    bool is_complex     = matrix_traits::is_float_complex(A.get_value_code());
    bool is_obj         = (A.get_value_code() == value_code::v_object);

    if (is_obj == true)
        throw error::object_value_type_not_allowed("schur_decomposition");

    Integer N           = A.rows();
       
    bool isv            = A.all_finite();
    
    if (isv == false)
    {
        m_TA            = details::make_nan_matrix(N, N, A.get_value_code());
        
        if (with_U)
            m_U_factor  = details::make_nan_matrix(N, N, A.get_value_code());

        m_eig           = details::make_nan_matrix(N, 1, A.get_value_code());
        m_is_nan        = true;        
        return;
    }
    else
    {
        m_is_nan        = false;
    };

    bool A_is_triu      = is_triu(A);
    bool A_is_qtriu     = A.get_struct().is_qtriu() && (is_complex == false);
    bool A_is_tril      = is_tril(A);
    bool A_is_qtril     = A.get_struct().is_qtril() && (is_complex == false);
    
    if (A.structural_nnz() == 0 || A_is_triu || A_is_qtriu)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        // remember that qtriu matrices need not be in the schur canonical form,
        // this problem is fixed later in qtriu2schur function called by set_eig()
        if (with_U)
            m_U_factor  = speye(N,N, vt);

        m_TA    = A;

        if (A_is_triu)
            m_TA.add_struct(predefined_struct_type::triu);

        if(A_is_triu && A_is_tril)
            m_TA.add_struct(predefined_struct_type::diag);

        set_eig();
        
        return;
    }
    
    if (A_is_tril || A_is_qtril)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        if (with_U)
        {
            Matrix I    = speye(N,N, vt);
            m_U_factor  = I(colon(), colon(end, -1, 1)); // reverse eye

            struct_flag sf;
            sf.set_user(unitary_flag());
            m_U_factor.add_struct(sf);
        };

        m_TA        = A(colon(end, -1, 1),colon(end, -1, 1));

        if (A_is_tril)
            m_TA.add_struct(predefined_struct_type::triu);
        else
            m_TA.add_struct(predefined_struct_type::qtriu);
        
        set_eig();
        return;
    }

    details::schur_visitor::make<const Matrix&>(A,*this, alg, with_U);
};

void schur_decomposition::set_eig()
{
    if (m_TA.structural_nnz() == 0 || is_triu(m_TA))
    {
        m_eig = get_diag(m_TA);
        return;
    }

    matcl::value_code vTA   = m_TA.get_value_code();
    matcl::value_code v     = matrix_traits::unify_value_types(vTA, value_code::v_float);

    if (m_has_U)
    {        
        matcl::value_code vQ    = m_U_factor.get_value_code();
        v                       = matrix_traits::unify_value_types(v, vQ);
    };    

    if (m_TA.get_value_code() != v)
    {
        matcl::mat_code mc      = matrix_traits::get_matrix_type(v, m_TA.get_struct_code());
        m_TA                    = matcl::convert(m_TA, mc);
    };

    if (m_has_U && m_U_factor.get_value_code() != v)
    {
        matcl::mat_code mc      = matrix_traits::get_matrix_type(v, m_U_factor.get_struct_code());
        m_U_factor              = matcl::convert(m_U_factor, mc);
    };

    return select_val(zeros(m_TA.rows(),1), true);
};

void schur_decomposition::select(const Matrix& ind)
{    
    Integer M;
    if (m_TA.structural_nnz() == 0 || m_TA.get_struct().is_id())
    {
        details::schur_reorder_check(ind, m_TA.rows(),M);
        return;
    }

    if (m_is_nan == true)
    {
        details::schur_reorder_check(ind, m_TA.rows(),M);
        return;
    };

    return select_val(ind, false);
};

void schur_decomposition::select_val(const Matrix& ind, bool no_fast_exit)
{
    Integer N   = m_TA.rows();
    Integer M   = 0;

    Matrix I = convert(ind,mat_code::integer_dense);
    details::schur_reorder_check(I,N,M);
    
    // fast exit for trivial ind vectors like [0 0 0 0] [1 1 1] [1 1 0 0 0] etc.
    if (!no_fast_exit && details::schur_is_trivial_reorder(I))
        return;

    return details::schur_visitor_reorder::make<const Matrix&>(m_TA,ind,*this);
};

Matrix schur_decomposition::left_eigenvectors() const
{
    Integer N = m_TA.rows();
    return left_eigenvectors(iones(N,1));
};
Matrix schur_decomposition::right_eigenvectors() const
{
    Integer N = m_TA.rows();
    return right_eigenvectors(iones(N,1));
};

tuple<Matrix,Matrix> schur_decomposition::eigenvectors() const
{
    Integer N = m_TA.rows();
    return eigenvectors(iones(N,1));
};
Matrix schur_decomposition::left_eigenvectors(const Matrix& ind0) const
{
    //increase refcount
    Matrix ind(ind0);

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, true, false);
    return XL;
};
Matrix schur_decomposition::left_eigenvectors(Matrix&& ind0) const
{
    //increase refcount
    Matrix ind(std::move(ind0));

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, true, false);
    return XL;
};

Matrix schur_decomposition::right_eigenvectors(const Matrix& ind0) const
{
    //increase refcount
    Matrix ind(ind0);

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, false, true);
    return XR;
};
Matrix schur_decomposition::right_eigenvectors(Matrix&& ind0) const
{
    //increase refcount
    Matrix ind(std::move(ind0));

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, false, true);
    return XR;
};

tuple<Matrix,Matrix> schur_decomposition::eigenvectors(const Matrix& ind0) const
{
    //increase refcount
    Matrix ind(ind0);

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, true, true);
    return tuple<Matrix,Matrix>(XL,XR);
};
tuple<Matrix,Matrix> schur_decomposition::eigenvectors(Matrix&& ind0) const
{
    //increase refcount
    Matrix ind(std::move(ind0));

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, true, true);
    return tuple<Matrix,Matrix>(XL,XR);
};
void schur_decomposition::comp_eigenvectors(Matrix& ind, Matrix& XL, Matrix& XR, 
                            bool comp_left, bool comp_right) const
{
    Integer N = m_TA.rows();
    Integer M;
    details::schur_reorder_check(ind, m_TA.rows(),M);

    if (m_is_nan)
    {
        XL  = details::make_nan_matrix(N, M, m_TA.get_value_code());
        XR  = XL;
        return;
    };

    // fast exit
    if (M == 0)
    {
        XL = zeros(N,0,m_TA.get_value_code());
        XR = XL;
        return;
    };

    // T is diagonal, eigenvectors are orthonormal
    if (is_diag(m_TA) == true)
    {
        Matrix I = find(ind == 1);

        if (m_has_U)
        {
            XL  = m_U_factor(colon(), I);
            XR  = XL;
            return;
        };

        Matrix U = speye(N,N,m_TA.get_value_code());
        XL      = U(colon(), I);
        XR      = XL;
        return;
    };

    return details::schur_visitor_eigvec
        ::make<const Matrix&>(m_TA,ind,*this, XL, XR, comp_left, comp_right);
};

Matrix schur_decomposition::rcond_eig(const Matrix& VL, const Matrix& VR) const
{
    Integer N = m_TA.rows();
    return rcond_eig(VL, VR, iones(N,1));
};
Matrix schur_decomposition::rcond_vec() const
{
    Integer N = m_TA.rows();
    return rcond_vec(iones(N,1));
};

Matrix schur_decomposition::rcond_eig() const
{
    Matrix VL, VR;
    tie(VL, VR) = eigenvectors();
    return rcond_eig(VL, VR);
};
Matrix schur_decomposition::rcond_eig(const Matrix& ind) const
{
    Matrix VL, VR;
    tie(VL, VR) = eigenvectors(ind);
    return rcond_eig(VL, VR, ind);
};
Matrix schur_decomposition::rcond_eig(const Matrix& VL, const Matrix& VR, const Matrix& ind) const
{
    Integer N = m_TA.rows();
    Integer M;
    details::schur_reorder_check(ind, m_TA.rows(),M);

    if (VL.rows() != N || VL.cols() != M)
        throw error::error_schur_cond_vec(N, M, VL.rows(), VL.cols(), true);

    if (VR.rows() != N || VR.cols() != M)
        throw error::error_schur_cond_vec(N, M, VR.rows(), VR.cols(), false);

    if (m_is_nan)
        return details::make_nan_matrix(M, 1, m_TA.get_value_code());

    // fast exit
    if (M == 0)
    {
        matcl::value_code vt = matrix_traits::real_value_type(m_TA.get_value_code());
        return zeros(0,1, vt);
    };

    if (is_diag(m_TA) == true)
    {
        matcl::value_code vt = matrix_traits::real_value_type(m_TA.get_value_code());
        return ones(M,1, vt);
    };

    Matrix ret;
    details::schur_visitor_cond_eig::make<const Matrix&>(m_TA, ret, *this, VL, VR, ind);
    return ret;

};
Matrix schur_decomposition::rcond_vec(const Matrix& ind) const
{
    Integer M;
    details::schur_reorder_check(ind, m_TA.rows(),M);

    if (m_is_nan)
        return details::make_nan_matrix(M, 1, m_TA.get_value_code());

    // fast exit
    if (M == 0)
    {
        matcl::value_code vt = matrix_traits::real_value_type(m_TA.get_value_code());
        return zeros(0,1, vt);
    };

    // T is diagonal, eigenvectors are orthonormal
    if (is_diag(m_TA) == true)
    {
        Matrix D    = m_TA.diag();
        Matrix ret;
        details::estim_rcond_vec_diag_vis::make<const Matrix&>(D, ret, ind);
        return ret;
    };

    Matrix ret;
    details::schur_visitor_cond_vec::make<const Matrix&>(m_TA, ret, *this, ind);
    return ret;
};

Real schur_decomposition::rcond_projector(Integer M) const
{
    Integer N = m_TA.rows();

    if (M < 0 || M > N)
        throw error::error_schur_eig_cluster_size(M, N);

    if (m_is_nan)
        return constants::nan();

    if (is_diag(m_TA) == true)
        return 1.0;

    return details::schur_visitor_rcond_eig_cluster::make<const Matrix&>(m_TA, M);
};
Real schur_decomposition::separation(Integer M) const
{
    Integer N = m_TA.rows();

    if (M < 0 || M > N)
        throw error::error_schur_eig_cluster_size(M, N);

    if (m_is_nan)
        return constants::nan();

    if (is_diag(m_TA) == true)
    {
        Matrix D    = m_TA.diag();
        return details::estim_rcond_subspace_diag_vis::make<const Matrix&>(D, M);
    };

    return details::schur_visitor_rcond_subspace::make<const Matrix&>(m_TA, M);
};

Real schur_decomposition::pert_bound(Integer M) const
{
    Real p      = rcond_projector(M);
    Real SEP    = separation(M);

    //taken from
    // B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
    // Eigenvalues of a Regular Matrix Pair (A, B) and Condition
    // Estimation: Theory, Algorithms and Software,
    // Report UMINF - 94.04, Department of Computing Science, Umea
    // University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working
    // Note 87.

    if (p == 0.0)
        return 0.0;

    Real ret    = SEP / p;
    return ret;
};

linsolve_obj matcl::linsolve_schur(const Matrix& A, const Matrix& U, const Matrix& T, const options& opts)
{
    return linsolve_schur(A, unitary_matrix(U,true), T, opts);
};
linsolve_obj matcl::linsolve_schur(const Matrix& A, const unitary_matrix& U, const Matrix& T, 
                                   const options& opts)
{
    if (A.rows() != U.rows())
        throw error::invalid_schur_factors();
    if (U.cols() != T.rows())
        throw error::invalid_schur_factors();
    if (U.rows() < U.cols())
        throw error::invalid_schur_factors();

    if (matcl::get_ld(T,1) > 1)
        throw error::invalid_schur_factors();

    if (T.rows() != T.cols())
        throw error::square_matrix_required(T.rows(), T.cols());

    if (U.rows() != U.cols())
        throw error::square_matrix_required(U.rows(), U.cols());

    Integer N   = U.rows();

    if (U.rows() == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(T.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, U.get_type())));
    };

    bool isv            = T.all_finite() && U.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(T.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, T.get_type() )));
    };

    using data_ptr  = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_schur(A, U, T, opts)));
};

};

#pragma warning( pop )