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

#include "matcl-linalg/decompositions/schur.h"
#include "matcl-linalg/decompositions/gschur.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/decompositions/eig/schur_utils.h"
#include "matcl-linalg/decompositions/qr.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-core/utils/workspace.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/details/struct_flag_predefined.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

template<class Val>
struct gschur_str
{
    using VC    = typename md::complex_type<Val>::type;

    using Mat   = raw::Matrix<Val,struct_dense>;
    using Mat_C = raw::Matrix<VC,struct_dense>;

    static void eval(const Mat& A, const Mat& B, gschur_decomposition& gd, bool with_QZ)
    {
        using matcl::details::lap;

        //test for nan should already be performed

        const char* jobvsl  = with_QZ? "V" : "N";
        const char* jobvsr  = with_QZ? "V" : "N";

        Integer N           = A.rows();
        Integer Nqz         = with_QZ? N : 0;        

        Mat Q(A.get_type(), Nqz, Nqz);
        Mat Z(A.get_type(), Nqz, Nqz);
    
        Mat TA  = A.make_unique();
        Mat TB  = B.make_unique();
    
        TA.set_struct(struct_flag());
        TB.set_struct(struct_flag());

        Mat_C   alpha(A.get_type(), N, 1);
        Mat     beta(A.get_type(), N, 1);

        Integer Sdim;
        Integer info    = 0;
               
        Val work_query;
        Integer iwork_query = 0;

        //TODO
        lapack::gges(jobvsl, jobvsr, "N", nullptr, N, lap(TA.ptr()), TA.ld(), lap(TB.ptr()), TB.ld(), &Sdim,
                    lap(alpha.ptr()), lap(beta.ptr()), lap(Q.ptr()), Q.ld(), lap(Z.ptr()), Z.ld(), 
                    lap(&work_query), -1, nullptr, &info);
        //(void)dohess;
        
        //lapack::gges2(jobvsl, jobvsr, N, lap(TA.ptr()), TA.ld(), lap(TB.ptr()), TB.ld(), 
        //            lap(alpha.ptr()), lap(beta.ptr()), lap(Q.ptr()), Q.ld(), lap(Z.ptr()), Z.ld(), 
        //            lap(&work_query), -1, &iwork_query, -1, info);
        (void)Sdim;

        Integer lwork       = (Integer) real(work_query);
        Integer liwork      = iwork_query;

        using VTR_pod       = matcl::pod_type<Val>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        workspace WORK      = workspace(lwork);
        iworkspace IWORK    = iworkspace(liwork);
        Val* ptr_WORK       = reinterpret_cast<Val*>(WORK.ptr());
        Integer* ptr_IWORK  = IWORK.ptr();

        //TODO
        lapack::gges(jobvsl, jobvsr, "N", nullptr, N, lap(TA.ptr()), TA.ld(), lap(TB.ptr()), 
                TB.ld(), &Sdim, lap(alpha.ptr()), lap(beta.ptr()), lap(Q.ptr()), Q.ld(), lap(Z.ptr()),
                Z.ld(), lap(ptr_WORK), lwork, nullptr, &info);

        //lapack::gges2(jobvsl, jobvsr, N, lap(TA.ptr()), TA.ld(), lap(TB.ptr()), TB.ld(), 
        //            lap(alpha.ptr()), lap(beta.ptr()), lap(Q.ptr()), Q.ld(), lap(Z.ptr()), Z.ld(), lap(ptr_WORK),
        //            lwork, ptr_IWORK, liwork, info);
        (void)ptr_IWORK;

        bool do_not_set_factors = false;

        if (info)
        {
            // NON-converging case - try to fill sparse up matrices
            if (info > 0 && info <= N)
            {
                // produce random full unitary matrix
                // use it to fill in Am and Bm which just failed gges
                Matrix P;
                {
                    matcl::rand_state rsp = matcl::get_rand_state();
                    unitary_matrix UC;
                    init_genrand(123456789);

                    P       = rand_unitary(N).to_matrix();
                   
                    set_rand_state(rsp);
                }

                Matrix Am   = P * Matrix(A,false) * P;
                Matrix Bm   = P * Matrix(B,false) * P;

                Mat TA2     = convert(Am, Mat::matrix_code).get_impl_unique<Mat>();
                Mat TB2     = convert(Bm, Mat::matrix_code).get_impl_unique<Mat>();                

                lapack::gges(jobvsl, jobvsr, "N", nullptr, N, lap(TA2.ptr()), TA2.ld(), lap(TB2.ptr()), 
                        TB2.ld(), &Sdim, lap(alpha.ptr()), lap(beta.ptr()), lap(Q.ptr()), Q.ld(), lap(Z.ptr()),
                        Z.ld(), lap(ptr_WORK), lwork, nullptr, &info);

                //lapack::gges2(jobvsl, jobvsr, N, lap(TA2.ptr()), TA2.ld(), lap(TB2.ptr()), 
                //        TB2.ld(), lap(alpha.ptr()), lap(beta.ptr()), lap(Q.ptr()), Q.ld(), lap(Z.ptr()),
                //        Z.ld(), lap(ptr_WORK), lwork, ptr_IWORK, liwork, info);
            
                // retrieve factors
                do_not_set_factors = true;
            
                gd.m_TA     = Matrix(TA2,true);
                gd.m_TB     = Matrix(TB2,true);

                if (with_QZ)
                {
                    gd.m_Q_factor  = mmul(P, Matrix(Q,true), trans_type::conj_trans);
                    gd.m_Z_factor  = mmul(P, Matrix(Z,true));
                };
            
                if (info) // give up
                {
                    gd.clear();
                    throw error::error_gschur();
                }
            }
            else
            {
                gd.clear();
                throw error::error_gschur();
            }
        };
    
        // if this was an info > 0 call to gges, these matrices have alreade been set!!
        if (!do_not_set_factors)
        {
            gd.m_TA     = Matrix(TA,true);
            gd.m_TB     = Matrix(TB,true);

            if (with_QZ)
            {
                gd.m_Q_factor  = Matrix(Q,true);
                gd.m_Z_factor  = Matrix(Z,true);
            };
        }

        if (gd.test_factors() == false)
            return;

        gd.m_alpha  = Matrix(alpha,true);
        gd.m_beta   = Matrix(beta,true);
        gd.m_eig    = gd.comp_eig(gd.m_alpha, gd.m_beta);

        if (md::is_complex<Val>::value == false)
        {
            Integer ld  = get_ld(gd.m_TA, 1);

            if (ld == 1)
                gd.m_TA.add_struct(md::predefined_struct::get_qtriu(gd.m_TA.get_struct()));
            else
                gd.m_TA.add_struct(md::predefined_struct::get_triu(gd.m_TA.get_struct(),0, 
                                                                   is_real_matrix(gd.m_TA)));
        }
        else
        {
            gd.m_TA.add_struct(md::predefined_struct::get_triu(gd.m_TA.get_struct(),0, 
                                                               is_real_matrix(gd.m_TA)));
        };

        gd.m_TB.add_struct(md::predefined_struct::get_triu(gd.m_TB.get_struct(),0, is_real_matrix(gd.m_TB)));

        if (with_QZ)
        {
            struct_flag sf_u;
            sf_u.set_user(unitary_flag());
            gd.m_Q_factor.add_struct(sf_u);
            gd.m_Z_factor.add_struct(sf_u);
        };        
    };

    static void eval_eigenvec(const Mat& TA0, matcl::Matrix& ind, const gschur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        (void)TA0;

        const char* SIDE        = (comp_left && comp_right)? "B" : (comp_left? "L" : "R");
        Integer* SELECT         = ind.get_array_unique<Integer>();

        const Mat& TA           = convert(sd.m_TA, Mat::matrix_code).get_impl<Mat>();
        const Mat& TB           = convert(sd.m_TB, Mat::matrix_code).get_impl<Mat>();

        Integer N               = sd.m_TA.rows();
        const Val* ptr_TA       = TA.ptr();
        Integer LD_TA           = TA.ld();
        const Val* ptr_TB       = TB.ptr();
        Integer LD_TB           = TB.ld();

        Integer M               = make_complex_eigenvectors<Val>::calc_sel_size(TA, SELECT, N);
        Integer ML              = comp_left? M : 1;
        Integer MR              = comp_right? M : 1;

        Mat VL                  = Mat(TA.get_type(), N, ML);
        Mat VR                  = Mat(TA.get_type(), N, MR);

        Integer info;

        using VTR_pod           = matcl::pod_type<Val>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        Integer lwork           = 6*N;
        workspace WORK          = workspace(lwork);
        Val* ptr_WORK           = reinterpret_cast<Val*>(WORK.ptr());

        lapack::tgevc(SIDE, "S", SELECT, N, lap(ptr_TA), LD_TA, lap(ptr_TB), LD_TB, lap(VL.ptr()), VL.ld(),
                      lap(VR.ptr()), VR.ld(), M, M, lap(ptr_WORK), &info);

        if (info != 0)
            throw error::error_general("illegal argument passed to trevc");

        Matrix mat_VL           = Matrix(VL,false);
        Matrix mat_VR           = Matrix(VR,false);

        if (sd.m_has_QZ)
        {
            // form U*VL, Z*VR
            if (comp_left)
                mat_VL          = mmul(sd.m_Q_factor, mat_VL, trans_type::no_trans);

            if (comp_right)
                mat_VR          = mmul(sd.m_Z_factor, mat_VR, trans_type::no_trans);
        };

        make_complex_eigenvectors<Val>::eval(TA, ind, mat_VL, mat_VR, XL, XR, comp_left, comp_right);
        return;
    };

    static void eval_cond_eig(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& WL, const Matrix& WR)
    {
        (void)A;

        using VR                = typename md::real_type<Val>::type;
        using Mat_R             = raw::Matrix<VR,struct_dense>;

        const char* JOB         = "E";
        const char* HOWMNY      = "S";
        const Integer* SELECT   = ind.get_array<Integer>();

        const Mat& TA           = convert(sd.m_TA, Mat::matrix_code).get_impl<Mat>();
        const Mat& TB           = convert(sd.m_TB, Mat::matrix_code).get_impl<Mat>();

        Integer N               = TA.rows();
        Integer M               = make_complex_eigenvectors<Val>::calc_sel_size(TA, SELECT, N);

        const Val* ptr_TA       = TA.ptr();
        Integer LD_TA           = TA.ld();
        const Val* ptr_TB       = TB.ptr();
        Integer LD_TB           = TB.ld();

        Mat_R S(TA.get_type(), M, 1);
        VR* ptr_S               = S.ptr();

        Mat WL_R(TA.get_type());
        Mat WR_R(TA.get_type());

        make_complex_eigenvectors<Val>::complex_vectors_to_real(TA, WL, WR, WL_R, WR_R, ind);

        Val* ptr_VL             = WL_R.ptr();
        Val* ptr_VR             = WR_R.ptr();
        Integer ld_VL           = WL_R.ld();
        Integer ld_VR           = WR_R.ld();

        VR* ptr_null_r          = nullptr;
        Integer* ptr_null_i     = nullptr;

        Val work_query;
        Integer info;        

        lapack::tgsna(JOB, HOWMNY, SELECT, N, lap(ptr_TA), LD_TA, lap(ptr_TB), LD_TB, lap(ptr_VL), ld_VL, 
                      lap(ptr_VR), ld_VR, lap(ptr_S), lap(ptr_null_r), M, M, 
                      lap(&work_query), -1, ptr_null_i, &info);

        Integer lwork           = (Integer)real(work_query);

        using VTR_pod           = matcl::pod_type<Val>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        workspace WORK          = workspace(lwork);
        Val* ptr_WORK           = reinterpret_cast<Val*>(WORK.ptr());

        lapack::tgsna(JOB, HOWMNY, SELECT, N, lap(ptr_TA), LD_TA, lap(ptr_TB), LD_TB, lap(ptr_VL), ld_VL, 
                      lap(ptr_VR), ld_VR, lap(ptr_S), lap(ptr_null_r), M, M, 
                      lap(ptr_WORK), lwork, ptr_null_i, &info);

        if (info != 0)
            throw error::error_general("illegal argument passed to tgsna");

        ret = matcl::Matrix(S,false);
        return;
    };

    static void eval_cond_vec(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind)
    {
        (void)A;

        using VR                = typename md::real_type<Val>::type;
        using Mat_R             = raw::Matrix<VR,struct_dense>;

        const char* JOB         = "V";
        const char* HOWMNY      = "S";
        const Integer* SELECT   = ind.get_array<Integer>();

        const Mat& TA           = convert(sd.m_TA, Mat::matrix_code).get_impl<Mat>();
        const Mat& TB           = convert(sd.m_TB, Mat::matrix_code).get_impl<Mat>();

        Integer N               = TA.rows();
        Integer M               = make_complex_eigenvectors<Val>::calc_sel_size(TA, SELECT, N);

        const Val* ptr_TA       = TA.ptr();
        Integer LD_TA           = TA.ld();
        const Val* ptr_TB       = TB.ptr();
        Integer LD_TB           = TB.ld();

        Mat_R SEP(TA.get_type(), M, 1);
        VR* ptr_SEP             = SEP.ptr();

        Val* ptr_null           = nullptr;
        VR* ptr_null_r          = nullptr;

        using iworkspace        = matcl::pod_workspace<Integer>;
        Integer liwork          = N + 6;
        iworkspace IWORK        = iworkspace(liwork);
        Integer* ptr_IWORK      = IWORK.ptr();

        Val work_query;
        Integer info;

        lapack::tgsna( JOB, HOWMNY, SELECT, N, lap(ptr_TA), LD_TA, lap(ptr_TB), LD_TB, lap(ptr_null), 1, 
                      lap(ptr_null), 1, lap(ptr_null_r), lap(ptr_SEP), M, M, lap(&work_query), -1, ptr_IWORK, &info);

        Integer lwork           = (Integer)real(work_query);

        using VTR_pod           = matcl::pod_type<Val>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        workspace WORK          = workspace(lwork);
        Val* ptr_WORK           = reinterpret_cast<Val*>(WORK.ptr());

        lapack::tgsna( JOB, HOWMNY, SELECT, N, lap(ptr_TA), LD_TA, lap(ptr_TB), LD_TB, lap(ptr_null), 1, 
                      lap(ptr_null), 1, lap(ptr_null_r), lap(ptr_SEP), M, M, lap(ptr_WORK), lwork, ptr_IWORK, &info);

        if (info != 0)
            throw error::error_general("illegal argument passed to tgsna");

        ret = matcl::Matrix(SEP,false);
        return;
    };
    static void eval_rcond(const Mat& A, const gschur_decomposition& sd, Integer M, Integer type, 
                           Real& ret1, Real& ret2)
    {
        (void)A;
        using VR        = typename md::real_type<Val>::type;

        Integer IJOB    = type; 

        const Mat& TA   = convert(sd.m_TA, Mat::matrix_code).get_impl<Mat>();
        const Mat& TB   = convert(sd.m_TB, Mat::matrix_code).get_impl<Mat>();

        Integer N           = TA.rows();
        const Val* ptr_TA   = TA.ptr();
        const Val* ptr_TB   = TB.ptr();

        Integer ld_TA   = TA.ld();
        Integer ld_TB   = TB.ld();

        //check if M-th eigenvalue is complex
        if (M > 0 && M < N)
        {
            if (ptr_TA[M + (M-1)*ld_TA] != Val(0.0))
                M       = M + 1;
        };

        VR PL, PR;
        VR DIF[2];

        Val     work_query;
        Integer iwork_query;
        Integer info;

        lapack::tgsen_cond(IJOB, N, M, lap(ptr_TA), ld_TA, lap(ptr_TB), ld_TB, PL, PR, DIF, 
                           lap(&work_query), -1, &iwork_query, -1, info);

        Integer lwork   = (Integer)real(work_query);
        Integer liwork  = iwork_query;

        using VTR_pod           = matcl::pod_type<Val>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        workspace WORK          = workspace(lwork);
        Val* ptr_WORK           = reinterpret_cast<Val*>(WORK.ptr());

        using iworkspace        = matcl::pod_workspace<Integer>;
        iworkspace IWORK        = iworkspace(liwork);
        Integer* ptr_IWORK      = IWORK.ptr();

        lapack::tgsen_cond(IJOB, N, M, lap(ptr_TA), ld_TA, lap(ptr_TB), ld_TB, PL, PR, DIF, 
                           lap(ptr_WORK), lwork, ptr_IWORK, liwork, info);

        if (info != 0)
            throw error::error_general("illegal argument passed to tgsen_cond");

        if (type == 1)
        {
            ret1    = PL;
            ret2    = PR;
        }
        else
        {
            ret1    = DIF[0];
            ret2    = DIF[1];
        }

        return;
    };
};

template<class V, class S>
struct gschur2_str
{
    using Mat   = raw::Matrix<V,S>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval_eigenvec(const Mat& TA, matcl::Matrix& ind, const gschur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(TA);
        return gschur2_str<V,struct_dense>::eval_eigenvec(Ac,ind,sd,XL,XR,comp_left,comp_right);
    };
    static void eval_cond_eig(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& VL, const Matrix& VR)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return gschur2_str<V,struct_dense>::eval_cond_eig(ret, Ac,sd, ind, VL, VR);
    };
    static void eval_cond_vec(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return gschur2_str<V,struct_dense>::eval_cond_vec(ret, Ac,sd, ind);
    };    
    static void eval_rcond(const Mat& A, const gschur_decomposition& sd, Integer M, Integer type, 
                           Real& ret1, Real& ret2)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return gschur2_str<V,struct_dense>::eval_rcond(Ac,sd,M,type,ret1,ret2);
    };

};
template<class V>
struct gschur2_str<V,struct_dense>
{
    using Mat   = raw::Matrix<V,struct_dense>;

    static void eval_eigenvec(const Mat& TA, matcl::Matrix& ind, const gschur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        return gschur_str<V>::eval_eigenvec(TA,ind,sd,XL,XR,comp_left,comp_right);
    };

    static void eval_cond_eig(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& WL, const Matrix& WR)
    {
        return gschur_str<V>::eval_cond_eig(ret, A, sd, ind, WL, WR);
    };

    static void eval_cond_vec(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind)
    {
        return gschur_str<V>::eval_cond_vec(ret, A, sd, ind);
    };
    static void eval_rcond(const Mat& A, const gschur_decomposition& sd, Integer M, Integer type, 
                           Real& ret1, Real& ret2)
    {
        return gschur_str<V>::eval_rcond(A,sd,M,type,ret1,ret2);
    };
};


template<class V1, class V2, class S1, class S2>
struct gschur_impl
{
    static_assert(std::is_same<V1,V2>::value == false 
                  || std::is_same<S1,struct_dense>::value == false
                  || std::is_same<S2,struct_dense>::value == false, "invalid template specialization");

    using Mat1  = raw::Matrix<V1,S1>;
    using Mat2  = raw::Matrix<V2,S2>;

    using VR    = typename md::unify_types<V1,V2>::type;
    using Mat_R = raw::Matrix<VR,struct_dense>;

    static void eval(const Mat1& A, const Mat2& B, gschur_decomposition& gd, bool with_QZ)
    {
        Mat_R Ac    = raw::converter<Mat_R, Mat1>::eval(A);
        Mat_R Bc    = raw::converter<Mat_R, Mat2>::eval(B);

        return gschur_impl<VR,VR,struct_dense,struct_dense>::eval(Ac,Bc,gd,with_QZ);
    };
};

template<class V1>
struct gschur_impl<V1,V1,struct_dense,struct_dense>
{
    using Mat   = raw::Matrix<V1,struct_dense>;

    static void eval(const Mat& A, const Mat& B, gschur_decomposition& gd, bool with_QZ)
    {
        return gschur_str<V1>::eval(A,B,gd,with_QZ);
    };
};

template<>
struct gschur_impl<Integer,Integer,struct_dense,struct_dense>
{
    using Mat   = raw::Matrix<Integer,struct_dense>;
    using Mat_R = raw::Matrix<Real,struct_dense>;

    static void eval(const Mat& A, const Mat& B, gschur_decomposition& gd, bool with_QZ)
    {
        Mat_R Ac    = raw::converter<Mat_R, Mat>::eval(A);
        Mat_R Bc    = raw::converter<Mat_R, Mat>::eval(B);

        return gschur_impl<Real,Real,struct_dense,struct_dense>::eval(Ac,Bc,gd,with_QZ);
    };
};
template<>
struct gschur_impl<Object,Object,struct_dense,struct_dense>
{
    using Mat   = raw::Matrix<Object,struct_dense>;

    static void eval(const Mat&, const Mat&, gschur_decomposition&, bool)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    };
};

template<class V, class S>
struct gschur2_impl
{
    using Mat   = raw::Matrix<V,S>;

    static void eval_eigenvec(const Mat& A, matcl::Matrix& ind, const gschur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        return gschur2_str<V,S>::eval_eigenvec(A,ind,sd,XL,XR,comp_left,comp_right);
    };
    static void eval_cond_eig(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& VL, const Matrix& VR)
    {
        return gschur2_str<V,S>::eval_cond_eig(ret,A,sd,ind,VL,VR);
    };
    static void eval_cond_vec(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind)
    {
        return gschur2_str<V,S>::eval_cond_vec(ret,A,sd,ind);
    };
    static void eval_rcond(const Mat& A, const gschur_decomposition& sd, Integer M, Integer type, 
                           Real& ret1, Real& ret2)
    {
        return gschur2_str<V,S>::eval_rcond(A,sd,M,type,ret1,ret2);
    };
};

template<class S>
struct gschur2_impl<Integer,S>
{
    using Mat   = raw::Matrix<Integer,S>;
    using Mat_D = raw::Matrix<Real,struct_dense>;

    static void eval_eigenvec(const Mat& A, matcl::Matrix& ind, const gschur_decomposition& sd, 
                              Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return gschur2_impl<Real,struct_dense>::eval_eigenvec(Ac,ind,sd,XL,XR,comp_left,comp_right);
    };
    static void eval_cond_eig(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind, 
                              const Matrix& VL, const Matrix& VR)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return gschur2_impl<Real,struct_dense>::eval_cond_eig(ret, Ac,sd, ind, VL, VR);
    };
    static void eval_cond_vec(Matrix& ret, const Mat& A, const gschur_decomposition& sd, const matcl::Matrix& ind)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return gschur2_impl<Real,struct_dense>::eval_cond_vec(ret, Ac,sd, ind);
    };    
    static void eval_rcond(const Mat& A, const gschur_decomposition& sd, Integer M, Integer type, 
                           Real& ret1, Real& ret2)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return gschur2_impl<Real,struct_dense>::eval_rcond(Ac,sd,M,type,ret1,ret2);
    };
};

struct gschur_vis : public extract_type2_switch<void,gschur_vis, mr::val_type_corrector_diag_dense>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, gschur_decomposition& gd, bool with_QZ)
    {
        using V1    = typename T1::value_type;
        using S1    = typename T1::struct_type;

        using V2    = typename T2::value_type;
        using S2    = typename T2::struct_type;

        return details::gschur_impl<V1,V2,S1,S2>::eval(A,B,gd,with_QZ);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, gschur_decomposition& gd, bool with_QZ)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        using Dense_2 = raw::Matrix<T2,struct_dense>;

        Dense_1 Ac(ti::get_ti(A),A,1,1);
        Dense_2 Bc(ti::get_ti(B),B,1,1);

        return eval_mat_mat<Dense_1,Dense_2>(Ac,Bc, gd, with_QZ);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& A, const T2& B, gschur_decomposition& gd, bool with_QZ)
    {
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<T1,Dense_2>(A,Bc, gd, with_QZ);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& A, const T2& B, gschur_decomposition& gd, bool with_QZ)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<Dense_1,T2>(Ac,B, gd, with_QZ);
    };
};

struct gschur_visitor_eigvec : public extract_type_switch<void, gschur_visitor_eigvec,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, matcl::Matrix& ind, const gschur_decomposition& sd,
                    Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::gschur2_impl<V,S>::eval_eigenvec(mat, ind, sd, XL, XR, comp_left, comp_right);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, matcl::Matrix& ind, const gschur_decomposition& sd,
                    Matrix& XL, Matrix& XR, bool comp_left, bool comp_right)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, ind, sd,XL,XR,comp_left,comp_right);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, matcl::Matrix&, const gschur_decomposition&,
                     Matrix&, Matrix&, bool, bool)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, matcl::Matrix&, const gschur_decomposition&,
                            Matrix&, Matrix&, bool, bool)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    }
};

struct gschur_visitor_cond_vec : public extract_type_switch<void, gschur_visitor_cond_vec,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, const gschur_decomposition& sd,
                     const matcl::Matrix& ind)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::gschur2_impl<V,S>::eval_cond_vec(ret, mat, sd, ind);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, Matrix& ret, const gschur_decomposition& sd,
                     const matcl::Matrix& ind)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, ret, sd, ind);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const gschur_decomposition&,
                     const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, Matrix&, const gschur_decomposition&,
                     const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    }
};
struct gschur_visitor_cond_eig : public extract_type_switch<void, gschur_visitor_cond_eig,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, const gschur_decomposition& sd,
                     const Matrix& VL, const Matrix& VR, const matcl::Matrix& ind)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::gschur2_impl<V,S>::eval_cond_eig(ret, mat, sd, ind, VL, VR);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, Matrix& ret, const gschur_decomposition& sd,
                     const Matrix& VL, const Matrix& VR, const matcl::Matrix& ind)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, ret, sd, VL, VR, ind);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const gschur_decomposition&,
                     const Matrix&, const Matrix&, const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, Matrix&, const gschur_decomposition&,
                     const Matrix&, const Matrix&, const matcl::Matrix&)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    }
};
struct gschur_visitor_rcond : public extract_type_switch<void, gschur_visitor_rcond,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, const gschur_decomposition& gd, 
                     Integer M, Integer type, Real& ret1, Real& ret2)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return details::gschur2_impl<V,S>::eval_rcond(mat, gd, M, type, ret1, ret2);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& scal, const gschur_decomposition& gd, 
                            Integer M, Integer type, Real& ret1, Real& ret2)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(scal), scal, 1, 1);

        return eval<Mat>(h, m, gd, M, type, ret1, ret2);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const gschur_decomposition&,
                     Integer, Integer, Real&, Real&)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    }

    static void eval_scalar(const Matrix&, const Object&, const gschur_decomposition&,
                            Integer, Integer, Real&, Real&)
    {
        throw error::object_value_type_not_allowed("gschur_decomposition");
    }
};

};};

namespace matcl
{

gschur_decomposition::gschur_decomposition()
{
    clear();
};

gschur_decomposition::~gschur_decomposition()
{};

void gschur_decomposition::clear()
{
    m_Q_factor  = zeros(0,0);
    m_Z_factor  = m_Q_factor;
    m_TA        = m_Q_factor;
    m_TB        = m_Q_factor;
    m_eig       = zeros(0,1);
    m_alpha     = m_eig;
    m_beta      = m_eig;
    m_has_QZ    = false;
    m_is_nan    = false;
};

Matrix gschur_decomposition::Q() const
{
    if (m_has_QZ)
        return m_Q_factor;
    else
        throw error::error_gschur_QZ_not_computed();
};
Matrix gschur_decomposition::Z() const
{
    if (m_has_QZ)
        return m_Z_factor;
    else
        throw error::error_gschur_QZ_not_computed();
};
Matrix gschur_decomposition::TA() const
{
    return m_TA;
};
Matrix gschur_decomposition::TB() const
{
    return m_TB;
};
Matrix gschur_decomposition::eig() const
{
    return m_eig;
};
Matrix gschur_decomposition::alpha() const
{
    return m_alpha;
};
Matrix gschur_decomposition::beta() const
{
    return m_beta;
};

gschur_decomposition::gschur_decomposition(const Matrix &A0, const Matrix &B0, bool with_QZ)
{
    Matrix A(A0);
    Matrix B(B0);

    clear();
    compute(A, B, with_QZ);
};
gschur_decomposition::gschur_decomposition(Matrix &&A0, const Matrix &B0, bool with_QZ)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    clear();
    compute(A, B, with_QZ);
};
gschur_decomposition::gschur_decomposition(const Matrix &A0, Matrix &&B0, bool with_QZ)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    clear();
    compute(A, B, with_QZ);
};
gschur_decomposition::gschur_decomposition(Matrix &&A0, Matrix && B0, bool with_QZ)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    clear();
    compute(A, B, with_QZ);
};

gschur_decomposition& gschur_decomposition::operator()(const Matrix &A0, const Matrix &B0, bool with_QZ)
{
    Matrix A(A0);
    Matrix B(B0);

    clear();
    compute(A, B, with_QZ);
    return *this;
};
gschur_decomposition& gschur_decomposition::operator()(Matrix &&A0, const Matrix &B0, bool with_QZ)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    clear();
    compute(A, B, with_QZ);
    return *this;
};
gschur_decomposition& gschur_decomposition::operator()(const Matrix &A0, Matrix &&B0, bool with_QZ)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    clear();
    compute(A, B, with_QZ);
    return *this;
};
gschur_decomposition& gschur_decomposition::operator()(Matrix && A0, Matrix && B0, bool with_QZ)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    clear();
    compute(A, B, with_QZ);
    return *this;
};
gschur_decomposition::gschur_decomposition(const Matrix &Q, const Matrix &Z, 
                                 const Matrix& TA, const Matrix& TB)
{
    clear();
    set(Q,Z,TA,TB, true);
};
gschur_decomposition& gschur_decomposition::set_factors(const Matrix &Q, const Matrix &Z, 
                                             const Matrix &TA, const Matrix &TB)
{
    clear();
    set(Q,Z,TA,TB, true);
    return *this;
};

gschur_decomposition& gschur_decomposition::set_factors(const Matrix &TA, const Matrix &TB)
{
    clear();
    Matrix Q,Z;
    set(Q,Z,TA,TB, false);
    return *this;
};

void gschur_decomposition::compute(const Matrix &A, const Matrix &B, bool with_QZ)
{
    if (!A.is_square() || !B.is_square() || A.rows() != B.rows())
        throw error::error_size_geig(A.rows(),A.cols(), B.rows(), B.cols());

    bool is_complex_A   = matrix_traits::is_float_complex(A.get_value_code());
    bool is_complex_B   = matrix_traits::is_float_complex(B.get_value_code());
    bool is_complex     = is_complex_A || is_complex_B;

    bool is_obj         = (A.get_value_code() == value_code::v_object)
                        ||(B.get_value_code() == value_code::v_object);

    m_has_QZ            = with_QZ;

    if (is_obj == true)
        throw error::object_value_type_not_allowed("gschur");

    Integer N           = A.rows();
    
    bool isv            = A.all_finite() && B.all_finite();
    
    value_code vt0      = matrix_traits::unify_value_types(A.get_value_code(), B.get_value_code());
    value_code vt       = matrix_traits::unify_value_types(vt0, value_code::v_float);

    if (isv == false)
    {
        m_TA            = details::make_nan_matrix(N, N, vt);
        m_TB            = m_TA;

        if (with_QZ)
        {
            m_Q_factor  = m_TA;
            m_Z_factor  = m_TA;
        };

        m_alpha     = details::make_nan_matrix(N, 1, vt);
        m_beta      = m_alpha;
        m_eig       = m_alpha;

        m_is_nan    = true;
        return;
    }
    else
    {
        m_is_nan    = false;
    };

    Integer A_ldiags    = get_ld(A,1);
    Integer A_udiags    = get_ud(A,1);

    bool A_is_tril      = A_udiags == 0;
    bool A_is_triu      = A_ldiags == 0;
    bool B_is_triu      = is_triu(B);
    bool B_is_tril      = is_tril(B);

    bool is_fac_B       = B_is_triu;
    bool is_fac_A       = A_is_triu  ||(is_complex == false && A.get_struct().is_qtriu());

    if(is_fac_B == true && is_fac_A == true)
    {
        m_TA        = convert(A, matrix_traits::get_matrix_type(vt, A.get_struct_code()));
        m_TB        = convert(B, matrix_traits::get_matrix_type(vt, B.get_struct_code()));

        if (A_is_triu)
            m_TA.add_struct(predefined_struct_type::triu);
        else 
            m_TA.add_struct(predefined_struct_type::qtriu);

        m_TB.add_struct(predefined_struct_type::triu);

        if (with_QZ)
        {
            m_Q_factor  = speye(N,N,vt);
            m_Z_factor  = m_Q_factor;
        };
        
        if (A_is_triu)
        {
            m_alpha = get_diag(m_TA);
            m_beta  = get_diag(m_TB);
            m_eig   = comp_eig(m_alpha, m_beta);
        }
        else
        {
            set_alpha_beta_eig();
        };

        return;
    }

    bool is_lfac_B  = B_is_tril;
    bool is_tril_A  = A_is_tril;
    bool is_lfac_A  = is_tril_A  ||(is_complex == false && A.get_struct().is_qtril());

    if(is_lfac_B && is_lfac_A)        
    {
        m_TA        = A(colon(end, -1, 1),colon(end, -1, 1));
        m_TB        = B(colon(end, -1, 1),colon(end, -1, 1));

        m_TA        = convert(m_TA, matrix_traits::get_matrix_type(vt, m_TA.get_struct_code()));
        m_TB        = convert(m_TB, matrix_traits::get_matrix_type(vt, m_TB.get_struct_code()));

        if (A_is_tril)
            m_TA.add_struct(predefined_struct_type::triu);
        else 
            m_TA.add_struct(predefined_struct_type::qtriu);

        m_TB.add_struct(predefined_struct_type::triu);

        if (is_tril_A)
        {
            m_alpha = get_diag(m_TA);
            m_beta  = get_diag(m_TB);
            m_eig   = comp_eig(m_alpha, m_beta);
        }
        else
        {
            set_alpha_beta_eig();
        };

        if (with_QZ)
        {
            Matrix I    = speye(N,N, vt);
            m_Q_factor  = I(colon(), colon(colon::cend(), -1, 1)); // reverse eye
            m_Z_factor  = m_Q_factor;

            struct_flag sf;
            sf.set_user(unitary_flag());
            m_Q_factor.add_struct(sf);
            m_Z_factor.add_struct(sf);
        };

        return;
    }

    // B trivial/easy
    if (B.structural_nnz() == 0 || B.get_struct().is_id())
    {
        Matrix Ac   = convert(A, matrix_traits::get_matrix_type(vt, A.get_struct_code()));

        schur_decomposition schur_obj;
        schur_obj(std::move(Ac), schur_sym_alg::dc, with_QZ);

        m_TA        = schur_obj.TA();
        m_TB        = convert(B, matrix_traits::get_matrix_type(vt, B.get_struct_code()));

        if (with_QZ)
        {
            m_Q_factor  = schur_obj.U();
            m_Z_factor  = schur_obj.U();
        };

        m_TB.add_struct(predefined_struct_type::triu);
        
        m_alpha     = schur_obj.eig();
        m_beta      = get_diag(m_TB);
        m_eig       = comp_eig(m_alpha, m_beta);

        test_factors();
        return;
    }
    else if (is_unitary(B.get_struct()))
    {
        Matrix ABprim   = mmul(A,B, trans_type::no_trans, trans_type::conj_trans);

        schur_decomposition schur_obj(ABprim, schur_sym_alg::dc, with_QZ);

        m_TA          = schur_obj.TA();
        m_TB          = speye(N, vt);
        m_alpha       = schur_obj.eig();
        m_beta        = get_diag(m_TB);
        m_eig         = schur_obj.eig();

        if (with_QZ)
        {
            m_Q_factor  = schur_obj.U();        
            m_Z_factor  = mmul(B, m_Q_factor, trans_type::conj_trans);
        };

        test_factors();
        return;
    }

    // A trivial/easy
    if (A.structural_nnz() == 0 || A.get_struct().is_id())
    {
        Matrix Bc   = convert(B,matrix_traits::get_matrix_type(vt,B.get_struct_code()));

        schur_decomposition schur_obj;            
        schur_obj(std::move(Bc), schur_sym_alg::dc, with_QZ);

        Matrix P        = schur_obj.TA();
        
        unitary_matrix UC;
        tie(UC, m_TB)   = qr2(std::move(P));
        m_TA            = mmul(A, UC, trans_type::no_trans, trans_type::conj_trans);

        if (with_QZ)
        {
            Matrix U    = schur_obj.U();
            m_Q_factor  = U * UC;
            m_Z_factor  = U;
        };

        if (is_complex_A == false)
            m_TA.add_struct(predefined_struct_type::qtriu);
        else
            m_TA.add_struct(predefined_struct_type::triu);

        set_alpha_beta_eig();

        test_factors();
        return;
    }
    else if (is_unitary(A.get_struct()))
    {
        Matrix BAprim   = mmul(B, A, trans_type::no_trans, trans_type::conj_trans);

        schur_decomposition schur_obj(std::move(BAprim), schur_sym_alg::dc, with_QZ);

        Matrix P        = schur_obj.TA();

        unitary_matrix UC;        
        tie(UC, m_TB)   = qr2(std::move(P));

        Matrix C        = UC.to_matrix();
        m_TA            = ctrans(C);

        if (with_QZ)
        {
            Matrix U    = schur_obj.U();     

            m_Q_factor  = U * C;
            m_Z_factor  = mmul(A,U,trans_type::conj_trans);
        };
        
        if (is_complex_A == false)
            m_TA.add_struct(predefined_struct_type::qtriu);
        else
            m_TA.add_struct(predefined_struct_type::triu);

        set_alpha_beta_eig();
        
        test_factors();
        return;
    }

    return details::gschur_vis::make(A,B,*this,with_QZ);
};

bool gschur_decomposition::test_factors()
{    
    bool isv_QZ     = (m_has_QZ == true) && m_Q_factor.all_finite() && m_Z_factor.all_finite()
                    || m_has_QZ == false;

    bool isv        = m_TA.all_finite() && m_TB.all_finite() && isv_QZ;

    Integer N       = m_TA.rows();

    if (isv == false)
    {
        value_code vt   = m_TA.get_value_code();

        m_TA            = details::make_nan_matrix(N, N, vt);
        m_TB            = m_TA;
        m_alpha         = details::make_nan_matrix(N, 1, vt);
        m_beta          = m_alpha;
        m_eig           = m_alpha;

        if (m_has_QZ)
        {
            m_Q_factor  = m_TA;
            m_Z_factor  = m_TA;
        };

        m_is_nan        = true;

        return false;
    }
    else
    {
        m_is_nan        = false;
    };
    
    return true;
};

void gschur_decomposition::set_alpha_beta_eig()
{
    Integer N = m_TA.rows();
    
    if (N == 0 || N == 1 || m_TA.structural_nnz() + m_TB.structural_nnz() == 0 || is_triu(m_TA))
    {
        m_alpha = get_diag(m_TA);
        m_beta  = get_diag(m_TB);
        m_eig   = comp_eig(m_alpha, m_beta);
        return;
    }

    matcl::value_code vTA   = m_TA.get_value_code();
    matcl::value_code vTB   = m_TB.get_value_code();
    matcl::value_code vt    = matrix_traits::unify_value_types(vTA, vTB);

    if (m_has_QZ)
    {
        value_code vQ       = m_Q_factor.get_value_code();
        value_code vZ       = m_Z_factor.get_value_code();
        vt                  = matrix_traits::unify_value_types(vt, vQ);
        vt                  = matrix_traits::unify_value_types(vt, vZ);
    };

    vt                      = matrix_traits::unify_value_types(vt, value_code::v_float);

    mat_code mt             = matrix_traits::get_matrix_type(vt,struct_code::struct_dense);

    m_TA                    = convert(m_TA, mt);
    m_TB                    = convert(m_TB, mt);

    if (m_has_QZ)
    {
        m_Q_factor          = convert(m_Q_factor, mt);
        m_Z_factor          = convert(m_Z_factor, mt);
    };

    switch (vt)
    {
        case value_code::v_integer:
        case value_code::v_real:
            return select_val<Real>(matcl::zeros(m_TA.rows(),1), true);
        case value_code::v_float:
            return select_val<Float>(matcl::zeros(m_TA.rows(),1), true);     
        case value_code::v_float_complex:
            return select_val<Float_complex>(matcl::zeros(m_TA.rows(),1), true);     
        case value_code::v_complex:
            return select_val<Complex>(matcl::zeros(m_TA.rows(),1), true);     
        case value_code::v_object:
        {
            throw error::object_value_type_not_allowed("gschur_decomposition");
        }
    }
}

void gschur_decomposition::set(const Matrix &Q, const Matrix &Z, const Matrix &TA, const Matrix &TB, bool with_QZ)
{
    //check structure    
    if (!TA.is_square() || !TB.is_square() || TA.rows() != TB.rows() || is_triu(TB) == false)
    {
        throw error::error_gschur_incorrect_factors();
    };

    if (with_QZ)
    {
        if (Q.rows() != Z.rows() || Q.cols() != TA.rows() || Z.cols() != TA.rows()
            || Q.rows() < Q.cols() || Z.rows() < Z.cols())
        {
            throw error::error_gschur_incorrect_factors();
        };
    };

    Integer ld  = get_ld(TA,1);

    if (ld > 1)
    {
        //S is not quasi upper triangular
        throw error::error_gschur_incorrect_factors();
    };

    bool is_complex_TA  = matrix_traits::is_float_complex(TA.get_value_code());
    bool is_complex_TB  = matrix_traits::is_float_complex(TB.get_value_code());
    bool is_complex_Q   = with_QZ? matrix_traits::is_float_complex(Q.get_value_code()) : false;
    bool is_complex_Z   = with_QZ? matrix_traits::is_float_complex(Z.get_value_code()) : false;

    bool compl  =  is_complex_TA || is_complex_TB || is_complex_Q || is_complex_Z;
    bool obj    =  (with_QZ && Q.get_value_code() == value_code::v_object)
                || (with_QZ && Z.get_value_code() == value_code::v_object)
                || (TA.get_value_code() == value_code::v_object)
                || (TB.get_value_code() == value_code::v_object);
    
    if (compl == true && ld > 0)
    {
        throw error::error_gschur_incorrect_factors();
    };

    if (obj == true)
        throw error::object_value_type_not_allowed("gschur_decomposition");

    value_code vc   = matrix_traits::unify_value_types(TA.get_value_code(),TB.get_value_code());

    if (with_QZ)
    {
        vc          = matrix_traits::unify_value_types(vc,Q.get_value_code());
        vc          = matrix_traits::unify_value_types(vc,Z.get_value_code());
    };

    m_TA            = convert_value(TA,vc);
    m_TB            = convert_value(TB,vc);
    m_Q_factor      = with_QZ? convert_value(Q,vc) : Q;
    m_Z_factor      = with_QZ? convert_value(Z,vc) : Z;
    m_has_QZ        = with_QZ;

    if (test_factors() == false)
        return;
    
    set_alpha_beta_eig();
}

Matrix gschur_decomposition::convert_value(const Matrix& mat, value_code vc) const
{
    mat_code mc     = matrix_traits::get_matrix_type(vc,mat.get_struct_code());
    return convert(mat,mc);
};

template<class Compl>
Matrix gschur_decomposition::calc_eig_inf(const Matrix& alpha, const Matrix& I) const
{
    using VR        = typename md::real_type<Compl>::type;

    Matrix a_inf    = alpha(I);
    Integer k       = a_inf.length();
    value_code vc   = matrix_traits::value_code<Compl>::value;
    Matrix eig_inf  = zeros(k, 1, vc);

    const Compl* a_ptr  = a_inf.get_array<Compl>();
    Compl* eig_ptr      = eig_inf.get_array_unique<Compl>();

    for (Integer i = 0; i < k; ++i)
    {
        Compl a     = a_ptr[i];
        VR ar       = real(a);
        VR ai       = imag(a);

        Compl e;

        if (ar == 0. && ai == 0.)
            e       = constants::nan<VR>();
        else if (ar == 0. && ai != 0.)
            e       = Compl(VR(0.), ai * constants::inf<VR>());
        else if (ar != 0. && ai == 0.)
            e       = Compl(ar * constants::inf<VR>(), VR(0.));
        else
            e       = Compl(ar * constants::inf<VR>(), ai * constants::inf<VR>());

        eig_ptr[i]  = e;
    };

    return eig_inf;
};

Matrix gschur_decomposition::comp_eig(const Matrix& alpha, const Matrix& beta) const
{
    Matrix I            = find(beta == 0.);
    bool is_alpha_compl = matrix_traits::is_float_complex(alpha.get_value_code()) == true;

    if (I.is_empty() == true || is_alpha_compl == false)
    {
        Matrix ret  = div(alpha,beta);
        return ret;
    };
    
    //alpha is complex
    Matrix eig_inf;
    if (alpha.get_value_code() == value_code::v_complex)
        eig_inf     = calc_eig_inf<Complex>(alpha,I);
    else
        eig_inf     = calc_eig_inf<Float_complex>(alpha,I);

    Matrix ret      = div(alpha,beta);
    ret(I)          = eig_inf;

    return ret;
};

void gschur_decomposition::select(const Matrix& ind)
{
    Integer M;
    details::schur_reorder_check(ind, m_TA.rows(), M);

    // FAST EXIT
    Integer N = m_TA.rows();
    if (N == 0 || N == 1 || m_TA.structural_nnz() + m_TA.structural_nnz() == 0
            || (m_TA.get_struct().is_id() && m_TB.get_struct().is_id()) )
    {
        return;
    }

    if (m_is_nan == true)
        return;

    value_code vc   = matrix_traits::unify_value_types(m_TA.get_value_code(),m_TB.get_value_code());

    if (m_has_QZ)
    {
        vc          = matrix_traits::unify_value_types(vc,m_Q_factor.get_value_code());
        vc          = matrix_traits::unify_value_types(vc,m_Z_factor.get_value_code());
    };

    switch (vc)
    {
        case value_code::v_integer:
        case value_code::v_real:
            return select_val<Real>(ind, false);
        case value_code::v_float:
            return select_val<Float>(ind, false);
        case value_code::v_float_complex:
            return select_val<Float_complex>(ind, false);
        case value_code::v_complex:
            return select_val<Complex>(ind, false);
        case value_code::v_object:
        {
            throw error::object_value_type_not_allowed("gschur_decomposition");
        }
    }
};

template<class Val>
void gschur_decomposition::select_val(const Matrix& ind, const bool no_fast_exit)
{
    using constants::nan;
    
    if (test_factors() == false)
        return;

    Integer N   = m_TA.rows();

    Integer M_sel;    
    details::schur_reorder_check(ind, N, M_sel);
    
    // fast exit for trivial ind vectors like [0 0 0 0] [1 1 1] [1 1 0 0 0] etc.
    if (!no_fast_exit && details::schur_is_trivial_reorder(ind))
        return;    

    using Mat = raw::Matrix<Val,struct_dense>;

    m_TA            = convert(m_TA, Mat::matrix_code);
    m_TB            = convert(m_TB, Mat::matrix_code);

    Mat mat_TA      = m_TA.get_impl_unique<Mat>();
    Mat mat_TB      = m_TB.get_impl_unique<Mat>();

    mat_TA.set_struct(struct_flag());
    mat_TB.set_struct(struct_flag());

    Mat mat_Q       = Mat(ti::ti_type<Val>());
    Mat mat_Z       = Mat(ti::ti_type<Val>());

    if (m_has_QZ)
    {
        m_Q_factor  = convert(m_Q_factor, Mat::matrix_code);
        m_Z_factor  = convert(m_Z_factor, Mat::matrix_code);

        mat_Q.assign_to_fresh(m_Q_factor.get_impl_unique<Mat>());
        mat_Z.assign_to_fresh(m_Z_factor.get_impl_unique<Mat>());

        mat_Q.set_struct(struct_flag());
        mat_Z.set_struct(struct_flag());
    };

    Val *r_TA       = mat_TA.ptr();
    Val *r_TB       = mat_TB.ptr();
    Integer TA_ld   = mat_TA.ld();
    Integer TB_ld   = mat_TB.ld();

    Val *r_Q        = mat_Q.ptr();
    Val *r_Z        = mat_Z.ptr();
    Integer Q_ld    = mat_Q.ld();
    Integer Z_ld    = mat_Z.ld();
    Integer Q_r     = mat_Q.rows();

    using matcl::details::lap;

    //transform quasi-triangular matrix to schur form
    matcl::lapack::qtriu2gschur(m_TA.rows(), Q_r, lap(r_TA), TA_ld, lap(r_TB), TB_ld,  m_has_QZ, m_has_QZ, 
                    lap(r_Q), Q_ld, lap(r_Z), Z_ld);

    using VC    = typename md::complex_type<Val>::type;
    using VR    = typename md::real_type<Val>::type;
    using Mat_C = raw::Matrix<VC,struct_dense>;

    Mat_C   alpha(ti::ti_type<Val>(),N,1);
    Mat     beta(ti::ti_type<Val>(),N,1);

    Integer M;
    Integer info    = 0;
    Integer wantq   = m_has_QZ;
    Integer wantz   = m_has_QZ;

    Val work_query;
    Integer iwork_query;
    
    //TODO: avoid copying
    raw::integer_dense I_r = convert(ind, mat_code::integer_dense).get_impl_unique<raw::integer_dense>();

    //TODO
    //lapack::tgsen3(wantq, wantz, lap(I_r.ptr()), N, Q_r, lap(r_TA), TA_ld, lap(r_TB), TB_ld,
    //                lap(alpha.ptr()), lap(beta.ptr()), lap(r_Q), Q_ld, lap(r_Z), Z_ld, &M, 
    //                lap(&work_query), -1, &info);

    lapack::tgsen(0, wantq, wantz, lap(I_r.ptr()), N, lap(r_TA), TA_ld, lap(r_TB), TB_ld,
                    lap(alpha.ptr()), lap(beta.ptr()), lap(r_Q), Q_ld, lap(r_Z), Z_ld, &M,
                    nullptr, nullptr, nullptr, lap(&work_query), -1, lap(&iwork_query), -1, &info);    

    Integer lwork   = (Integer) real(work_query);
    lwork           = lwork + 1; // + 1 to work around bug in MKL lapack (old lapack 3.1.0)

    using VTR               = pod_type<Val>;
    using workspace         = matcl::pod_workspace<VTR>;
    workspace WORK          = workspace(lwork);
    Val* ptr_WORK           = reinterpret_cast<Val*>(WORK.ptr());

    //TODO
    //lapack::tgsen3(wantq, wantz, lap(I_r.ptr()), N, Q_r, lap(r_TA), TA_ld, lap(r_TB), TB_ld,
    //                lap(alpha.ptr()), lap(beta.ptr()), lap(r_Q), Q_ld, lap(r_Z), Z_ld, &M, 
    //                lap(ptr_WORK), lwork, &info);

    Integer liwork  = iwork_query + 1;

    using iworkspace        = matcl::pod_workspace<Integer>;
    iworkspace IWORK        = iworkspace(liwork);
    Integer* ptr_IWORK      = IWORK.ptr();

    lapack::tgsen(0, wantq, wantz, lap(I_r.ptr()), N, lap(r_TA), TA_ld, lap(r_TB), TB_ld,
                    lap(alpha.ptr()), lap(beta.ptr()), lap(r_Q), Q_ld, lap(r_Z), Z_ld, &M, 
                    nullptr, nullptr, nullptr, lap(ptr_WORK), lwork, lap(ptr_IWORK), liwork, &info);

    if (info)
    {
        clear();
        throw error::error_gschur_sel_comp();
    };

    m_alpha = Matrix(alpha,false);
    m_beta  = Matrix(beta,false);
    m_eig   = comp_eig(m_alpha,m_beta);
    
    m_TA    = Matrix(mat_TA, true);
    m_TB    = Matrix(mat_TB, true);

    if (md::is_complex<Val>::value == false)
    {
        Integer ld = get_ld(m_TA, 1);

        if (ld == 1)
            m_TA.add_struct(predefined_struct_type::qtriu);
        else
            m_TA.add_struct(predefined_struct_type::triu);        
    }
    else
    {
        m_TA.add_struct(predefined_struct_type::triu);  
    };

    m_TB.add_struct(predefined_struct_type::triu);  

    if (m_has_QZ)
    {
        struct_flag sf_u;
        sf_u.set_user(unitary_flag());

        m_Q_factor  = Matrix(mat_Q, true);
        m_Z_factor  = Matrix(mat_Z, true);

        m_Q_factor.set_struct(sf_u);
        m_Z_factor.set_struct(sf_u);    
    };

    return;
};

Matrix gschur_decomposition::left_eigenvectors() const
{
    Integer N = m_TA.rows();
    return left_eigenvectors(iones(N,1));
};
Matrix gschur_decomposition::right_eigenvectors() const
{
    Integer N = m_TA.rows();
    return right_eigenvectors(iones(N,1));
};
tuple<Matrix,Matrix> gschur_decomposition::eigenvectors() const
{
    Integer N = m_TA.rows();
    return eigenvectors(iones(N,1));
};

Matrix gschur_decomposition::left_eigenvectors(const Matrix& ind0) const
{
    //increase refcount
    Matrix ind(ind0);

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, true, false);
    return XL;
};
Matrix gschur_decomposition::left_eigenvectors(Matrix&& ind0) const
{
    //increase refcount
    Matrix ind(std::move(ind0));

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, true, false);
    return XL;
};
Matrix gschur_decomposition::right_eigenvectors(const Matrix& ind0) const
{
    //increase refcount
    Matrix ind(ind0);

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, false, true);
    return XR;
};

Matrix gschur_decomposition::right_eigenvectors(Matrix&& ind0) const
{
    //increase refcount
    Matrix ind(std::move(ind0));

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, false, true);
    return XR;
};
tuple<Matrix,Matrix> gschur_decomposition::eigenvectors(const Matrix& ind0) const
{
    //increase refcount
    Matrix ind(ind0);

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, true, true);
    return tuple<Matrix,Matrix>(XL,XR);
};

tuple<Matrix,Matrix> gschur_decomposition::eigenvectors(Matrix&& ind0) const
{
    //increase refcount
    Matrix ind(std::move(ind0));

    Matrix XL, XR;
    comp_eigenvectors(ind, XL, XR, true, true);
    return tuple<Matrix,Matrix>(XL,XR);
};

void gschur_decomposition::comp_eigenvectors(Matrix& ind, Matrix& XL, Matrix& XR, 
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

    return details::gschur_visitor_eigvec
        ::make<const Matrix&>(m_TA,ind,*this, XL, XR, comp_left, comp_right);
};

Matrix gschur_decomposition::rcond_eig(const Matrix& VL, const Matrix& VR) const
{
    Integer N = m_TA.rows();
    return rcond_eig(VL, VR, iones(N,1));
}
Matrix gschur_decomposition::rcond_eig() const
{
    Matrix VL, VR;
    tie(VL, VR) = eigenvectors();
    return rcond_eig(VL, VR);
};

Matrix gschur_decomposition::rcond_eig(const Matrix& ind) const
{
    Matrix VL, VR;
    tie(VL, VR) = eigenvectors(ind);
    return rcond_eig(VL, VR, ind);
};

Matrix gschur_decomposition::rcond_vec() const
{
    Integer N = m_TA.rows();
    return rcond_vec(iones(N,1));
};
Matrix gschur_decomposition::rcond_vec(const Matrix& ind) const
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

    Matrix ret;
    details::gschur_visitor_cond_vec::make<const Matrix&>(m_TA, ret, *this, ind);
    return ret;
};

Matrix gschur_decomposition::rcond_eig(const Matrix& VL, const Matrix& VR, const Matrix& ind) const
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

    Matrix ret;
    details::gschur_visitor_cond_eig::make<const Matrix&>(m_TA, ret, *this, VL, VR, ind);
    return ret;
};

tuple<Real,Real> gschur_decomposition::rcond_projector(Integer M) const
{
    Real PL, PR;
    rcond_est(M, 1, PL, PR);
    return tuple<Real,Real>(PL,PR);
};
tuple<Real,Real> gschur_decomposition::separation_est(Integer M) const
{
    Real Difu, Difl;
    rcond_est(M, 2, Difu, Difl);
    return tuple<Real,Real>(Difu,Difl);
};
tuple<Real,Real> gschur_decomposition::separation(Integer M) const
{
    Real Difu, Difl;
    rcond_est(M, 3, Difu, Difl);
    return tuple<Real,Real>(Difu,Difl);
};

Real gschur_decomposition::pert_bound(Integer M, bool use_est) const
{
    Real PL, PR, Difu, Difl;
    tie(PL, PR) = rcond_projector(M);

    if (use_est)
        tie(Difu,Difl) = separation_est(M);
    else
        tie(Difu,Difl) = separation(M);

    Real p = std::sqrt(1.0/(PL*PL)+1.0/(PR*PR)) + 2.0 * std::max(1.0/PL,1.0/PR);

    if (p == 0.0)
        return 0.0;

    Real x = std::min(Difu,Difl) / p;
    return x;
};

void gschur_decomposition::rcond_est(Integer M, Integer type, Real& ret1, Real& ret2) const
{
    Integer N = m_TA.rows();

    if (M < 0 || M > N)
        throw error::error_schur_eig_cluster_size(M, N);

    if (m_is_nan)
    {
        ret1 = constants::nan();
        ret2 = constants::nan();
        return;
    }

    return details::gschur_visitor_rcond::make<const Matrix&>(m_TA, *this, M, type, ret1, ret2);
};

};

#pragma warning( pop )