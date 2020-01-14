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

#include "matcl-linalg/decompositions/hess.h"
#include "matcl-linalg/decompositions/qr.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/general/control_exceptions.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/inplace.h"
#include "matcl-linalg/decompositions/householder_q.h"
#include "matcl-linalg/decompositions/qr.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-internals/base/sort.h"
#include "matcl-internals/base/utils.h"
#include "matcl-linalg/linear_eq/linsolve_objects_decomp.h"
#include "matcl-linalg/graph/matcl_graph.h"
#include "matcl-linalg/decompositions/eig/schur_utils.h"
#include "matcl-linalg/utils/optim_params.h"

namespace matcl { namespace details
{

template<class V> struct hess_un_or_he_sy_selector 
{
    using Mat = raw::Matrix<V, struct_dense>;

    static void eval_trd(char uplo, Integer N, Mat& Qc, Mat& tau, Mat& DE, Mat& work, 
                         Integer lwork, Integer* info2)
    {
        V* ptr_DE       = DE.ptr();
        Integer DE_LD   = DE.ld();
        V* ptr_D        = ptr_DE;
        V* ptr_E        = ptr_DE + DE_LD;

        lapack::sytrd(&uplo, N, lap(Qc.ptr()), Qc.ld(), lap(ptr_D), lap(ptr_E), 
                    lap(tau.ptr()), lap(work.ptr()), lwork, lap(info2));
    };

    static bool is_sym_her(const Mat& A)
    {
        return A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));
    }
};

template<>
struct hess_un_or_he_sy_selector<Complex>
{
    using VR        = Real;
    using Mat       = raw::Matrix<Complex, struct_dense>;
    using de_Mat    = raw::Matrix<Real, struct_dense>;

    static void eval_trd(char uplo, Integer N, Mat& Qc, Mat& tau, de_Mat& DE, Mat& work, 
                         Integer lwork, Integer* info2)
    {
        VR* ptr_DE      = DE.ptr();
        Integer DE_LD   = DE.ld();
        VR* ptr_D       = ptr_DE;
        VR* ptr_E       = ptr_DE + DE_LD;

        lapack::zhetrd(&uplo, N, lap(Qc.ptr()), Qc.ld(), lap(ptr_D), lap(ptr_E), 
                        lap(tau.ptr()), lap(work.ptr()), lwork, lap(info2));
    };

    static bool is_sym_her(const Mat& A)
    {
        return A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));
    }
};

template<>
struct hess_un_or_he_sy_selector<Float_complex>
{
    using VR        = Float;
    using Mat       = raw::Matrix<Float_complex, struct_dense>;
    using de_Mat    = raw::Matrix<Float, struct_dense>;

    static void eval_trd(char uplo, Integer N, Mat& Qc, Mat& tau, de_Mat& DE, Mat& work, 
                         Integer lwork, Integer* info2)
    {
        VR* ptr_DE      = DE.ptr();
        Integer DE_LD   = DE.ld();
        VR* ptr_D       = ptr_DE;
        VR* ptr_E       = ptr_DE + DE_LD;

        lapack::chetrd(&uplo, N, lap(Qc.ptr()), Qc.ld(), lap(ptr_D), lap(ptr_E), 
                    lap(tau.ptr()), lap(work.ptr()), lwork, lap(info2));
    };

    static bool is_sym_her(const Mat& A)
    {
        return A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));
    }
};

struct hess_eig_data
{
    const Matrix*   m_E;
    Matrix&         m_VL;
    Matrix&         m_VR;
    const Matrix*   m_init_L;
    const Matrix*   m_init_R;
    bool            m_eval_L;
    bool            m_eval_R;

    hess_eig_data(const Matrix& E, Matrix& VL, Matrix& VR, const Matrix& init_L, const Matrix& init_R, 
                  bool eval_L, bool eval_R)
        :m_E(&E), m_VL(VL), m_VR(VR), m_init_L(&init_L), m_init_R(&init_R), m_eval_L(eval_L), m_eval_R(eval_R)
    {};

    hess_eig_data(const hess_eig_data&) = delete;
};

//--------------------------------------------------------------------------------
//                                  HESS_EIG
//--------------------------------------------------------------------------------
template<class V, bool Mat_real>
struct init_WI
{};
template<class V>
struct init_WI<V,true>
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<VR,struct_dense>;

    static void eval(VR* ptr_WR, VR* ptr_WI, Integer M, const Matrix& E, bool& has_complex_eig, 
                     Integer& error)
    {
        const VR* ptr_E = E.get_array<VR>();
        error           = 0;
        has_complex_eig = false;

        for (Integer i = 0; i < M; ++i)
            ptr_WR[i]   = ptr_E[i];

        for (Integer i = 0; i < M; ++i)
            ptr_WI[i]   = VR(0.0);
    };
};

template<class V>
struct init_WI<V,false>
{
    using VR    = typename md::real_type<V>::type;
    using VC    = typename md::complex_type<V>::type;
    using Mat   = raw::Matrix<VC,struct_dense>;

    static void eval(VR* ptr_WR, VR* ptr_WI, Integer M, const Matrix& E, bool& has_complex_eig, Integer& error)
    {
        const VC* ptr_E = E.get_array<VC>();
        error           = 0;
        has_complex_eig = false;

        for (Integer i = 0; i < M; ++i)
        {
            VR re           = real(ptr_E[i]);
            VR im           = imag(ptr_E[i]);

            ptr_WR[i]       = re;
            ptr_WI[i]       = im;

            if (im != VR(0.0))
            {
                has_complex_eig = true;

                if (i == M-1)
                {
                    error   = i+1 + 1;
                    return;
                };

                VR re2      = real(ptr_E[i+1]);
                VR im2      = imag(ptr_E[i+1]);
                im          = -im;

                if (re2 != re || im2 != im)
                {
                    error   = i+1+ 1;
                    return;
                };

                ptr_WR[i+1] = re;
                ptr_WI[i+1] = im;

                ++i;
            };
        };        
    };
};

template<class V, bool Mat_real>
struct init_WI_compl
{};

template<class V>
struct init_WI_compl<V,true>
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<VR,struct_dense>;

    static const bool Mat_real  = true;

    static void eval(V* ptr_W, Integer M, const Matrix& E)
    {
        const VR* ptr_E = E.get_array<VR>();

        for (Integer i = 0; i < M; ++i)
            ptr_W[i]    = V(ptr_E[i]);
    };
};

template<class V>
struct init_WI_compl<V,false>
{
    using Mat   = raw::Matrix<V,struct_dense>;

    static void eval(V* ptr_W, Integer M, const Matrix& E)
    {
        const V* ptr_E = E.get_array<V>();

        for (Integer i = 0; i < M; ++i)
            ptr_W[i]    = ptr_E[i];
    };
};

template<class V, bool Mat_real>
struct init_W
{};

template<class V>
struct init_W<V,true>
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<VR,struct_dense>;

    static const bool Mat_real  = true;

    static void eval(VR* ptr_W, Integer M, const Matrix& E, Integer& error)
    {
        const VR* ptr_E = E.get_array<VR>();
        error = 0;

        for (Integer i = 0; i < M; ++i)
            ptr_W[i]    = ptr_E[i];
    };
};

template<class V>
struct init_W<V,false>
{
    using VR    = typename md::real_type<V>::type;
    using VC    = typename md::complex_type<V>::type;
    using Mat   = raw::Matrix<VR,struct_dense>;

    static const bool Mat_real  = true;

    static void eval(VR* ptr_W, Integer M, const Matrix& E, Integer& error)
    {
        const VC* ptr_E = E.get_array<VC>();

        error = 0;

        for (Integer i = 0; i < M; ++i)
        {
            VR im       = imag(ptr_E[i]);

            if (error == 0 && im != VR(0.0))
                error = i+1;

            ptr_W[i]    = real(ptr_E[i]);
        };
    };
};

template<class V, bool Is_mat_real>
struct init_VLR_storage
{};

template<class V>
struct init_VLR_storage<V,true>
{
    using VR    = V;
    using Mat_R = raw::Matrix<VR,struct_dense>;
    using Mat   = raw::Matrix<V,struct_dense>;

    static const bool Mat_real  = true;

    static void eval(Mat& VLR, const Matrix& init_mat, const VR* ptr_WI)
    {
        const Mat_R& rep    = convert(init_mat, Mat_R::matrix_code).get_impl<Mat_R>();
        const VR* ptr_rep   = rep.ptr();
        Integer rep_LD      = rep.ld();

        V* ptr_VLR          = VLR.ptr();
        Integer VLR_ld      = VLR.ld();

        Integer N           = VLR.rows();
        Integer M           = rep.cols();

        for (Integer i = 0; i < M; ++i)
        {
            bool is_compl   = ptr_WI[i] != VR();

            for (Integer j = 0; j < N; ++j)
                ptr_VLR[j]  = ptr_rep[j];

            ptr_VLR         += VLR_ld;
            ptr_rep         += rep_LD;

            if (is_compl == true)
            {
                ++i;

                if (i == M)
                {
                    //not enough place for storing imaginary part,
                    //this case should be impossible
                    return;
                };
                
                for (Integer j = 0; j < N; ++j)
                    ptr_VLR[j]  = VR(0.0);

                ptr_VLR     += VLR_ld;
                ptr_rep     += rep_LD;                
            };
        };
    };
};

template<class V>
struct init_VLR_storage<V,false>
{
    using VR    = V;
    using VC    = typename md::complex_type<V>::type;
    using Mat_C = raw::Matrix<VC,struct_dense>;
    using Mat   = raw::Matrix<V,struct_dense>;

    static const bool Mat_real  = false;

    static void eval(Mat& VLR, const Matrix& init_mat, const VR* ptr_WI)
    {
        const Mat_C& rep    = convert(init_mat, Mat_C::matrix_code).get_impl<Mat_C>();
        const VC* ptr_rep   = rep.ptr();
        Integer rep_LD      = rep.ld();

        V* ptr_VLR          = VLR.ptr();
        Integer VLR_ld      = VLR.ld();

        Integer N           = VLR.rows();
        Integer M           = rep.cols();

        for (Integer i = 0; i < M; ++i)
        {
            bool is_compl   = ptr_WI[i] != VR();

            for (Integer j = 0; j < N; ++j)
                ptr_VLR[j]  = real(ptr_rep[j]);

            ptr_VLR         += VLR_ld;      

            if (is_compl == true)
            {
                ++i;

                if (i == M)
                {
                    //not enough place for storing imaginary part,
                    //this case should be impossible
                    return;
                };

                for (Integer j = 0; j < N; ++j)
                    ptr_VLR[j]  = imag(ptr_rep[j]);

                ptr_rep     += rep_LD;
                ptr_VLR     += VLR_ld;
            };

            ptr_rep         += rep_LD;
        };
    };
};

template<class V>
struct hess_eig_tridiag
{
    using VR    = typename md::real_type<V>::type;
    using VC    = typename md::complex_type<V>::type;

    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_I = raw::Matrix<Integer,struct_dense>;
    using Mat_C = raw::Matrix<VC,struct_dense>;    
    using Mat_R = raw::Matrix<VR,struct_dense>;    

    static void eval(const Matrix& D0, const Matrix& D1, hess_eig_data& data)
    {
        Integer N           = D0.length();
        Integer M           = data.m_E->length();
        bool eval_L         = data.m_eval_L;
        bool eval_R         = data.m_eval_R;
        bool is_D1_compl    = matrix_traits::is_float_complex(D1.get_value_code()) == true;

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        using rworkspace    = matcl::pod_workspace<VR>;

        iworkspace work_BL  = iworkspace(M);
        rworkspace WORK     = rworkspace(N+5*N);
        iworkspace work_Wp  = iworkspace(2*M);        
        iworkspace IWORK    = iworkspace(N);
        workspace WORK_U    = workspace(is_D1_compl? N : 0);

        Integer* ptr_BL     = work_BL.ptr();
        VR* ptr_W           = WORK.ptr();
        Integer* ptr_Wp     = work_Wp.ptr();
        Integer* ptr_Wp_w   = ptr_Wp + M;
        VR* ptr_WORK        = WORK.ptr() + N;
        Integer* ptr_IWORK  = IWORK.ptr();
        V* ptr_U            = reinterpret_cast<V*>(WORK_U.ptr());

        // diagonal elements of hermitian matrices are real
        const VR* ptr_D     = D0.get_array<VR>();

        Matrix D1r          = make_tridiag_subdiag_real(D1, ptr_U, is_D1_compl);
        const VR* ptr_E     = D1r.get_array<VR>();

        //init block info; we assume that there exists only one block
        for (Integer i = 0; i < M; ++i)
            ptr_BL[i]       = 1;

        bool E_is_real      = matrix_traits::is_float_real(data.m_E->get_value_code()) == true;
        Integer error;

        if (E_is_real)
        {
            static const bool Mat_real  = true;
            init_W<V,Mat_real>::eval(ptr_W, M, *data.m_E, error);
        }
        else
        {
            static const bool Mat_real  = false;
            init_W<V,Mat_real>::eval(ptr_W, M, *data.m_E, error);
        };

        if (error != 0)
            throw error::error_hess_eig_her_complex_eig(error);

        //sort W
        for (Integer i = 0; i < M; ++i)
            ptr_Wp[i] = i+1;

        matcl::utils::sort_q(ptr_W, ptr_Wp, M);
        lapack::perm2int(M, ptr_Wp, ptr_Wp_w, false);

        Mat_I ret_fail(ti::ti_type<Integer>(), M, 1);
        Mat   ret_VLR(ti::ti_type<V>(), N, M);

        Integer* ptr_fail   = ret_fail.ptr();
        V* ptr_VLR          = ret_VLR.ptr();
        Integer ld_VLD      = ret_VLR.ld();

        Integer info;

        //we assume, that there exists only one block, ISPLIT should be an array of size
        //N, but only the first element should be accessed
        Integer ISPLIT      = N;
        lapack::stein(N, lap(ptr_D), lap(ptr_E), M, lap(ptr_W), ptr_BL, &ISPLIT, lap(ptr_VLR), ld_VLD, 
                      lap(ptr_WORK), ptr_IWORK, ptr_fail, &info );

        if (info < 0)
            throw error::error_general("invalid argument passed to hsein");

        //reorder eigenvectors
        lapack::laswpc(N, lap(ptr_VLR), ld_VLD, 1, M, ptr_Wp, 1);

        // rescale eigenvalues
        if (is_D1_compl == true)
        {
            V* ptr_VLR2 = ptr_VLR;

            for (Integer i = 0; i < M; ++i)
            {
                for (Integer j = 0; j < N-1; ++j)
                    ptr_VLR2[j] = ptr_VLR2[j] * ptr_U[j];

                ptr_VLR2    += ld_VLD;
            };
        };

        data.m_VL = Matrix(ret_VLR,true);
        data.m_VR = Matrix(ret_VLR,true);

        if (info > 0)
        {
            Matrix fail     = Matrix(ret_fail,false);
            throw error::hess_eig_failed(data.m_VL, data.m_VR, fail, fail, eval_L, eval_R, info);
        };

        return;
    };
};

template<class V, bool Is_real = md::is_float_real_scalar<V>::value>
struct hess_eig_dense
{
    using VR    = V;
    using VC    = typename md::complex_type<V>::type;

    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_I = raw::Matrix<Integer,struct_dense>;
    using Mat_C = raw::Matrix<VC,struct_dense>;    

    static void eval(const Mat& A, hess_eig_data& data)
    {
        bool eval_L         = data.m_eval_L;
        bool eval_R         = data.m_eval_R;
        Integer N           = A.rows();
        Integer M           = data.m_E->length();
        bool has_init_L     = data.m_init_L->is_empty() == false;
        bool has_init_R     = data.m_init_R->is_empty() == false;
        bool has_init       = has_init_L || has_init_R;

        const char* SIDE    = (eval_L && eval_R)? "B" : (eval_L? "L" : "R");
        const char* EIGSRC  = "N";
        const char* INITV   = has_init? "U" : "N";

        const V* ptr_H      = A.ptr();
        Integer LDH         = A.ld();

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        using rworkspace    = matcl::pod_workspace<VR>;

        iworkspace work_sel = iworkspace(N);
        rworkspace work_WR  = rworkspace(N);
        rworkspace work_WI  = rworkspace(N);
        workspace WORK      = workspace((N+2)*N);

        Integer* ptr_sel    = work_sel.ptr();
        VR* ptr_WR          = work_WR.ptr();
        VR* ptr_WI          = work_WI.ptr();
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());
        Integer VLR_LD      = N;

        init_sel(ptr_sel, N, M);

        bool E_is_real      = matrix_traits::is_float_real(data.m_E->get_value_code()) == true;
        Integer error;
        bool has_complex_eig;

        if (E_is_real)
        {
            static const bool Mat_real  = true;
            init_WI<V,Mat_real>::eval(ptr_WR, ptr_WI, M, *data.m_E, has_complex_eig, error);
        }
        else
        {
            static const bool Mat_real  = false;
            init_WI<V,Mat_real>::eval(ptr_WR, ptr_WI, M, *data.m_E, has_complex_eig, error);
        };

        if (error != 0)
            throw error::error_hess_eig_invalid_conj_pair(error, error > M);

        Integer M_L         = eval_L? M : 0;
        Integer M_R         = eval_R? M : 0;

        Mat_I ret_fail_L(ti::ti_type<Integer>(), M_L, 1);
        Mat_I ret_fail_R(ti::ti_type<Integer>(), M_R, 1);

        Integer* ptr_fail_L = ret_fail_L.ptr();
        Integer* ptr_fail_R = ret_fail_R.ptr();

        //create matrices, make real representation of complex vectors
        Mat ret_VL          = init_VLR(has_init_L, N, M, *data.m_init_L, ptr_WI);
        Mat ret_VR          = init_VLR(has_init_R, N, M, *data.m_init_R, ptr_WI);

        ret_VL.get_struct().reset();
        ret_VR.get_struct().reset();

        V* ptr_VL           = reinterpret_cast<V*>(ret_VL.ptr());
        V* ptr_VR           = reinterpret_cast<V*>(ret_VR.ptr());

        Integer info;

        lapack::hsein( SIDE, EIGSRC, INITV, ptr_sel, N, lap(ptr_H), LDH, lap(ptr_WR), lap(ptr_WI),
                    lap(ptr_VL), VLR_LD, lap(ptr_VR), VLR_LD, M, M, lap(ptr_WORK), ptr_fail_L, ptr_fail_R, &info );

        if (info < 0)
            throw error::error_general("invalid argument passed to hsein");

        if (has_complex_eig == false)
        {
            data.m_VL = Matrix(ret_VL,true);
            data.m_VR = Matrix(ret_VR,true);
        }
        else
        {
            Mat_C ret_VL2   = make_complex_VLR(ret_VL, ptr_WI, eval_L, M);
            Mat_C ret_VR2   = make_complex_VLR(ret_VR, ptr_WI, eval_R, M);

            data.m_VL = Matrix(ret_VL2,true);
            data.m_VR = Matrix(ret_VR2,true);
        };

        if (info > 0)
        {
            // created matrices ret_fail_L, ret_fail_R are OK
            Matrix fail_L   = Matrix(ret_fail_L,false);
            Matrix fail_R   = Matrix(ret_fail_R,false);

            throw error::hess_eig_failed(data.m_VL, data.m_VR, fail_L, fail_R, eval_L, eval_R, info);
        };

        return;
    };

    static void init_sel(Integer* ptr_sel, Integer N, Integer M)
    {
        for (Integer i = 0; i < M; ++i)
            ptr_sel[i] = 1;

        for (Integer i = M; i < N; ++i)
            ptr_sel[i] = 0;
    };

    static Mat init_VLR(bool has_init, Integer N, Integer M, const Matrix& init_mat,
                        const VR* ptr_WI)
    {
        if (has_init == false)
            return Mat(ti::ti_type<V>(), N, M);

        value_code this_vc  = matrix_traits::value_code<V>::value;
        mat_code this_mc    = matrix_traits::get_matrix_type(this_vc, struct_code::struct_dense);
        mat_code init_mc    = init_mat.get_matrix_code();

        if (init_mc == this_mc)
            return init_mat.get_impl<Mat>().make_unique();

        Mat VLR(ti::ti_type<V>(), N, M);

        bool init_is_real   = matrix_traits::is_float_real(init_mat.get_value_code()) == true;

        if (init_is_real)
            init_VLR_storage<V,true>::eval(VLR, init_mat, ptr_WI);
        else
            init_VLR_storage<V,false>::eval(VLR, init_mat, ptr_WI);   

        return VLR;
    };

    static Mat_C make_complex_VLR(const Mat& ret_VL, const VR* ptr_WI, bool eval_L, Integer M)
    {
        if (!eval_L)
            return Mat_C(ret_VL.get_type());

        Integer N           = ret_VL.rows();

        Mat_C ret           = Mat_C(ret_VL.get_type(), N, M);

        Integer ret_ld      = ret.ld();
        Integer VL_ld       = ret_VL.ld();

        VC* ptr_ret         = ret.ptr();
        const V* ptr_VL     = ret_VL.ptr();
        const V* ptr_VL2    = ret_VL.ptr() + VL_ld;

        for (Integer i = 0; i < M; ++i)
        {
            bool im         = ptr_WI[i] != V(0.0);

            if (im == false)
            {
                for (Integer j = 0; j < N; ++j)
                    ptr_ret[j]  = VC(ptr_VL[j], VR(0.0));

                ptr_ret     += ret_ld;
                ptr_VL      += VL_ld;
                ptr_VL2     += VL_ld;
            }
            else
            {
                ++i;

                for (Integer j = 0; j < N; ++j)
                    ptr_ret[j]  = VC(ptr_VL[j], ptr_VL2[j]);

                ptr_ret         += ret_ld;                

                if (i == M)
                {
                    //this should be impossible
                    return ret;
                };

                for (Integer j = 0; j < N; ++j)
                    ptr_ret[j]  = VC(ptr_VL[j], -ptr_VL2[j]);

                ptr_ret += ret_ld;
                ptr_VL  += VL_ld * 2;
                ptr_VL2 += VL_ld * 2;
            };
        };

        return ret;
    };
};

template<class V>
struct hess_eig_dense<V,false>
{
    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_I = raw::Matrix<Integer,struct_dense>;

    using VR                = typename md::real_type<V>::type;

    static void eval(const Mat& A, hess_eig_data& data)
    {
        bool eval_L         = data.m_eval_L;
        bool eval_R         = data.m_eval_R;
        Integer N           = A.rows();
        Integer M           = data.m_E->length();
        bool has_init_L     = data.m_init_L->is_empty() == false;
        bool has_init_R     = data.m_init_R->is_empty() == false;
        bool has_init       = has_init_L || has_init_R;

        Integer M_L         = eval_L? M : 0;
        Integer M_R         = eval_R? M : 0;

        const char* SIDE    = (eval_L && eval_R)? "B" : (eval_L? "L" : "R");
        const char* EIGSRC  = "N";
        const char* INITV   = has_init? "U" : "N";

        const V* ptr_H      = A.ptr();
        Integer LDH         = A.ld();

        Mat_I ret_fail_L(ti::ti_type<Integer>(), M_L, 1);
        Mat_I ret_fail_R(ti::ti_type<Integer>(), M_R, 1);

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        using rworkspace    = matcl::pod_workspace<VR>;

        iworkspace work_sel = iworkspace(N);
        workspace work_W    = workspace(N);
        workspace WORK      = workspace(N*N);
        rworkspace RWORK    = rworkspace(N);

        Integer* ptr_sel    = work_sel.ptr();
        V* ptr_W            = reinterpret_cast<V*>(work_W.ptr());
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());
        VR* ptr_RWORK       = RWORK.ptr();
        Integer* ptr_fail_L = ret_fail_L.ptr();
        Integer* ptr_fail_R = ret_fail_R.ptr();
        Integer VLR_LD      = N;

        init_sel(ptr_sel, N, M);

        bool E_is_real              = matrix_traits::is_float_real(data.m_E->get_value_code()) == true;
        if (E_is_real)
        {
            static const bool Mat_real  = true;
            init_WI_compl<V, Mat_real>::eval(ptr_W, M, *data.m_E);
        }
        else
        {
            static const bool Mat_real  = false;
            init_WI_compl<V,Mat_real>::eval(ptr_W, M, *data.m_E);
        };

        //create matrices, make real representation of complex vectors
        Mat ret_VL          = init_VLR(has_init_L, N, M, *data.m_init_L);
        Mat ret_VR          = init_VLR(has_init_R, N, M, *data.m_init_R);

        ret_VL.get_struct().reset();
        ret_VR.get_struct().reset();

        V* ptr_VL           = reinterpret_cast<V*>(ret_VL.ptr());
        V* ptr_VR           = reinterpret_cast<V*>(ret_VR.ptr());

        Integer info;

        lapack::hsein( SIDE, EIGSRC, INITV, ptr_sel, N, lap(ptr_H), LDH, lap(ptr_W), lap(ptr_VL), 
                      VLR_LD, lap(ptr_VR), VLR_LD, M, M, lap(ptr_WORK), lap(ptr_RWORK), ptr_fail_L, 
                      ptr_fail_R, &info );

        if (info < 0)
            throw error::error_general("invalid argument passed to hsein");

        data.m_VL = Matrix(ret_VL,true);
        data.m_VR = Matrix(ret_VR,true);

        if (info > 0)
        {
            Matrix fail_L   = Matrix(ret_fail_L,false);
            Matrix fail_R   = Matrix(ret_fail_R,false);

            throw error::hess_eig_failed(data.m_VL, data.m_VR, fail_L, fail_R, eval_L, eval_R, info);
        };

        return;
    };

    static void init_sel(Integer* ptr_sel, Integer N, Integer M)
    {
        for (Integer i = 0; i < M; ++i)
            ptr_sel[i] = 1;

        for (Integer i = M; i < N; ++i)
            ptr_sel[i] = 0;
    };

    static Mat init_VLR(bool has_init, Integer N, Integer M, const Matrix& init_mat)
    {
        if (has_init == false)
            return Mat(ti::ti_type<V>(), N, M);

        value_code this_vc  = matrix_traits::value_code<V>::value;
        mat_code this_mc    = matrix_traits::get_matrix_type(this_vc, struct_code::struct_dense);
        mat_code init_mc    = init_mat.get_matrix_code();        

        if (init_mc == this_mc)
            return init_mat.get_impl<Mat>().make_unique();

        Matrix mat_VLR = convert(init_mat, Mat::matrix_code);
        return mat_VLR.get_impl_unique<Mat>();
    };
};

//--------------------------------------------------------------------------------
//                                  DENSE
//--------------------------------------------------------------------------------
template<class V> struct lapack_hrd_trd_maker
{
    using VR        = typename md::real_type<V>::type;
    using Mat       = raw::Matrix<V,struct_dense>;
    using de_Mat    = raw::Matrix<VR,struct_dense>;

    static void eval_nonsym(Integer N, Mat& Ac, Mat& tau)
    {
        Integer info1 = 0;
        Mat work1(Ac.get_type(), 1, 1);
    
        lapack::gehrd(N, 1, N, lap(Ac.ptr()), Ac.ld(), lap(tau.ptr()), lap(work1.ptr()), -1, lap(&info1));
        Integer lwork1 = (Integer) real(work1.ptr()[0]);

        work1.reset_unique(lwork1, 1);
            
        lapack::gehrd(N, 1, N, lap(Ac.ptr()), Ac.ld(), lap(tau.ptr()),  lap(work1.ptr()), lwork1, lap(&info1));

        if (info1 != 0)
        {
            std::ostringstream msg;
            msg << "lapack routine gehrd returned info = " << info1;
            throw error::error_general(msg.str());
        }
    }

    static void eval_sym(Integer N, Mat& Ac, Mat& tau, de_Mat& DE)
    {
        Integer info1 = 0;
        Mat work1(Ac.get_type(), 1, 1);
    
        hess_un_or_he_sy_selector<V>::eval_trd('L', N, Ac, tau, DE, work1, -1, &info1);
        Integer lwork1 = (Integer) real(work1.ptr()[0]);

        work1.reset_unique(lwork1, 1);
            
        hess_un_or_he_sy_selector<V>::eval_trd('L', N, Ac, tau, DE, work1, lwork1, &info1);  

        if (info1 != 0)
        {
            std::ostringstream msg;
            msg << "lapack routine syhrd/hehrd returned info = " << info1;
            throw error::error_general(msg.str());
        }
    }
};

template<class V, class S> 
struct hess_str {};

template<class V>
struct hess_str<V, struct_dense>
{
    using VR        = typename md::real_type<V>::type;
    using Mat       = raw::Matrix<V, struct_dense>;
    using de_Mat    = raw::Matrix<VR, struct_dense>;
    using i_Mat     = raw::Matrix<Integer, struct_dense>;
    
    static void eval(Matrix& ret, const Mat& A)
    {
        unitary_matrix ret_U;
        return eval(A, ret_U, ret, false);
    };

    static void eval(const Mat& A, unitary_matrix& ret_U, Matrix& ret_H, bool with_U = true)
    {
        bool is_symher  = hess_un_or_he_sy_selector<V>::is_sym_her(A);

        if (is_symher == true)
            return eval_sym(A, ret_U, ret_H, with_U);

        Integer N   = A.cols();
        Mat tau(A.get_type(),N-1,1);

        Mat Ac      = A.make_unique();
        
        Ac.set_struct(struct_flag());

        lapack_hrd_trd_maker<V>::eval_nonsym(N, Ac, tau);

        if (with_U == true)
        {
            bool isv = Ac.all_finite() && tau.all_finite();

            if (isv == false)
            {
                ret_U = unitary_matrix::from_nan(N,N,matrix_traits::value_code<V>::value);
            }
            else
            {
                Mat Qc = Ac.copy(); 

                using householder_impl = std::shared_ptr<householder_q<V>>;
                householder_impl Q_impl(new details::householder_q<V>(N, Qc, tau, N, 1));
                ret_U = unitary_matrix(Q_impl);
            };
        };        

        matcl::raw::inplace::make_triu<Mat>::eval(ret_H,Ac,-1);

        ret_H.add_struct(predefined_struct_type::hessu);

        return;
    };

    static void eval_sym(const Mat& A, unitary_matrix& ret_U, Matrix& ret_H, bool with_U)
    {
        Integer N   = A.cols();
        Mat tau(A.get_type(),N-1,1);

        Mat Ac      = A.make_unique();

        Ac.set_struct(struct_flag());

        de_Mat DE(A.get_type(), N, 3);
        lapack_hrd_trd_maker<V>::eval_sym(N, Ac, tau, DE);        

        if (with_U)
        {
            bool isv    = Ac.all_finite() && tau.all_finite();

            if (isv == false)
            {
                ret_U = unitary_matrix::from_nan(N,N,matrix_traits::value_code<V>::value);
            }
            else
            {
                using householder_impl = std::shared_ptr<householder_q<V>>;
                householder_impl Q_impl(new details::householder_q<V>(N, Ac, tau, N, 1));
                ret_U = unitary_matrix(Q_impl);
            };
        };               

        VR* ptr_DE      = DE.ptr();
        Integer DE_ld   = DE.ld();

        {
            VR* ptr_E1  = ptr_DE + 1*DE_ld;
            VR* ptr_E2  = ptr_DE + 2*DE_ld;

            for (Integer i = 0; i < N-1; ++i)
            {
                ptr_E2[i]   = ptr_E1[i];
            };
        };

        Matrix d        = (mat_row(), 0, -1, 1);
        ret_H           = bdiags(Matrix(DE,false), d, N, N);

        ret_H.add_struct(predefined_struct_type::hessu);
        ret_H.add_struct(predefined_struct_type::hessl);
        ret_H.add_struct(predefined_struct_type::sym);

        return;
    };

    static void eval_eig(const Mat& A, hess_eig_data& data)
    {
        bool is_symher  = A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));

        if (is_symher)
        {
            Matrix mat_A    = Matrix(A,false);
            Integer LD      = get_ld(mat_A, 1);
            Integer UD      = get_ud(mat_A, 1);

            // can call tridiagonal solver
            if (LD == 1 && UD == 1)
            {
                Matrix D0   = mat_A.diag(0);
                Matrix D1   = mat_A.diag(1);

                hess_str<V,struct_banded>::eval_tridiagonal(D0, D1, data);
                return;
            };
        };

        return hess_eig_dense<V>::eval(A, data);
    };
};

template<class Val>
struct gen_hess_dense
{
    using Mat   = raw::Matrix<Val,struct_dense>;
    using VR    = typename real_type<Val>::type;

    static void eval2(const Mat& A, const Mat& B, hess_gen2_ret& ret)
    {
        Integer N   = A.rows();
        Mat Ac      = A.make_unique();
        Mat Bc      = B.make_unique();

        Ac.get_struct().reset();
        Bc.get_struct().reset();

        Val* ptr_A  = Ac.ptr();
        Val* ptr_B  = Bc.ptr();
        Integer LDA = Ac.ld();
        Integer LDB = Bc.ld();

        Mat Q(A.get_type(), N, N);
        Mat Z(A.get_type(), N, N);

        Val* ptr_Q  = Q.ptr();
        Val* ptr_Z  = Z.ptr();
        Integer LDQ = Q.ld();
        Integer LDZ = Z.ld();

        Integer info = 0;

        Val work_query;
        Integer iwork_query;

        lapack::gghrdf("V", "V", N, lap(ptr_A), LDA, lap(ptr_B), LDB, lap(ptr_Q), LDQ, 
                      lap(ptr_Z), LDZ, lap(&work_query), -1, &iwork_query, -1, info);

        Integer lwork   = (Integer)real(work_query);
        Integer liwork  = iwork_query;

        using workspace     = md::workspace2<Val>;
        using iworkspace    = md::workspace2<Integer>;
        workspace WORK      = workspace(ti::ti_empty(), lwork);
        iworkspace IWORK    = iworkspace(ti::ti_empty(), liwork);

        Val* ptr_WORK       = WORK.ptr();
        Integer* ptr_IWORK  = IWORK.ptr();

        lapack::gghrdf("V", "V", N, lap(ptr_A), LDA, lap(ptr_B), LDB, lap(ptr_Q), LDQ, 
                      lap(ptr_Z), LDZ, lap(ptr_WORK), lwork, ptr_IWORK, liwork, info);

        if (info != 0)
            throw error::error_general("invalid argument is passed to gghrdf");

        Matrix ret_A    = Matrix(Ac,false);
        Matrix ret_B    = Matrix(Bc,false);

        ret_A.set_struct(predefined_struct_type::hessu);
        ret_B.set_struct(predefined_struct_type::triu);

        struct_flag sf_u;
        sf_u.set_user(unitary_flag());

        Q.set_struct(sf_u);
        Z.set_struct(sf_u);

        ret = hess_gen2_ret(ret_A, ret_B, Matrix(Q,true), Matrix(Z,true));
        return;
    };

    static void eval(const Mat& A, const Mat& B, hess_gen_ret& ret)
    {
        Integer N   = A.rows();
        Mat Ac      = A.make_unique();
        Mat Bc      = B.make_unique();

        Ac.get_struct().reset();
        Bc.get_struct().reset();

        Val* ptr_A  = Ac.ptr();
        Val* ptr_B  = Bc.ptr();
        Integer LDA = Ac.ld();
        Integer LDB = Bc.ld();

        Val* ptr_Q  = nullptr;
        Val* ptr_Z  = nullptr;
        Integer LDQ = 1;
        Integer LDZ = 1;

        Integer info = 0;

        Val work_query;
        Integer iwork_query;

        lapack::gghrdf("N", "N", N, lap(ptr_A), LDA, lap(ptr_B), LDB, lap(ptr_Q), LDQ, 
                      lap(ptr_Z), LDZ, lap(&work_query), -1, &iwork_query, -1, info);

        Integer lwork   = (Integer)real(work_query);
        Integer liwork  = iwork_query;

        using workspace     = md::workspace2<Val>;
        using iworkspace    = md::workspace2<Integer>;
        workspace WORK      = workspace(ti::ti_empty(), lwork);
        iworkspace IWORK    = iworkspace(ti::ti_empty(), liwork);

        Val* ptr_WORK       = WORK.ptr();
        Integer* ptr_IWORK  = IWORK.ptr();

        lapack::gghrdf("N", "N", N, lap(ptr_A), LDA, lap(ptr_B), LDB, lap(ptr_Q), LDQ, 
                      lap(ptr_Z), LDZ, lap(ptr_WORK), lwork, ptr_IWORK, liwork, info);

        if (info != 0)
            throw error::error_general("invalid argument is passed to gghrdf");

        Matrix ret_A    = Matrix(Ac,true);
        Matrix ret_B    = Matrix(Bc,true);

        ret_A.set_struct(predefined_struct_type::hessu);
        ret_B.set_struct(predefined_struct_type::triu);

        ret = hess_gen_ret(ret_A, ret_B);
        return;
    };
};

template<class V>
struct hess_str<V, struct_sparse>
{
    using Mat = raw::Matrix<V,struct_sparse>;

    static void eval_sym(const Mat& A, unitary_matrix& ret_U, Matrix& ret_H, bool with_U)
    {
        Matrix mA       = Matrix(A,false);
        permvec p       = order_rcm(mA);
        
        Matrix mA2      = mA(p,p);
        mA2.add_struct(predefined_struct_type::her);

        Integer LD      = get_ld(mA2, -1);
        Integer N       = A.rows();
        Integer LD_max  = (Integer)(linalg_optim_params::crossover_bd_hess_sym() * N);

        if (LD > LD_max)
        {
            using MC    = raw::Matrix<V,struct_dense>;
            MC Ac       = raw::converter<MC,Mat>::eval(A);
            hess_str<V,struct_dense>::eval_sym(Ac, ret_U, ret_H, with_U);
            return;
        };

        using MC        = raw::Matrix<V,struct_banded>;
        const Mat& A2   = mA2.get_impl<Mat>();
        MC Ac           = raw::converter<MC,Mat>::eval(A2);
        
        hess_str<V,struct_banded>::eval_sym(Ac, ret_U, ret_H, with_U);

        if (with_U == false)
            return;

        if (ret_U.all_finite() == false)
            return;

        Matrix U        = ret_U.to_matrix()(p.invperm(), colon());
        ret_U           = unitary_matrix(U, false);
    };

    static void eval(Matrix& ret, const Mat& A)
    {
        bool is_symher  = A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));

        if (is_symher)
        {
            unitary_matrix U;
            return eval_sym(A, U, ret, false);
        };

        using MC = raw::Matrix<V,struct_dense>;
        MC Ac = raw::converter<MC,Mat>::eval(A);
        hess_str<V,struct_dense>::eval(ret, Ac);
    };

    static void eval(const Mat& A, unitary_matrix& ret_U, Matrix& ret_H)
    {
        bool is_symher  = A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));

        if (is_symher)
            return eval_sym(A, ret_U, ret_H, true);

        using MC    = raw::Matrix<V,struct_dense>;
        MC Ac       = raw::converter<MC,Mat>::eval(A);
        hess_str<V,struct_dense>::eval(Ac, ret_U, ret_H);
    };

    static void eval_eig(const Mat& A, hess_eig_data& data)
    {
        bool is_symher  = A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));

        if (is_symher)
        {
            Matrix mat_A    = Matrix(A,false);
            Integer LD      = get_ld(mat_A, 1);

            // can call tridiagonal solver
            if (LD <= 1)
            {
                Matrix D0   = mat_A.diag(0);
                Matrix D1   = mat_A.diag(1);

                hess_str<V,struct_banded>::eval_tridiagonal(D0, D1, data);
                return;
            };
        };

        using MC    = raw::Matrix<V,struct_dense>;
        MC Ac       = raw::converter<MC,Mat>::eval(A);
        hess_str<V,struct_dense>::eval_eig(Ac, data);
    };
};

template<class V>
struct hess_str<V, struct_banded>
{
    using VR        = typename md::real_type<V>::type;
    using Mat       = raw::Matrix<V,struct_banded>;
    using Mat_D     = raw::Matrix<V,struct_dense>;
    using de_Mat    = raw::Matrix<VR, struct_dense>;

    static void eval(Matrix& ret, const Mat& A)
    {
        bool is_symher  = A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));

        if (is_symher)
        {
            unitary_matrix U;
            return eval_sym(A, U, ret, false);
        };

        using MC = raw::Matrix<V,struct_dense>;
        MC Ac = raw::converter<MC,Mat>::eval(A);
        hess_str<V,struct_dense>::eval(ret, Ac);
    };

    static void eval(const Mat& A, unitary_matrix& ret_U, Matrix& ret_H)
    {
        bool is_symher  = A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));

        if (is_symher)
            return eval_sym(A, ret_U, ret_H, true);

        using MC    = raw::Matrix<V,struct_dense>;
        MC Ac       = raw::converter<MC,Mat>::eval(A);
        hess_str<V,struct_dense>::eval(Ac, ret_U, ret_H);
    };

    static void eval_sym(const Mat& A, unitary_matrix& ret_U, Matrix& ret_H, bool with_U)
    {
        Mat Ac          = A.make_unique();
        Ac.get_struct().reset();

        if (Ac.has_diag(0) == false)
            throw error::band_matrix_with_main_diag_required(Ac.first_diag(), Ac.last_diag());

        Integer N       = Ac.cols();
        Integer A_LD    = Ac.number_subdiagonals();    
        Integer ld_A    = Ac.ld();

        de_Mat DE(A.get_type(), N, 3);
        Mat_D Q         = (with_U)? Mat_D(A.get_type(), N, N) :  Mat_D(A.get_type(), 1, 1);

        Integer ld_Q    = Q.ld();

        using VTR_pod       = matcl::pod_type<V>;
        using workspace     = matcl::pod_workspace<VTR_pod>;
        workspace WORK      = workspace(N);
        V* ptr_WORK         = reinterpret_cast<V*>(WORK.ptr());

        VR* ptr_DE          = DE.ptr();
        Integer DE_LD       = DE.ld();
        VR* ptr_D           = ptr_DE;
        VR* ptr_E           = ptr_DE + DE_LD;

        V* ptr_A            = Ac.rep_ptr() + Ac.first_elem_diag(0);
        V* ptr_Q            = Q.ptr();

        const char* U       = (with_U)? "V" : "N";
        Integer info        = 0;
        lapack::sbtrd(U, "L", N, A_LD, lap(ptr_A), ld_A, lap(ptr_D), lap(ptr_E), lap(ptr_Q), 
                      ld_Q, lap(ptr_WORK), info);
            
        if (info != 0)
        {
            std::ostringstream msg;
            msg << "lapack routine sbtrd returned info = " << info;
            throw error::error_general(msg.str());
        }

        if (with_U)
        {
            Matrix mat_Q    = Matrix(Q,true);

            struct_flag sf_u;
            sf_u.set_user(unitary_flag());

            mat_Q.set_struct(sf_u);

            ret_U           = unitary_matrix(mat_Q, true);
        };
        
        {
            VR* ptr_E1  = ptr_DE + 1*DE_LD;
            VR* ptr_E2  = ptr_DE + 2*DE_LD;

            for (Integer i = 0; i < N-1; ++i)
            {
                ptr_E2[i]   = ptr_E1[i];
            };
        };

        Matrix d        = (mat_row(), 0, -1, 1);
        ret_H           = bdiags(Matrix(DE,false), d, N, N);

        ret_H.add_struct(predefined_struct_type::hessu);
        ret_H.add_struct(predefined_struct_type::hessl);
        ret_H.add_struct(predefined_struct_type::sym);

        return;
    };

    static void eval_eig(const Mat& A, hess_eig_data& data)
    {
        bool is_symher  = A.get_struct().is_hermitian(A.rows() == A.cols(), is_real_matrix(A));

        if (is_symher)
        {
            Matrix mat_A    = Matrix(A,false);
            Integer LD      = get_ld(mat_A, 1);
            Integer UD      = get_ud(mat_A, 1);

            // can call tridiagonal solver
            if (LD == 1 && UD == 1)
            {
                Matrix D0   = mat_A.diag(0);
                Matrix D1   = mat_A.diag(1);

                hess_str<V,struct_banded>::eval_tridiagonal(D0, D1, data);
                return;
            };
        };

        using MC    = raw::Matrix<V,struct_dense>;
        MC Ac       = raw::converter<MC,Mat>::eval(A);
        hess_str<V,struct_dense>::eval_eig(Ac, data);
    };

    static void eval_tridiagonal(const Matrix& D0, const Matrix& D1, hess_eig_data& data)
    {
        return hess_eig_tridiag<V>::eval(D0, D1, data);
    };
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct hess_impl
{
    using M = raw::Matrix<V, S>;

    static void eval(Matrix& ret, const M& A)
    {
        return hess_str<V,S>::eval(ret, A);
    };
    static void eval(const M& A, unitary_matrix& ret_U, Matrix& ret_H)
    {
        return hess_str<V,S>::eval(A, ret_U, ret_H);
    };

    static void eval_eig(const M& A, hess_eig_data& data)
    {
        return hess_str<V,S>::eval_eig(A, data);
    };
};

template<class S>
struct hess_impl<Integer, S>
{
    using M = raw::Matrix<Integer, S>;
    static void eval(Matrix& ret, const M& A)
    {
        using MC    = raw::Matrix<Real,struct_dense>;
        MC AC       = raw::converter<MC,M>::eval(A);
        return hess_impl<Real,struct_dense>::eval(ret, AC);
    };

    static void eval(const M& A, unitary_matrix& ret_U, Matrix& ret_H)
    {
        using MC    = raw::Matrix<Real,struct_dense>;
        MC AC       = raw::converter<MC,M>::eval(A);
        hess_impl<Real,struct_dense>::eval(AC, ret_U, ret_H);
    };

    static void eval_eig(const M& A, hess_eig_data& data)
    {
        using MC    = raw::Matrix<Real,struct_dense>;
        MC AC       = raw::converter<MC,M>::eval(A);
        hess_impl<Real,struct_dense>::eval_eig(AC, data);
    };
};

template<class S>
struct hess_impl<Object, S>
{
    using M = raw::Matrix<Object, S>;
    static void eval(Matrix&, const M&)
    {
        throw error::object_value_type_not_allowed("hess");
    };
    static void eval(const M& , unitary_matrix& , Matrix&)
    {
        throw error::object_value_type_not_allowed("hess");
    };
    static void eval_eig(const M&, hess_eig_data&)
    {
        throw error::object_value_type_not_allowed("hess");
    };
};

template<class V1, class S>
struct gen_hess2_impl{};

template<class V1>
struct gen_hess2_impl<V1,struct_dense>
{
    using Mat1  = raw::Matrix<V1,struct_dense>;

    static void eval(const Mat1& A, const Mat1& B, hess_gen2_ret& ret)
    {
        return gen_hess_dense<V1>::eval2(A, B, ret);
    };
};

template<>
struct gen_hess2_impl<Object,struct_dense>
{
    using Mat1  = raw::Matrix<Object,struct_dense>;

    static void eval(const Mat1&, const Mat1&, hess_gen2_ret&)
    {
        throw error::object_value_type_not_allowed("hess_gen2");
    };
};

template<>
struct gen_hess2_impl<Integer,struct_dense>
{
    using Mat1  = raw::Matrix<Integer,struct_dense>;
    using Mat_R = raw::Matrix<Real,struct_dense>;

    static void eval(const Mat1& A, const Mat1& B, hess_gen2_ret& ret)
    {
        Mat_R AR    = raw::converter<Mat_R,Mat1>::eval(A);
        Mat_R BR    = raw::converter<Mat_R,Mat1>::eval(B);
        return gen_hess2_impl<Real,struct_dense>::eval(AR,BR,ret);
    };
};

template<class V1, class S1>
struct gen_hess_impl{};

template<class V1>
struct gen_hess_impl<V1,struct_dense>
{
    using Mat1  = raw::Matrix<V1,struct_dense>;

    static void eval(const Mat1& A, const Mat1& B, hess_gen_ret& ret)
    {
        return gen_hess_dense<V1>::eval(A, B, ret);
    };
};

template<>
struct gen_hess_impl<Object,struct_dense>
{
    using Mat1  = raw::Matrix<Object,struct_dense>;

    static void eval(const Mat1&, const Mat1&, hess_gen_ret&)
    {
        throw error::object_value_type_not_allowed("hess_gen");
    };
};

template<>
struct gen_hess_impl<Integer,struct_dense>
{
    using Mat1  = raw::Matrix<Integer,struct_dense>;
    using Mat_R = raw::Matrix<Real,struct_dense>;

    static void eval(const Mat1& A, const Mat1& B, hess_gen_ret& ret)
    {
        Mat_R AR    = raw::converter<Mat_R,Mat1>::eval(A);
        Mat_R BR    = raw::converter<Mat_R,Mat1>::eval(B);
        return gen_hess_impl<Real,struct_dense>::eval(AR,BR,ret);
    };
};

struct unary_visitor_hess1 : public extract_type_switch<void,unary_visitor_hess1,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret)
    {
        using constants::nan;
        
        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        using VR    = typename md::real_type_int_real<V>::type;

        bool isv    = mat.all_finite();
        
        if (isv == false)
        {
            ret = details::make_nan_matrix<VR>(mat.rows(), mat.cols());
            return;
        };

        return details::hess_impl<V, S>::eval(ret, mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, Matrix& ret)
    {
        using constants::nan;
        bool isv = is_finite(v);
        
        using VR = typename md::real_type_int_real<T>::type;

        if (isv == false)
        {
            ret = nan<VR>();
            return;
        }

        ret = handle;
        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&)
    {
        throw error::object_value_type_not_allowed("hess");
    }
    static void eval_scalar(const Matrix&, const Object&, Matrix&)
    {
        throw error::object_value_type_not_allowed("hess");
    }
};

struct unary_visitor_hess2 : public extract_type_switch<void, unary_visitor_hess2,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, hess_ret& ret)
    {
        using constants::nan;

        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        using VR    = typename md::real_type_int_real<V>::type;
        
        bool isv    = mat.all_finite();

        if (isv == false)
        {            
            Matrix mat_nan      = make_nan_matrix<VR>(mat.rows(), mat.cols());
            unitary_matrix Q    = unitary_matrix::from_nan(mat.rows(), mat.cols(), mat_nan.get_value_code());
            ret = hess_ret(Q, mat_nan);
            return;
        };

        unitary_matrix u;
        Matrix h;
        details::hess_impl<V,S>::eval(mat, u, h);

        ret = hess_ret(u, h) ;
        return;
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, hess_ret& ret)
    {
        using constants::nan;
        bool isv = is_finite(v);

        using VR    = typename md::real_type_int_real<T>::type;

        if (isv == false)
        {
            unitary_matrix Q = unitary_matrix::from_nan(1,1,matrix_traits::value_code<VR>::value);
            ret = hess_ret(Q, nan<VR>());
            return;
        }

        ret = hess_ret(unitary_matrix(VR(1.0),false), handle); // (u, h)
        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, hess_ret&)
    {
        throw error::object_value_type_not_allowed("hess");
    }

    static void eval_scalar(const Matrix&, const Object&, hess_ret&)
    {
        throw error::object_value_type_not_allowed("hess");
    }
};

struct hess_gen2_vis : public extract_type2_switch<void,hess_gen2_vis, mr::val_type_corrector_diag_dense>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, hess_gen2_ret& ret)
    {
        using constants::nan;

        using V1    = typename T1::value_type;
        using S1    = typename T1::struct_type;

        using V2    = typename T2::value_type;
        using S2    = typename T2::struct_type;

        using VR1   = typename md::unify_types<V1,V2>::type;
        using VR2   = typename md::unify_types<VR1,Float>::type;
        using VR    = typename md::real_type<VR2>::type;

        bool isv    = A.all_finite() && B.all_finite();

        if (isv == false)
        {
            Integer N   = A.rows();
            Matrix mat_nan  = make_nan_matrix_not_object<VR>(N, N, "hess_gen2");
            ret = hess_gen2_ret(mat_nan,mat_nan,mat_nan,mat_nan);
            return;
        };

        return details::gen_hess2_impl<V1,S1>::eval(A,B,ret);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, hess_gen2_ret& ret)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<Dense_1,Dense_2>(Ac,Bc, ret);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& A, const T2& B, hess_gen2_ret& ret)
    {
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<T1,Dense_2>(A,Bc, ret);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& A, const T2& B, hess_gen2_ret& ret)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<Dense_1,T2>(Ac,B, ret);
    };
};

struct hess_gen_vis : public extract_type2_switch<void,hess_gen_vis, mr::val_type_corrector_diag_dense>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, hess_gen_ret& ret)
    {
        using constants::nan;

        using V1    = typename T1::value_type;
        using S1    = typename T1::struct_type;

        using V2    = typename T2::value_type;
        using S2    = typename T2::struct_type;

        using VR1   = typename md::unify_types<V1,V2>::type;
        using VR2   = typename md::unify_types<VR1,Float>::type;
        using VR    = typename md::real_type<VR2>::type;

        bool isv    = A.all_finite() && B.all_finite();

        if (isv == false)
        {
            Integer N   = A.rows();
            Matrix mat_nan  = make_nan_matrix_not_object<VR>(N, N, "hess_gen2");
            ret = hess_gen_ret(mat_nan,mat_nan);
            return;
        };

        return details::gen_hess_impl<V1,S1>::eval(A,B,ret);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, hess_gen_ret& ret)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<Dense_1,Dense_2>(Ac,Bc, ret);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& A, const T2& B, hess_gen_ret& ret)
    {
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<T1,Dense_2>(A,Bc, ret);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& A, const T2& B, hess_gen_ret& ret)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<Dense_1,T2>(Ac,B, ret);
    };
};

struct hess_eig_vis : public extract_type_switch<void, hess_eig_vis,true>
{
    template<class V>
    static value_code make_value_code_E(value_code v_E)
    {
        using VR    = typename md::real_type_int_real<V>::type;

        value_code this_vc      = matrix_traits::value_code<V>::value;
        value_code this_vc_re   = matrix_traits::value_code<VR>::value;
        value_code v_E_re       = matrix_traits::real_value_type(v_E);
        bool is_compl_E         = matrix_traits::is_float_complex(v_E);

        bool is_compl           = md::is_complex<V>::value;

        if (is_compl == true)
            return this_vc;

        if (v_E_re == this_vc_re)
            return v_E;
        
        if (is_compl_E)
            return this_vc_re;

        value_code v_ret        = matrix_traits::complex_value_type(this_vc_re);
        return v_ret;
    };

    template<class T>
    static void eval(const Matrix&, const T& mat, hess_eig_data& data)
    {
        using constants::nan;
        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        using VR    = typename md::real_type_int_real<V>::type;

        //convert value types; conversion should be performed before casting,
        //because casting can generate inf values
        value_code v_E0     = data.m_E->get_value_code();
        value_code v_iL0    = data.m_init_L->get_value_code();
        value_code v_iR0    = data.m_init_R->get_value_code();
        value_code v_E      = make_value_code_E<V>(v_E0);
        value_code v_iL     = make_value_code_E<V>(v_iL0);
        value_code v_iR     = make_value_code_E<V>(v_iR0);

        bool has_init_L     = data.m_init_L->is_empty() == false;
        bool has_init_R     = data.m_init_R->is_empty() == false;

        Matrix loc_E, loc_iL, loc_iR;

        if (v_E != v_E0)
        {
            mat_code m_E    = matrix_traits::get_matrix_type(v_E, struct_code::struct_dense);
            loc_E           = convert(*data.m_E, m_E);
            data.m_E        = &loc_E;
        };

        if (v_iL0 != v_iL && has_init_L)
        {
            mat_code m_iL   = matrix_traits::get_matrix_type(v_iL, struct_code::struct_dense);
            loc_iL          = convert(*data.m_init_L, m_iL);
            data.m_init_L   = &loc_iL;
        };

        if (v_iR0 != v_iR && has_init_R)
        {
            mat_code m_iR   = matrix_traits::get_matrix_type(v_iR, struct_code::struct_dense);
            loc_iR          = convert(*data.m_init_R, m_iR);
            data.m_init_R   = &loc_iR;
        };

        bool isv    = mat.all_finite();
        bool isve   = data.m_E->all_finite();
        bool isv_iL = (has_init_L == false) || data.m_init_L->all_finite();
        bool isv_iR = (has_init_R == false) || data.m_init_R->all_finite();

        if (isv == false || isve == false || isv_iL == false || isv_iR == false)
        {
            Integer N       = mat.rows();
            Integer M       = data.m_E->length();
            Matrix mat_nan  = make_nan_matrix<VR>(N, M);

            data.m_VL = mat_nan;
            data.m_VR = mat_nan;
            return;
        };

        details::hess_impl<V,S>::eval_eig(mat, data);
        return;
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, hess_eig_data& data)
    {
        using Mat = raw::Matrix<T,struct_dense>;
        Mat m(ti::get_ti(v), v, 1, 1);
        return eval<Mat>(handle, m, data);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, hess_eig_data&)
    {
        throw error::object_value_type_not_allowed("hess_eig");
    }

    static void eval_scalar(const Matrix&, const Object&, hess_eig_data&)
    {
        throw error::object_value_type_not_allowed("hess_eig");
    }
};

static void hess_eig_impl(const Matrix& H, const Matrix& E, Matrix& VL, Matrix& VR, 
                          const Matrix& init_L, const Matrix& init_R, bool eval_L, bool eval_R)
{
    // test input arguments
    if (!H.is_square())
        throw error::square_matrix_required(H.rows(), H.cols());

    Integer H_ldiags    = get_ld(H,1);

    if (H_ldiags > 1)
        throw error::error_upper_hess_matrix_required();

    Integer N = H.rows();

    if (E.is_vector() == false || E.length() > N)
        throw error::error_hess_eig_invalid_eig(E.rows(), E.cols(), N);

    Integer M   = E.length();

    if (init_L.is_empty() == false)
    {
        if (init_L.rows() != N || init_L.cols() != N)
            throw error::error_hess_eig_invalid_init(init_L.rows(), init_L.cols(), N, M, true);
    };

    if (init_R.is_empty() == false)
    {
        if (init_R.rows() != N || init_R.cols() != N)
            throw error::error_hess_eig_invalid_init(init_R.rows(), init_R.cols(), N, M, false);
    };

    if (eval_L == true && eval_R == true)
    {
        //if evalating both sides, then both initial vectors must be empty or both nonempty
        if (init_L.is_empty() == false && init_R.is_empty() == true)
            throw error::error_hess_eig_invalid_init(init_R.rows(), init_R.cols(), N, M, false);

        if (init_R.is_empty() == false && init_L.is_empty() == true)
            throw error::error_hess_eig_invalid_init(init_L.rows(), init_L.cols(), N, M, true);
    };

    // fast exit
    if (M == 0)
    {
        VL  = spzeros(M,N,0,H.get_value_code());
        VR  = VL;
        return;
    };

    hess_eig_data hed(E, VL, VR, init_L, init_R, eval_L, eval_R);
    return hess_eig_vis::make<const Matrix&>(H, hed);
};

};};

namespace matcl
{

static void hess_impl(Matrix& ret, const Matrix& A)
{
    if (!A.is_square())
        throw error::error_hess_nonsq();

    if (A.structural_nnz() == 0  || get_ld(A,1) <= 1)
    {
        ret = A;
        ret.get_struct().add_ldiags(struct_flag::one);
        return;
    };

    if (get_ud(A,-1) <= 1)
    {
        Matrix ret_A = A(colon(colon::cend(), -1, 1), colon(colon::cend(), -1, 1));
        ret = ret_A;
        ret.get_struct().add_ldiags(struct_flag::one);
        return;
    };

    return details::unary_visitor_hess1::make<const Matrix&>(A,ret);
};

Matrix matcl::hess(const Matrix& A0)
{
    //increase refcount;
    Matrix A(A0);

    Matrix ret;
    hess_impl(ret,A);
    return ret;
};
Matrix matcl::hess(Matrix&& A0)
{
    Matrix A(std::move(A0));

    Matrix ret;
    hess_impl(ret,A);
    return ret;
};

static void hess2_impl(hess_ret& ret, const Matrix& A)
{
    if (!A.is_square())
        throw error::error_hess_nonsq();

    if (A.structural_nnz() == 0 || get_ld(A,1) <= 1)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        A.get_struct().add_ldiags(struct_flag::one);
        ret = hess_ret(unitary_matrix(speye(A.rows(),A.rows(), vt),false), A);
        return;
    }

    if(get_ud(A,-1) <= 1)
    {
        Matrix ret_A = A(colon(colon::cend(), -1, 1), colon(colon::cend(), -1, 1));

        matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        // reversed eye
        Matrix I        = speye(A.rows(),A.rows(), vt);
        Matrix ret_U    = I(colon(), colon(colon::cend(), -1, 1)); 

        struct_flag sf;
        sf.set_user(unitary_flag());

        ret_U.add_struct(sf);
        
        ret_A.get_struct().add_ldiags(struct_flag::one);
        ret = hess_ret(unitary_matrix(ret_U,false), ret_A);
        return;
    }

    return details::unary_visitor_hess2::make<const Matrix&>(A,ret);
};

hess_ret matcl::hess2(const Matrix& A0)
{
    //increase refcount
    Matrix A(A0);

    hess_ret ret;
    hess2_impl(ret,A);
    return ret;
};
hess_ret matcl::hess2(Matrix&& A0)
{
    //increase refcount
    Matrix A(std::move(A0));

    hess_ret ret;
    hess2_impl(ret,A);
    return ret;
};

static void hess_gen_impl(hess_gen_ret& ret, const Matrix& A, const Matrix& B)
{
    Integer N_A = A.rows();
    Integer N_B = B.rows();

    if (A.cols() != N_A || B.cols()!= N_B || N_A != N_B)
        throw error::error_gen_hess(A.rows(), A.cols(), B.rows(), B.cols());

    Integer N       = N_A;
    Integer A_ld    = get_ld(A,1);
    Integer B_ld    = get_ld(A,0);

    if (A_ld <= 1 && B_ld == 0)
    {
        //nothing to do
        ret = hess_gen_ret(A,B);
        return;
    };

    Integer A_ud    = get_ud(A,1);
    Integer B_ud    = get_ud(A,0);

    if (A_ud <= 1 && B_ud == 0)
    {
        //just revert matrices
        Matrix A2   = A(colon(N,-1,1),colon(N,-1,1));
        Matrix B2   = B(colon(N,-1,1),colon(N,-1,1));

        ret = hess_gen_ret(A2,B2);
        return;
    };

    return details::hess_gen_vis::make(A,B,ret);
};

hess_gen_ret matcl::hess_gen(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    hess_gen_ret ret;
    hess_gen_impl(ret,A,B);
    return ret;
}
hess_gen_ret matcl::hess_gen(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    hess_gen_ret ret;
    hess_gen_impl(ret,A,B);
    return ret;
}
hess_gen_ret matcl::hess_gen(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    hess_gen_ret ret;
    hess_gen_impl(ret,A,B);
    return ret;
}
hess_gen_ret matcl::hess_gen(Matrix&& A0, Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    hess_gen_ret ret;
    hess_gen_impl(ret,A,B);
    return ret;
}

static void hess_gen2_impl(hess_gen2_ret& ret, const Matrix& A, const Matrix& B)
{
    Integer N_A = A.rows();
    Integer N_B = B.rows();

    if (A.cols() != N_A || B.cols()!= N_B || N_A != N_B)
        throw error::error_gen_hess(A.rows(), A.cols(), B.rows(), B.cols());

    Integer N       = N_A;
    Integer A_ld    = get_ld(A,1);
    Integer B_ld    = get_ld(A,0);

    if (A_ld <= 1 && B_ld == 0)
    {
        matcl::value_code vtA   = matrix_traits::real_value_type(A.get_value_code());
        matcl::value_code vtB   = matrix_traits::real_value_type(A.get_value_code());
        matcl::value_code vtAB  = matrix_traits::unify_value_types(vtA, vtB);
        matcl::value_code vt    = matrix_traits::unify_value_types(vtAB, value_code::v_float);

        //nothing to do

        Matrix Q    = speye(N,N,vt);
        Matrix Z    = Q;

        ret = hess_gen2_ret(A,B,Q,Z);
        return;
    };

    Integer A_ud    = get_ud(A,1);
    Integer B_ud    = get_ud(A,0);

    if (A_ud <= 1 && B_ud == 0)
    {
        matcl::value_code vtA   = matrix_traits::real_value_type(A.get_value_code());
        matcl::value_code vtB   = matrix_traits::real_value_type(A.get_value_code());
        matcl::value_code vtAB  = matrix_traits::unify_value_types(vtA, vtB);
        matcl::value_code vt    = matrix_traits::unify_value_types(vtAB, value_code::v_float);

        //just revert matrices
        Matrix A2   = A(colon(N,-1,1),colon(N,-1,1));
        Matrix B2   = B(colon(N,-1,1),colon(N,-1,1));

        // reversed eye
        Matrix I        = speye(N, N, vt);
        Matrix Q        = I(colon(), colon(colon::cend(), -1, 1)); 
        Matrix Z        = Q;

        ret = hess_gen2_ret(A2,B2,Q,Z);
        return;
    };

    return details::hess_gen2_vis::make(A,B,ret);
};

hess_gen2_ret matcl::hess_gen2(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    hess_gen2_ret ret;
    hess_gen2_impl(ret,A,B);
    return ret;
};
hess_gen2_ret matcl::hess_gen2(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    hess_gen2_ret ret;
    hess_gen2_impl(ret,A,B);
    return ret;
};
hess_gen2_ret matcl::hess_gen2(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    hess_gen2_ret ret;
    hess_gen2_impl(ret,A,B);
    return ret;
};
hess_gen2_ret matcl::hess_gen2(Matrix&& A0, Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    hess_gen2_ret ret;
    hess_gen2_impl(ret,A,B);
    return ret;
};

Matrix matcl::hess_left_eig(const Matrix& H, const Matrix& E, const Matrix& init_L0)
{
    Matrix init_L(init_L0);

    Matrix VL, VR;
    Matrix init_R = zeros(0,0);
    details::hess_eig_impl(H, E, VL, VR, init_L, init_R, true, false);
    return VL;
};
Matrix matcl::hess_left_eig(const Matrix& H, const Matrix& E, Matrix&& init_L0)
{
    Matrix init_L(std::move(init_L0));

    Matrix VL, VR;
    Matrix init_R = zeros(0,0);
    details::hess_eig_impl(H, E, VL, VR, init_L, init_R, true, false);
    return VL;
};
Matrix matcl::hess_right_eig(const Matrix& H, const Matrix& E, const Matrix& init_R0)
{
    Matrix init_R(init_R0);

    Matrix VL, VR;
    Matrix init_L = zeros(0,0);
    details::hess_eig_impl(H, E, VL, VR, init_L, init_R, false, true);
    return VR;
};
Matrix matcl::hess_right_eig(const Matrix& H, const Matrix& E, Matrix&& init_R0)
{
    Matrix init_R(std::move(init_R0));

    Matrix VL, VR;
    Matrix init_L = zeros(0,0);
    details::hess_eig_impl(H, E, VL, VR, init_L, init_R, false, true);
    return VR;
};
mat_tup_2 matcl::hess_eig(const Matrix& H, const Matrix& E, const Matrix& init_L0, const Matrix& init_R0)
{
    Matrix init_L(init_L0);
    Matrix init_R(init_R0);

    Matrix VL, VR;
    details::hess_eig_impl(H, E, VL, VR, init_L, init_R, true, true);
    return mat_tup_2(VL,VR);
};

mat_tup_2 matcl::hess_eig(const Matrix& H, const Matrix& E, const Matrix& init_L0, Matrix&& init_R0)
{
    Matrix init_L(init_L0);
    Matrix init_R(std::move(init_R0));

    Matrix VL, VR;
    details::hess_eig_impl(H, E, VL, VR, init_L, init_R, true, true);
    return mat_tup_2(VL,VR);
};

mat_tup_2 matcl::hess_eig(const Matrix& H, const Matrix& E, Matrix&& init_L0, const Matrix& init_R0)
{
    Matrix init_L(std::move(init_L0));
    Matrix init_R(init_R0);

    Matrix VL, VR;
    details::hess_eig_impl(H, E, VL, VR, init_L, init_R, true, true);
    return mat_tup_2(VL,VR);
};

mat_tup_2 matcl::hess_eig(const Matrix& H, const Matrix& E, Matrix&& init_L0, Matrix&& init_R0)
{
    Matrix init_L(std::move(init_L0));
    Matrix init_R(std::move(init_R0));

    Matrix VL, VR;
    details::hess_eig_impl(H, E, VL, VR, init_L, init_R, true, true);
    return mat_tup_2(VL,VR);
};

hess_ret matcl::hess2triu(const Matrix& H, bool from_left)
{
    return from_left? qr_givens(H) : rq_givens(H);
};

hess_ret matcl::hess2triu(Matrix&& H, bool from_left)
{
    return from_left? qr_givens(std::move(H)) : rq_givens(std::move(H));
};

hess_ret matcl::lhess2tril(const Matrix& H, bool from_left)
{
    return from_left? ql_givens(H) : lq_givens(H);
};

hess_ret matcl::lhess2tril(Matrix&& H, bool from_left)
{
    return from_left? ql_givens(std::move(H)) : lq_givens(std::move(H));
};

linsolve_obj matcl::linsolve_hess(const Matrix& A, const Matrix& U, const Matrix& H,
                                  const options& opts)
{
    return linsolve_hess(A, unitary_matrix(U,true), H, opts);
};

linsolve_obj matcl::linsolve_hess(const Matrix& A, const unitary_matrix& U, const Matrix& H,
                                  const options& opts)
{
    if (A.rows() != U.rows())
        throw error::invalid_hess_factors();
    if (U.rows() != U.cols())
        throw error::invalid_hess_factors();
    if (U.cols() != H.rows())
        throw error::invalid_hess_factors();
    if (H.rows() != H.cols())
        throw error::invalid_hess_factors();

    Integer ld  = matcl::get_ld(H,1);
    if (ld > 1)
        throw error::invalid_hess_factors();

    Integer N   = U.rows();

    if (U.rows() == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(H.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, U.get_type())));
    };

    bool isv            = H.all_finite() && U.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(U.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(H.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, H.get_type() )));
    };

    using data_ptr = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_hess(A, U, H, opts)));
};

};