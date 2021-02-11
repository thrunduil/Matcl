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

#include "matcl-linalg/decompositions/householder_q.h"
#include "matcl-linalg/decompositions/qr_utils.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-linalg/utils/optim_params.h"
#include "matcl-linalg/general/linalg_exception.h"

#pragma warning(push)
#pragma warning(disable:4127) //conditional expression is constant

namespace matcl { namespace details
{

//------------------------------------------------------------------------
//                          lapack_xyygqr_maker
//------------------------------------------------------------------------
template<class Val> 
struct qr_un_or_selector{};

template<> 
struct qr_un_or_selector<Real>
{
    using Mat = raw::Matrix<Real, struct_dense>;

    static void eval(Integer M, Integer N, Integer K, Mat& Qc, const Mat& tau, Mat& work, 
                        Integer lwork, Integer* info2, Integer offset)
    {
        lapack::orgqr(M-offset, N-offset, K, lap(Qc.ptr() + offset), Qc.ld(), 
                        lap(tau.ptr()), lap(work.ptr()), lwork, lap(info2));
    };
};

template<> 
struct qr_un_or_selector<Float>
{
    using Mat = raw::Matrix<Float, struct_dense>;

    static void eval(Integer M, Integer N, Integer K, Mat& Qc, const Mat& tau, Mat& work, 
                        Integer lwork, Integer* info2, Integer offset)
    {
        lapack::orgqr(M-offset, N-offset, K, lap(Qc.ptr() + offset), 
                      Qc.ld(), lap(tau.ptr()), lap(work.ptr()), lwork, lap(info2));
    };
};

template<>
struct qr_un_or_selector<Complex>
{
    using Mat = raw::Matrix<Complex, struct_dense>;

    static void eval(Integer M, Integer N, Integer K, Mat& Qc, const Mat& tau, Mat& work, 
                        Integer lwork, Integer* info2, Integer offset)
    {
        lapack::ungqr(M-offset, N-offset, K, lap(Qc.ptr() + offset), 
                      Qc.ld(), lap(tau.ptr()), lap(work.ptr()), lwork, lap(info2));
    };
};

template<>
struct qr_un_or_selector<Float_complex>
{
    using Mat = raw::Matrix<Float_complex, struct_dense>;

    static void eval(Integer M, Integer N, Integer K, Mat& Qc, const Mat& tau, Mat& work, 
                        Integer lwork, Integer* info2, Integer offset)
    {
        lapack::ungqr(M-offset, N-offset, K, lap(Qc.ptr() + offset), 
                      Qc.ld(), lap(tau.ptr()), lap(work.ptr()), lwork, lap(info2));
    };
};

template<class V> 
struct lapack_xyygqr_maker
{
    using Mat = raw::Matrix<V, struct_dense>;

    static void eval(Integer M, Integer N, Integer K, Mat& Qc, const Mat& tau, Integer offset)
    {
        Integer info2 = 0;
        Mat work2(Qc.get_type(), 1, 1);

        qr_un_or_selector<V>::eval(M, N, K, Qc, tau, work2, -1, &info2, offset);

        Integer lwork2 = (Integer) real(work2.ptr()[0]);
        work2.reset_unique(lwork2, 1);
        
        qr_un_or_selector<V>::eval(M, N, K, Qc, tau, work2, lwork2, &info2, offset);
        
        if (info2 != 0)
        {
            std::ostringstream msg;
            msg << "ungqr/orgqr returned info = " << info2;

            throw error::error_general(msg.str());
        };
    }
};

//eval W = C*V
template<class Val_ret, class Val_1>
struct eval_gemv
{
    static void eval(Integer M, Integer K, const Val_ret* C, Integer C_ld, const Val_1* V, Val_ret* W)
    {
        if (K == 0 || M == 0)
            return;

        for (Integer i = 0; i < M; ++i)
            W[i]    = C[i];

        if (K == 1)
            return;

        C += C_ld;

        for (Integer j = 1; j < K; ++j)
        {
            Val_ret val = Val_ret(V[j]);

            for (Integer i = 0; i < M; ++i)
            {
                W[i]    = W[i] + C[i] * val;
            };

            C += C_ld;
        }
    };
};

template<class Val>
struct eval_gemv<Val,Val>
{
    static void eval(Integer M, Integer K, const Val* C, Integer C_ld, const Val* V, Val* W)
    {
        if (K == 0 || M == 0)
            return;

        for (Integer i = 0; i < M; ++i)
            W[i]    = C[i];

        C += C_ld;

        Val one(1.0);

        if (K > 1)
            lapack::gemv("N",M,K-1,*lap(&one),lap(C),C_ld, lap(V+1), 1, *lap(&one), lap(W), 1);
    };
};

// eval C := C + alpha * V * W'
template<class Val_ret, class Val_1>
struct eval_gerc
{
    static void eval(Integer M, Integer K, const Val_ret& alpha, const Val_ret* W, const Val_1* V, 
                     Val_ret* C, Integer C_ld)
    {
        if (M == 0 || K == 0)
            return;

        for (Integer i = 0; i < M; ++i)
        {
            C[i]    = C[i] + alpha * W[i];
        };

        C += C_ld;

        for (Integer j = 1; j < K; ++j)
        {
            Val_ret val = alpha * Val_ret(conj(V[j]));

            for (Integer i = 0; i < M; ++i)
            {
                C[i]    = C[i] + val * W[i];
            };

            C += C_ld;
        };
    };
};

template<class Val>
struct eval_gerc<Val,Val>
{
    static void eval(Integer M, Integer K, const Val& alpha, const Val* W, const Val* V, Val* C, Integer C_ld)
    {
        if (M == 0 || K == 0)
            return;

        for (Integer i = 0; i < M; ++i)
        {
            C[i]    = C[i] + alpha * W[i];
        };

        C += C_ld;

        if (K > 1)
            lapack::gerc(M,K-1,*lap(&alpha),lap(W),1,lap(V+1),1,lap(C),C_ld);
    };
};

template<class Val1, class Val_ret>
struct eval_larfb
{
    // workspace has size LD_work x IB, where LD_work >= NI if SIDE == 'L' or LD_work >= MI otherwise
    // ptr_work_convert has size IBxIB + LD2_work * IB, where LD2_work >= MI if SIDE == 'L' 
    // or LD2_work >= NI otherwise
    static void eval(const char* SIDE, const char* TRANS, const char* Forward, 
                     Integer MI, Integer NI, Integer IB, const Val1* ptr_V, Integer V_ld,
                     const Val1* ptr_T, Integer LDT, Val_ret* ptr_C, Integer C_ld,
                     Val_ret* ptr_work, Integer LD_work, Val_ret* ptr_work_convert)
    {
        bool is_real    = md::is_float_real_scalar<Val1>::value;
        bool is_real2   = md::is_float_real_scalar<Val_ret>::value;
        bool make_conj  = (is_real == false) && (TRANS[0] == 'T' || TRANS[0] == 't');

        Val_ret* T2     = ptr_work_convert;
        Val_ret* T2_s   = T2;
        Integer T2_ld   = IB;
        Integer T2_size = IB * IB;

        if (make_conj)
        {
            for (Integer j = 0; j < IB; ++j)
            {
                for (Integer i = 0; i <= j; ++i)
                {
                    T2[i]   = Val_ret(conj(ptr_T[i]));
                };

                T2      += T2_ld;
                ptr_T   += LDT;
            };
        }
        else
        {
            for (Integer j = 0; j < IB; ++j)
            {
                for (Integer i = 0; i <= j; ++i)
                {
                    T2[i]   = Val_ret(ptr_T[i]);
                };

                T2      += T2_ld;
                ptr_T   += LDT;
            };
        };

        Integer V2_ld;

        if (SIDE[0] == 'l' || SIDE[0] == 'L')
            V2_ld       = std::max(1,MI);
        else
            V2_ld       = std::max(1,NI);

        Val_ret* V2     = ptr_work_convert + T2_size;
        Val_ret* V2_s   = V2;        

        Integer V2_size = V2_ld * IB;
        (void)V2_size;

        if (make_conj)
        {
            for (Integer j = 0; j < V2_ld; ++j)
            {
                for (Integer i = j+1; i < IB; ++i)
                {
                    V2[i]   = Val_ret(conj(ptr_V[i]));
                };
                V2      += V2_ld;
                ptr_V   += V_ld;
            };
        }
        else
        {
            for (Integer j = 0; j < IB; ++j)
            {
                for (Integer i = j+1; i < V2_ld; ++i)
                {
                    V2[i]   = Val_ret(ptr_V[i]);
                };
                V2      += V2_ld;
                ptr_V   += V_ld;
            };
        };

        const char* TRANS2;
        if (TRANS[0] == 'N' || TRANS[0] == 'n')
            TRANS2 = "N";
        else if (is_real2 == true)
            TRANS2 = "T";
        else
            TRANS2 = "C";

        eval_larfb<Val_ret,Val_ret>::eval(SIDE, TRANS2, Forward, MI, NI, IB, V2_s, V2_ld,
                     T2_s, T2_ld, ptr_C, C_ld, ptr_work, LD_work, nullptr);
    };
};

template<class Val>
struct eval_larfb<Val,Val>
{
    static void make_conj_mat(Integer MI, Integer NI, Integer IB, const char* SIDE, const Val* ptr_T, 
                Integer LDT, const Val* ptr_V, Integer V_ld, Val*& T2_s, Integer& T2_ld, Val*& V2_s,
                Integer& V2_ld, Val* ptr_work_convert)
    {
        Val* T2         = ptr_work_convert;
        T2_s            = T2;
        T2_ld           = IB;
        Integer T2_size = IB * IB;

        for (Integer j = 0; j < IB; ++j)
        {
            for (Integer i = 0; i <= j; ++i)
            {
                T2[i]   = conj(ptr_T[i]);
            };

            T2      += T2_ld;
            ptr_T   += LDT;
        };

        if (SIDE[0] == 'l' || SIDE[0] == 'L')
            V2_ld       = std::max(1,MI);
        else
            V2_ld       = std::max(1,NI);

        Val* V2         = ptr_work_convert + T2_size;
        V2_s            = V2;        

        Integer V2_size = V2_ld * IB;
        (void)V2_size;

        for (Integer j = 0; j < IB; ++j)
        {
            for (Integer i = j+1; i < V2_ld; ++i)
            {
                V2[i]   = conj(ptr_V[i]);
            };

            V2      += V2_ld;
            ptr_V   += V_ld;
        };
    };

    // workspace has size LD_work x IB, where LD_work >= NI if SIDE == 'L' or LD_work >= MI otherwise
    // ptr_work_convert has size IBxIB + LD2_work * IB, where LD2_work >= MI if SIDE == 'L' 
    // or LD2_work >= NI otherwise if complex conjugation is required, otherwise not referred;
    // complex conjugation is required if Val is complex type and TRANS == 'T'
    static void eval(const char* SIDE, const char* TRANS, const char* Forward, 
                     Integer MI, Integer NI, Integer IB, const Val* ptr_V, Integer V_ld,
                     const Val* ptr_T, Integer LDT, Val* ptr_C, Integer C_ld, Val* ptr_work, Integer LD_work,
                     Val* ptr_work_convert)
    {
        bool make_conj = (md::is_float_real_scalar<Val>::value == false) 
                        && (TRANS[0] == 'T' || TRANS[0] == 't');

        using LV       = typename lapack_value_type<Val>::type;

        if (make_conj == true)
        {
            Val* T2;
            Val* V2;
            Integer T2_ld;
            Integer V2_ld;

            make_conj_mat(MI, NI, IB, SIDE, ptr_T, LDT, ptr_V, V_ld, T2, T2_ld, V2, V2_ld, ptr_work_convert);

            return lapack::larfb<LV>(SIDE, "C", Forward, "Columnwise", MI, NI, IB, lap(V2), V2_ld,
                     lap(T2), T2_ld, lap(ptr_C), C_ld, lap(ptr_work), LD_work);
        }
        else
        {
            return lapack::larfb<LV>(SIDE, TRANS, Forward, "Columnwise", MI, NI, IB, lap(ptr_V), V_ld,
                     lap(ptr_T), LDT, lap(ptr_C), C_ld, lap(ptr_work), LD_work);
        };
    };
};

template<class Val1, class Val2, class Struct>
struct householder_band_mult_struct{};

template<class Val1, class Val2, class Struct>
struct householder_mult_struct{};

//------------------------------------------------------------------------
//                          DENSE
//------------------------------------------------------------------------
template<class Val1, class Val2>
struct householder_mult_struct<Val1, Val2, struct_dense>
{
    using Mat           = raw::Matrix<Val2,struct_dense>;
    using umatrix       = householder_q<Val1>;

    static const bool is_real_umatrix   = md::is_float_real_scalar<Val1>::value;

    static void eval_block(Matrix& ret, bool left0, Integer NB, const Mat& X, 
                           const umatrix& data, trans_type t_unitary)
    {
        //block version based on DORMQR; ORMQR cannot be used drectly because ORMQR 
        //modifies X
        using VTR       = typename md::unify_types<Val1,Val2>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;
        using LV1       = typename lapack_value_type<Val1>::type;

        static const bool ISR = md::is_float_real_scalar<Val1>::value;

        Integer XM, XN, Ret_M, Ret_N;
        if (left0 == false)
        {
            XM      = data.reflector_length();
            XN      = X.cols();

            if (t_unitary == trans_type::no_trans)
            {
                Ret_M   = data.rows();
                Ret_N   = XN;
            }
            else
            {
                Ret_M   = data.cols();
                Ret_N   = XN;
            };
        }
        else
        {
            if (t_unitary == trans_type::no_trans)
            {
                XM      = X.rows();
                XN      = data.reflector_length();
                Ret_M   = XM;
                Ret_N   = data.m_mat_cols;
            }
            else
            {
                XM      = X.rows();
                XN      = data.reflector_length();
                Ret_M   = XM;
                Ret_N   = data.rows();
            };
        };

        ret_type C          = make_copy<VTR, Val2>::eval(X, XM, XN);
        VTR* ptr_C          = C.ptr();
        Integer C_ld        = C.ld();

        const Val1* ptr_V   = data.m_reflectors.ptr();
        Integer V_ld        = data.m_reflectors.ld();
        const Val1* ptr_tau = data.m_tau_vec.ptr();
        Integer offset      = data.m_offset;

        bool lap_left       = !left0;

        //apply offset
        Integer M           = C.rows();
        Integer N           = C.cols();        
        Integer K           = data.number_reflectors();

        if (lap_left == true)
        {
            M               = M - offset;
            ptr_V           = ptr_V + offset;
            ptr_C           = ptr_C + offset;
        }
        else
        {
            N               = N - offset;
            ptr_V           = ptr_V + offset;
            ptr_C           = ptr_C + offset * C_ld;
        };
        
        Integer I1, I2, I3;
        bool NOTRAN         = (t_unitary == trans_type::no_trans);        

        if ( (lap_left == true && NOTRAN == false) || ( lap_left == false && NOTRAN == true ) )
        {
            I1  = 1;
            I2  = K;
            I3  = NB;
        }
        else
        {
            I1  = ( ( K-1 ) / NB )*NB + 1;
            I2  = 1;
            I3  = -NB;
        };

        Integer NI = 0, JC = 0, MI = 0, IC = 0;
        Integer NQ = 0, NW = 0;

        if (lap_left == true)
        {
            NI  = N;
            JC  = 1;
            NQ  = M;
            NW  = N;
        }
        else
        {
            MI  = M;
            IC  = 1;
            NQ  = N;
            NW  = M;
        };

        const char* SIDE    = (lap_left == true)? "L" : "R";
        const char* TRANS;

        if (t_unitary == trans_type::no_trans)
            TRANS = "N";
        else if (t_unitary == trans_type::trans || ISR)
            TRANS = "T";
        else
            TRANS = "C";

        //local array
        static const Integer NBMAX  = linalg_optim_params::max_block_size_ORMQR;
        static const Integer LDT    = NBMAX+1;

        Val1 local_T[LDT * NBMAX];

        Integer LD_work     = (lap_left == true)? N : M;
        Integer LD2_work    = (lap_left == true)? M : N;
        Integer work_size   = LD_work * NB;
        Integer work2_size  = NB*NB + LD2_work * NB;
        bool make_conj      = (t_unitary == trans_type::trans) 
                            && md::is_float_real_scalar<Val1>::value == false;

        if ((make_conj == false) && std::is_same<Val1,Val2>::value)
            work2_size      = 0;

        using VTR_pod       = pod_type<VTR>;
        using workspace     = pod_workspace<VTR_pod>;
        workspace WORK      = matcl::pod_workspace<VTR_pod>(work_size);
        workspace WORK2     = matcl::pod_workspace<VTR_pod>(work2_size);
        VTR* ptr_w          = reinterpret_cast<VTR*>(WORK.ptr());
        VTR* ptr_w2         = reinterpret_cast<VTR*>(WORK2.ptr());

        for (Integer I = I1; ; I += I3)
        {
            if (I3 >= 0 && !(I <= I2))
                break;
            if (I3 < 0 && !(I >= I2))
                break;

            Integer IB  = std::min(NB, K-I+1 );

            // Form the triangular factor of the block reflector
            // H = H(i) H(i+1) . . . H(i+ib-1)

            lapack::larft<LV1>("Forward", "Columnwise", NQ-I+1, IB, lap(ptr_V + (I-1) + (I-1)*V_ld), 
                          V_ld, lap(ptr_tau + (I-1)), lap(local_T), LDT );

            if (lap_left == true)
            {
                // H or H**T is applied to C(i:m,1:n)
                MI  = M - I + 1;
                IC  = I;
            }
            else
            {
                // H or H**T is applied to C(1:m,i:n)
               NI   = N - I + 1;
               JC   = I;
            };

            // Apply H or H**T
            eval_larfb<Val1, VTR>::eval(SIDE, TRANS, "Forward", MI, NI, IB, 
                        ptr_V + (I-1) + (I-1)*V_ld, V_ld, local_T, LDT, 
                        ptr_C + IC-1 +(JC-1)*C_ld, C_ld, ptr_w, LD_work, ptr_w2);
        };

        C.assign_to_fresh(C.resize(Ret_M,Ret_N));

        struct_flag su;
        if (data.rows() == data.cols())
            su.set_user(unitary_flag());

        bool is_square  = C.rows() == C.cols();
        bool is_sq_X    = X.rows() == X.cols();

        if (left0 == false)
            C.set_struct(predefined_struct_ext::mult_struct(su,X.get_struct(),t_unitary,trans_type::no_trans,
                                is_real_umatrix, is_real_matrix(X), true, is_sq_X, is_square));
        else
            C.set_struct(predefined_struct_ext::mult_struct(X.get_struct(), su, trans_type::no_trans,t_unitary,
                                is_real_matrix(X), is_real_umatrix, is_sq_X, true, is_square));

        ret = Matrix(C,true);
        return;
    };

    static void eval_right(Matrix& ret, const Mat& X, const umatrix& data, trans_type t_unitary)
    {
        using VTR       = typename md::unify_types<Val1,Val2>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;

        Integer ret_N   = X.cols();
        Integer ret_M;

        if (t_unitary == trans_type::no_trans)
            ret_M       = data.rows();
        else
            ret_M       = data.cols();

        Integer num_refl    = data.number_reflectors();
        Integer refl_length = data.reflector_length();
        Integer V_ldiags    = data.m_ldiags;

        Integer NB          = linalg_optim_params::block_size_ORMQR();

        if (NB <= num_refl && V_ldiags * 2 >= NB )
        {
            return eval_block(ret, false, NB, X, data, t_unitary);
        };

        Integer V_ld        = data.m_reflectors.ld();
        const Val1* ptr_tau = data.m_tau_vec.ptr();

        ret_type Y          = make_copy<VTR, Val2>::eval(X, refl_length, ret_N);
        VTR* ptr_Y          = Y.ptr();
        Integer Y_ld        = Y.ld();

        const Val1* ptr_V0  = data.m_reflectors.ptr();

        //apply offset
        Integer offset      = data.m_offset;
        ptr_V0              = ptr_V0 + offset;
        ptr_Y               = ptr_Y + offset;

        if (t_unitary == trans_type::no_trans)
        {
            ptr_V0          = ptr_V0 + (num_refl - 1) * V_ld;

            for (Integer j = 0; j < ret_N; ++j)
            {
                const Val1* ptr_V       = ptr_V0;

                for (Integer i = num_refl - 1; i >= 0; --i)
                {
                    //apply i-th elementary reflector H(i) = I - tau v*v'

                    Val1 tau            = ptr_tau[i];

                    if (tau != Val1(0))
                    {
                        //eval dot v'*Y

                        VTR dot         = ptr_Y[i];

                        Integer max_k   = std::min(refl_length, i + V_ldiags + 1 ) - offset;

                        for (Integer k = i + 1; k < max_k; ++k)
                            dot         = dot + conj(ptr_V[k]) * ptr_Y[k];

                        // form Y := Y - dot * tau * v
                        VTR dot_t       = dot * tau;

                        ptr_Y[i]        = ptr_Y[i] - dot_t;

                        for (Integer k = i + 1; k < max_k; ++k)
                            ptr_Y[k]    = ptr_Y[k] - dot_t*ptr_V[k];
                    };

                    ptr_V   -= V_ld;
                };

                ptr_Y   += Y_ld;
            };

            struct_flag su;
            if (data.cols() == data.rows())
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(su,X.get_struct(),t_unitary,trans_type::no_trans,
                                is_real_umatrix, is_real_matrix(X), true, is_sq_X, is_square));

            ret = Matrix(Y,true);

            return;
        }
        else
        {
            bool conj_U     = (t_unitary == trans_type::conj_trans);

            for (Integer j = 0; j < ret_N; ++j)
            {
                const Val1* ptr_V       = ptr_V0;

                for (Integer i = 0; i < num_refl; ++i)
                {
                    //apply i-th elementary reflector H(i) = I - tau v*v'

                    Val1 tau            = conj_U? conj(ptr_tau[i]) : ptr_tau[i] ;

                    if (tau != Val1(0))
                    {
                        //eval dot v'*Y

                        VTR dot         = ptr_Y[i];
                        Integer max_k   = std::min(refl_length, i + V_ldiags + 1 ) - offset;

                        if (conj_U == true)
                        {
                            for (Integer k = i + 1; k < max_k; ++k)
                                dot         = dot + conj(ptr_V[k]) * ptr_Y[k];

                            // form Y := Y - dot * tau * v
                            VTR dot_t       = dot * tau;

                            ptr_Y[i]        = ptr_Y[i] - dot_t;

                            for (Integer k = i + 1; k < max_k; ++k)
                                ptr_Y[k]    = ptr_Y[k] - dot_t*ptr_V[k];
                        }
                        else
                        {
                            for (Integer k = i + 1; k < max_k; ++k)
                                dot         = dot + ptr_V[k] * ptr_Y[k];

                            // form Y := Y - dot * tau * v
                            VTR dot_t       = dot * tau;

                            ptr_Y[i]        = ptr_Y[i] - dot_t;

                            for (Integer k = i + 1; k < max_k; ++k)
                                ptr_Y[k]    = ptr_Y[k] - dot_t*conj(ptr_V[k]);
                        };
                    };

                    ptr_V   += V_ld;
                };

                ptr_Y   += Y_ld;
            };

            Y.assign_to_fresh(Y.resize(ret_M,ret_N));

            struct_flag su;
            if (data.cols() == data.rows())
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(su,X.get_struct(),t_unitary,trans_type::no_trans,
                                is_real_umatrix, is_real_matrix(X), true, is_sq_X, is_square));

            ret = Matrix(Y,true);

            return;
        };
    };
    static void eval_left(Matrix& ret, const Mat& X, const umatrix& data, trans_type t_unitary)
    {
        using VTR       = typename md::unify_types<Val1,Val2>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;

        Integer ret_M   = X.rows();        
        Integer ret_N;

        if (t_unitary == trans_type::no_trans)
            ret_N       = data.cols();
        else
            ret_N       = data.rows();

        Integer num_refl    = data.number_reflectors();
        Integer refl_length = data.reflector_length();
        Integer V_ld        = data.m_reflectors.ld();
        Integer V_ldiags    = data.m_ldiags;

        Integer NB          = linalg_optim_params::block_size_ORMQR();

        if (NB <= num_refl && V_ldiags * 2 >= NB )
            return eval_block(ret, true, NB, X, data, t_unitary);
        
        const Val1* ptr_tau = data.m_tau_vec.ptr();        

        using VTR_pod       = pod_type<VTR>;
        using workspace     = pod_workspace<VTR_pod>;
        workspace WORK      = matcl::pod_workspace<VTR_pod>(ret_M);
        VTR* ptr_w          = reinterpret_cast<VTR*>(WORK.ptr());

        ret_type Y          = make_copy<VTR, Val2>::eval(X, ret_M, data.reflector_length());
        VTR* ptr_Y          = Y.ptr();
        Integer Y_ld        = Y.ld();
        const Val1* ptr_V0  = data.m_reflectors.ptr();

        //apply offset
        Integer offset      = data.m_offset;
        ptr_V0              = ptr_V0 + offset;
        ptr_Y               = ptr_Y + offset * Y_ld;

        if (t_unitary == trans_type::no_trans)
        {            
            const Val1* ptr_V   = ptr_V0;
            for (Integer i = 0; i < num_refl; ++i)
            {
                //apply i-th elementary reflector H(i) = I - tau v*v'

                Val1 tau            = ptr_tau[i];

                if (tau != Val1(0))
                {
                    Integer max_k   = std::min(refl_length, i + V_ldiags + 1 ) - offset;
                    Integer K2      = max_k - i;

                    //eval w = Y*v
                    eval_gemv<VTR, Val1>::eval(ret_M, K2, ptr_Y + i*Y_ld, Y_ld, ptr_V + i, ptr_w);

                    // eval Y := Y - tau * w * v'
                    eval_gerc<VTR, Val1>::eval(ret_M, K2, VTR(-tau), ptr_w, ptr_V + i, ptr_Y + i*Y_ld, Y_ld);
                };

                ptr_V   += V_ld;
            };

            Y.assign_to_fresh(Y.resize(ret_M,ret_N));

            struct_flag su;
            if (data.rows() == data.cols())
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),su,trans_type::no_trans,t_unitary,
                                is_real_matrix(X), is_real_umatrix, is_sq_X, true, is_square));

            ret = Matrix(Y,true);

            return;
        }
        else
        {
            bool conj_U     = (t_unitary == trans_type::conj_trans)
                            || md::is_real_scalar<Val1>::value == true;            

            const Val1* ptr_V   = ptr_V0 + (num_refl-1) * V_ld;

            if (conj_U == false)
            {
                using V1_pod        = matcl::pod_type<Val1>;
                using workspace_V   = matcl::pod_workspace<V1_pod>;
                workspace_V WORK_V  = workspace_V(refl_length);
                Val1* ptr_Vw        = reinterpret_cast<Val1*>(WORK_V.ptr());

                for (Integer i = num_refl - 1; i >= 0; --i)
                {
                    //apply i-th elementary reflector H(i) = I - tau v*v'

                    Val1 tau            = conj_U? conj(ptr_tau[i]) : ptr_tau[i] ;

                    if (tau != Val1(0))
                    {
                        Integer max_k   = std::min(refl_length, i + V_ldiags + 1 ) - offset;
                        Integer K2      = max_k - i;

                        ptr_Vw[0]       = Val1(1.0);
                        for (Integer l = 1; l < K2; ++l)
                            ptr_Vw[l]   = conj(ptr_V[i+l]);

                        //eval w = Y*v
                        eval_gemv<VTR, Val1>::eval(ret_M, K2, ptr_Y + i*Y_ld, Y_ld, ptr_Vw, ptr_w);
                    
                        // eval Y := Y - tau * w * v'
                        eval_gerc<VTR, Val1>::eval(ret_M, K2, VTR(-tau), ptr_w, ptr_Vw, ptr_Y + i*Y_ld, Y_ld);
                    };

                    ptr_V   -= V_ld;
                };
            }
            else
            {
                for (Integer i = num_refl - 1; i >= 0; --i)
                {
                    //apply i-th elementary reflector H(i) = I - tau v*v'

                    Val1 tau            = conj_U? conj(ptr_tau[i]) : ptr_tau[i] ;

                    if (tau != Val1(0))
                    {
                        Integer max_k   = std::min(refl_length, i + V_ldiags + 1 ) - offset;
                        Integer K2      = max_k - i;

                        //eval w = Y*v
                        eval_gemv<VTR, Val1>::eval(ret_M, K2, ptr_Y + i*Y_ld, Y_ld, ptr_V + i, ptr_w);

                        // eval Y := Y - tau * w * v'
                        eval_gerc<VTR, Val1>::eval(ret_M, K2, VTR(-tau), ptr_w, ptr_V + i, ptr_Y + i*Y_ld, Y_ld);
                    };

                    ptr_V   -= V_ld;
                };
            };

            struct_flag su;
            if (data.cols() == data.rows())
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(), su, trans_type::no_trans,t_unitary,
                                is_real_matrix(X), is_real_umatrix, is_sq_X, true, is_square));

            ret = Matrix(Y,true);

            return;
        };
    };
};

template<class Val1, class Val2>
struct householder_band_mult_struct<Val1, Val2, struct_dense>
{
    using Mat           = raw::Matrix<Val2,struct_dense>;
    using umatrix       = householder_band_q<Val1>;

    static const bool is_real_umatrix = md::is_float_real_scalar<Val1>::value;

    static void eval_right(Matrix& ret, const Mat& X, const umatrix& data, trans_type t_unitary)
    {
        using VTR       = typename md::unify_types<Val1,Val2>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;
        using BM        = raw::Matrix<Val1,struct_banded>;

        Integer N   = X.cols();
        Integer M, K;

        if (t_unitary == trans_type::no_trans)
        {
            M   = data.rows();
            K   = data.cols();
        }
        else
        {
            M   = data.rows();
            K   = data.cols();
        };

        Integer MK          = data.number_reflectors();
        Integer V_ld        = data.m_reflectors.ld();
        const Val1* ptr_tau = data.m_tau_vec.ptr();
        const BM& V         = data.m_reflectors;

        if (V.has_diag(0) == false)
            throw error::band_matrix_with_main_diag_required(V.first_diag(), V.last_diag());

        Integer V_ldiags    = V.number_subdiagonals();

        if (t_unitary == trans_type::no_trans)
        {
            ret_type Y          = make_copy<VTR, Val2>::eval(X, M, X.cols());

            VTR* ptr_Y          = Y.ptr();
            Integer Y_ld        = Y.ld();
            const Val1* ptr_V0  = V.rep_ptr() + V.first_elem_diag(0) + (MK - 1) * (V_ld-1);

            for (Integer j = 0; j < N; ++j)
            {
                const Val1* ptr_V       = ptr_V0;

                for (Integer i = MK - 1; i >= 0; --i)
                {
                    //apply i-th elementary reflector H(i) = I - tau v*v'

                    Val1 tau            = ptr_tau[i];

                    if (tau != Val1(0))
                    {
                        //eval dot v'*Y

                        VTR dot         = ptr_Y[i];
                        Integer max_k   = std::min(M, i + V_ldiags + 1 );

                        for (Integer k = i + 1; k < max_k; ++k)
                            dot         = dot + conj(ptr_V[k]) * ptr_Y[k];

                        // form Y := Y - dot * tau * v
                        VTR dot_t       = dot * tau;

                        ptr_Y[i]        = ptr_Y[i] - dot_t;

                        for (Integer k = i + 1; k < max_k; ++k)
                            ptr_Y[k]    = ptr_Y[k] - dot_t*ptr_V[k];
                    };

                    ptr_V   -= V_ld-1;
                };

                ptr_Y   += Y_ld;
            };

            struct_flag su;
            if (M == K)
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(su,X.get_struct(),t_unitary,trans_type::no_trans,
                                is_real_umatrix, is_real_matrix(X), true, is_sq_X, is_square));

            ret = Matrix(Y,true);

            return;
        }
        else
        {
            bool conj_U     = (t_unitary == trans_type::conj_trans);
            ret_type Y      = make_copy<VTR, Val2>::eval(X, K, X.cols());
            VTR* ptr_Y      = Y.ptr();
            Integer Y_ld    = Y.ld();
            const Val1* ptr_V0  = V.rep_ptr() + V.first_elem_diag(0);

            for (Integer j = 0; j < N; ++j)
            {
                const Val1* ptr_V       = ptr_V0;

                for (Integer i = 0; i < MK; ++i)
                {
                    //apply i-th elementary reflector H(i) = I - tau v*v'

                    Val1 tau            = conj_U? conj(ptr_tau[i]) : ptr_tau[i] ;
                    Integer max_k       = std::min(K, i + V_ldiags + 1);

                    if (tau != Val1(0))
                    {
                        //eval dot v'*Y

                        VTR dot         = ptr_Y[i];

                        if (conj_U == true)                        
                        {
                            for (Integer k = i + 1; k < max_k; ++k)
                                dot         = dot + conj(ptr_V[k]) * ptr_Y[k];

                            // form Y := Y - dot * tau * v
                            VTR dot_t       = dot * tau;

                            ptr_Y[i]        = ptr_Y[i] - dot_t;

                            for (Integer k = i + 1; k < max_k; ++k)
                                ptr_Y[k]    = ptr_Y[k] - dot_t*ptr_V[k];
                        }
                        else
                        {
                            for (Integer k = i + 1; k < max_k; ++k)
                                dot         = dot + ptr_V[k] * ptr_Y[k];

                            // form Y := Y - dot * tau * v
                            VTR dot_t       = dot * tau;

                            ptr_Y[i]        = ptr_Y[i] - dot_t;

                            for (Integer k = i + 1; k < max_k; ++k)
                                ptr_Y[k]    = ptr_Y[k] - dot_t*conj(ptr_V[k]);
                        };
                    };

                    ptr_V   += V_ld - 1;
                };

                ptr_Y   += Y_ld;
            };

            Y.assign_to_fresh(Y.resize(M,N));

            struct_flag su;
            if (M == K)
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(su,X.get_struct(),t_unitary,trans_type::no_trans,
                                is_real_umatrix, is_real_matrix(X), true, is_sq_X, is_square));

            ret = Matrix(Y,true);

            return;
        };
    };

    static void eval_left(Matrix& ret, const Mat& X, const umatrix& data, trans_type t_unitary)
    {
        using VTR       = typename md::unify_types<Val1,Val2>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;
        using BM        = raw::Matrix<Val1,struct_banded>;

        Integer M   = X.rows();
        Integer K   = X.cols();
        Integer N;

        if (t_unitary == trans_type::no_trans)
            N       = data.cols();
        else
            N       = data.rows();

        Integer NK          = data.number_reflectors();
        Integer V_ld        = data.m_reflectors.ld();
        const Val1* ptr_tau = data.m_tau_vec.ptr();

        const BM& V         = data.m_reflectors;        

        if (V.has_diag(0) == false)
            throw error::band_matrix_with_main_diag_required(V.first_diag(), V.last_diag());

        Integer V_ldiags    = V.number_subdiagonals();

        using VTR_pod       = pod_type<VTR>;
        using workspace     = pod_workspace<VTR_pod>;
        workspace WORK      = matcl::pod_workspace<VTR_pod>(M);
        VTR* ptr_w          = reinterpret_cast<VTR*>(WORK.ptr());

        if (t_unitary == trans_type::no_trans)
        {
            ret_type Y          = make_copy<VTR, Val2>::eval(X, M, X.cols());

            VTR* ptr_Y          = Y.ptr();
            Integer Y_ld        = Y.ld();
            const Val1* ptr_V   = V.rep_ptr()+ V.first_elem_diag(0);            

            for (Integer i = 0; i < NK; ++i)
            {
                //apply i-th elementary reflector H(i) = I - tau v*v'

                Val1 tau            = ptr_tau[i];

                if (tau != Val1(0))
                {
                    Integer max_k   = std::min(K, i + V_ldiags + 1);
                    Integer K2      = max_k - i;

                    //eval w = Y*v
                    eval_gemv<VTR, Val1>::eval(M, K2, ptr_Y + i*Y_ld, Y_ld, ptr_V + i, ptr_w);

                    // eval Y := Y - tau * w * v'
                    eval_gerc<VTR, Val1>::eval(M, K2, VTR(-tau), ptr_w, ptr_V + i, ptr_Y + i*Y_ld, Y_ld);
                };

                ptr_V   += V_ld - 1;
            };

            Y.assign_to_fresh(Y.resize(M,N));

            struct_flag su;
            if (data.rows() == data.cols())
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),su,trans_type::no_trans,t_unitary,
                                is_real_matrix(X), is_real_umatrix, is_sq_X, true, is_square));

            ret = Matrix(Y,true);

            return;
        }
        else
        {
            bool conj_U     = (t_unitary == trans_type::conj_trans)
                            || md::is_real_scalar<Val1>::value == true;

            ret_type Y      = make_copy<VTR, Val2>::eval(X, M, N);

            VTR* ptr_Y      = Y.ptr();
            Integer Y_ld    = Y.ld();
            const Val1* ptr_V   = V.rep_ptr() + V.first_elem_diag(0) + (NK-1) * (V_ld-1);

            if (conj_U == false)
            {
                using V1_pod        = matcl::pod_type<Val1>;
                using workspace_V   = matcl::pod_workspace<V1_pod>;
                workspace_V WORK_V  = workspace_V(N);
                Val1* ptr_Vw        = reinterpret_cast<Val1*>(WORK_V.ptr());

                for (Integer i = NK - 1; i >= 0; --i)
                {
                    //apply i-th elementary reflector H(i) = I - tau v*v'

                    Val1 tau            = conj_U? conj(ptr_tau[i]) : ptr_tau[i] ;

                    if (tau != Val1(0))
                    {
                        Integer max_k   = std::min(N, i + V_ldiags + 1 );

                        Integer K2      = max_k - i;

                        ptr_Vw[0]       = Val1(1.0);
                        for (Integer l = 1; l < K2; ++l)
                            ptr_Vw[l]   = conj(ptr_V[i+l]);

                        //eval w = Y*v
                        eval_gemv<VTR, Val1>::eval(M, K2, ptr_Y + i*Y_ld, Y_ld, ptr_Vw, ptr_w);
                    
                        // eval Y := Y - tau * w * v'
                        eval_gerc<VTR, Val1>::eval(M, K2, VTR(-tau), ptr_w, ptr_Vw, ptr_Y + i*Y_ld, Y_ld);
                    };

                    ptr_V   -= V_ld - 1;
                };
            }
            else
            {
                for (Integer i = NK - 1; i >= 0; --i)
                {
                    //apply i-th elementary reflector H(i) = I - tau v*v'

                    Val1 tau            = conj_U? conj(ptr_tau[i]) : ptr_tau[i] ;

                    if (tau != Val1(0))
                    {
                        Integer max_k   = std::min(N, i + V_ldiags + 1 );
                        Integer K2      = max_k - i;

                        //eval w = Y*v
                        eval_gemv<VTR, Val1>::eval(M, K2, ptr_Y + i*Y_ld, Y_ld, ptr_V + i, ptr_w);
                    
                        // eval Y := Y - tau * w * v'
                        eval_gerc<VTR, Val1>::eval(M, K2, VTR(-tau), ptr_w, ptr_V + i, ptr_Y + i*Y_ld, Y_ld);
                    };

                    ptr_V   -= V_ld - 1;
                };
            };

            struct_flag su;
            if (data.rows() == data.cols())
                su.set_user(unitary_flag());

            bool is_square  = Y.rows() == Y.cols();
            bool is_sq_X    = X.rows() == X.cols();
            Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),su,trans_type::no_trans,t_unitary,
                                is_real_matrix(X), is_real_umatrix, is_sq_X, true, is_square));

            ret = Matrix(Y,true);

            return;
        };
    };
};

//------------------------------------------------------------------------
//                          SPARSE
//------------------------------------------------------------------------
template<class Val1, class Val2>
struct householder_mult_struct<Val1, Val2, struct_sparse>
{
    using Mat           = raw::Matrix<Val2,struct_sparse>;
    using umatrix       = householder_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return householder_mult_struct<Val1, Val2, struct_dense>::eval_right(ret,md,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return householder_mult_struct<Val1, Val2, struct_dense>::eval_left(ret,md,data,t_unitary);
    };
};

template<class Val1, class Val2>
struct householder_band_mult_struct<Val1, Val2, struct_sparse>
{
    using Mat           = raw::Matrix<Val2,struct_sparse>;
    using umatrix       = householder_band_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return householder_band_mult_struct<Val1, Val2, struct_dense>::eval_right(ret,md,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return householder_band_mult_struct<Val1, Val2, struct_dense>::eval_left(ret,md,data,t_unitary);
    };
};

//------------------------------------------------------------------------
//                          BAND
//------------------------------------------------------------------------
template<class Val1, class Val2>
struct householder_mult_struct<Val1, Val2, struct_banded>
{
    using Mat           = raw::Matrix<Val2,struct_banded>;
    using umatrix       = householder_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return householder_mult_struct<Val1, Val2, struct_dense>::eval_right(ret,md,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return householder_mult_struct<Val1, Val2, struct_dense>::eval_left(ret,md,data,t_unitary);
    };
};

template<class Val1, class Val2>
struct householder_band_mult_struct<Val1, Val2, struct_banded>
{
    using Mat           = raw::Matrix<Val2,struct_banded>;
    using umatrix       = householder_band_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return householder_band_mult_struct<Val1, Val2, struct_dense>::eval_right(ret,md,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return householder_band_mult_struct<Val1, Val2, struct_dense>::eval_left(ret,md,data,t_unitary);
    };
};

//------------------------------------------------------------------------
//                          householder_helper
//------------------------------------------------------------------------
template<class Val1, class Val2, class Struct>
struct householder_mult_impl
{
    using Mat           = raw::Matrix<Val2,Struct>;
    using umatrix       = householder_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        return householder_mult_struct<Val1,Val2,Struct>::eval_right(ret,mat,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        return householder_mult_struct<Val1,Val2,Struct>::eval_left(ret,mat,data,t_unitary);
    };
};

template<class Val1, class Struct>
struct householder_mult_impl<Val1, Integer, Struct>
{
    using Mat           = raw::Matrix<Integer,Struct>;
    using umatrix       = householder_q<Val1>;
    using umatrix_band  = householder_band_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using Mat_D     = raw::Matrix<Real,struct_dense>;
        Mat_D md        = raw::converter<Mat_D,Mat>::eval(mat);
        return householder_mult_impl<Val1,Real,struct_dense>::eval_right(ret, md, data, t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using Mat_D     = raw::Matrix<Real,struct_dense>;
        Mat_D md        = raw::converter<Mat_D,Mat>::eval(mat);
        return householder_mult_impl<Val1,Real,struct_dense>::eval_left(ret, md, data, t_unitary);
    };
};

template<class Val1, class Val2, class Struct>
struct householder_band_mult_impl
{
    using Mat           = raw::Matrix<Val2,Struct>;
    using umatrix       = householder_band_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        return householder_band_mult_struct<Val1,Val2,Struct>::eval_right(ret,mat,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        return householder_band_mult_struct<Val1,Val2,Struct>::eval_left(ret,mat,data,t_unitary);
    };
};

template<class Val1, class Struct>
struct householder_band_mult_impl<Val1, Integer, Struct>
{
    using Mat           = raw::Matrix<Integer,Struct>;
    using umatrix       = householder_band_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using Mat_D     = raw::Matrix<Real,struct_dense>;
        Mat_D md        = raw::converter<Mat_D,Mat>::eval(mat);
        return householder_band_mult_impl<Val1,Real,struct_dense>::eval_right(ret, md, data, t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using Mat_D     = raw::Matrix<Real,struct_dense>;
        Mat_D md        = raw::converter<Mat_D,Mat>::eval(mat);
        return householder_band_mult_impl<Val1,Real,struct_dense>::eval_left(ret, md, data, t_unitary);
    };
};

template<class Val>
struct householder_mult_right : public extract_type_switch<void, householder_mult_right<Val>,true>
{
    using umatrix   = householder_q<Val>;

    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using VM    = typename T::value_type;
        using ST    = typename T::struct_type;

        return householder_mult_impl<Val,VM,ST>::eval_right(ret, mat, data, t_unitary);
    };
    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret, t_unitary, data);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_right");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_right");
    };
};

template<class Val>
struct householder_mult_left : public extract_type_switch<void, householder_mult_left<Val>,true>
{
    using umatrix   = householder_q<Val>;

    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using VM    = typename T::value_type;
        using ST    = typename T::struct_type;

        return householder_mult_impl<Val,VM, ST>::eval_left(ret, mat, data, t_unitary);
    };
    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret, t_unitary, data);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_left");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_left");
    };
};

template<class Val>
struct householder_q_convert_vis : public extract_type_switch<void, householder_q_convert_vis<Val>,true>
{
    using umatrix   = householder_q<Val>;
    using data_ptr  = typename unitary_matrix_data::data_ptr;

    template<class T>
    static void eval(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        using VM    = typename T::value_type;

        ret = data.convert_impl<VM>();
    };
    template<class T>
    static void eval_scalar(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        ret = data.convert_impl<T>();
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::convert");
    };
    static void eval_scalar(const Matrix&, const Object&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::convert");
    };
};
template<class Val>
struct householder_band_q_convert_vis : public extract_type_switch<void, householder_band_q_convert_vis<Val>,true>
{
    using umatrix   = householder_band_q<Val>;
    using data_ptr  = typename unitary_matrix_data::data_ptr;

    template<class T>
    static void eval(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        using VM    = typename T::value_type;

        ret = data.convert_impl<VM>();
    };
    template<class T>
    static void eval_scalar(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        ret = data.convert_impl<T>();
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::convert");
    };
    static void eval_scalar(const Matrix&, const Object&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::convert");
    };
};

template<class Val>
struct householder_band_mult_right : public extract_type_switch<void, householder_band_mult_right<Val>,true>
{
    using umatrix   = householder_band_q<Val>;

    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using VM    = typename T::value_type;
        using ST    = typename T::struct_type;

        return householder_band_mult_impl<Val,VM, ST>::eval_right(ret, mat, data, t_unitary);
    };
    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret, t_unitary, data);
    };
    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_right");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_right");
    };
};

template<class Val>
struct householder_band_mult_left : public extract_type_switch<void, householder_band_mult_left<Val>,true>
{
    using umatrix   = householder_band_q<Val>;

    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using VM    = typename T::value_type;
        using ST    = typename T::struct_type;

        return householder_band_mult_impl<Val,VM,ST>::eval_left(ret, mat, data, t_unitary);
    };
    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret, t_unitary, data);
    };
    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_left");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_left");
    };
};

template<class Val>
struct householder_helper
{
    using Mat       = raw::Matrix<Val,struct_dense>;
    using umatrix   = householder_q<Val>;

    static void to_matrix(Matrix& ret, const umatrix& data);
    static void eval_mult_right(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                                const umatrix& data);
    static void eval_mult_left(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                                const umatrix& data);
    static void create_reflectors(Mat& Ac, Mat& tau);
};

template<class Val>
struct householder_band_helper
{
    using umatrix = householder_band_q<Val>;

    static void to_matrix(Matrix& ret, const umatrix& data);
    static void eval_mult_right(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                                const umatrix& data);
    static void eval_mult_left(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                                const umatrix& data);
};

template<class Val>
void householder_helper<Val>::to_matrix(Matrix& ret, const umatrix& data)
{
    using Mat = raw::Matrix<Val,struct_dense>;

    Integer M       = data.rows();
    Integer N       = data.cols();

    Integer K       = data.number_reflectors();
    Integer off     = data.m_offset;

    Mat Qc(data.m_reflectors.get_type());
    if (data.m_reflectors.rows() == M && data.m_reflectors.cols() == N)
    {
        Qc.assign_to_fresh(data.m_reflectors.copy());
    }
    else
    {
        Qc.assign_to_fresh(data.m_reflectors.resize(M,N).make_unique());
    };

    lapack_xyygqr_maker<Val>::eval(M, N, K, Qc, data.m_tau_vec, off);

    if (off != 0)
    {
        Integer LD  = Qc.ld();
        Val* ptr    = Qc.ptr() + (N-1)*LD;
        Val* ptr_p  = Qc.ptr() + (N-2)*LD;

        for (Integer j = N - 1; j >= off; --j)
        {
            for (Integer i = 0; i < off; ++i)
                ptr[i] = Val(0.0);

            for (Integer i = off; i < M; ++i)
                ptr[i] = ptr_p[i];

            ptr     -= LD;
            ptr_p   -= LD;
        };

        for (Integer j = off - 1; j >= 0; --j)
        {
            for (Integer i = 0; i < M; ++i)
                ptr[i] = Val(0.0);

            ptr[j]  = Val(1.0);
            ptr     -= LD;
        };
    };
        
    Matrix Qc_mat   = Matrix(Qc, true);
    Qc_mat.set_struct(struct_flag()); 

    if (Qc.rows() == Qc.cols())
    {
        struct_flag sf_u;
        sf_u.set_user(unitary_flag());
        Qc_mat.add_struct(sf_u); 
    };
 
    ret = Qc_mat;
    return;
};

template<class Val>
void householder_helper<Val>::create_reflectors(Mat& Ac, Mat& tau)
{
    Integer M       = Ac.rows();
    Integer N       = Ac.cols();
    Integer K       = std::min(M,N);

    Val* ptr_X      = Ac.ptr();
    Val* ptr_tau    = tau.ptr();

    Integer X_ld    = Ac.ld();
    Integer order   = M;

    for (Integer i = 0; i < K; ++i)
    {
        lapack::larfg(order, lap(ptr_X + 0), lap(ptr_X + 1), 1, lap(ptr_tau + i));
        ptr_X       += X_ld + 1;
        --order;
    };    
};


template<class Val>
void householder_band_helper<Val>::to_matrix(Matrix& ret, const umatrix& data)
{
    using Mat_B = raw::Matrix<Val,struct_banded>;
    using Mat_D = raw::Matrix<Val,struct_dense>;

    Integer K   = data.number_reflectors();
    Mat_D Qc    = raw::converter<Mat_D,Mat_B>::eval(data.m_reflectors);
    lapack_xyygqr_maker<Val>::eval(data.rows(), data.cols(), K, Qc, data.m_tau_vec, 0);
        
    Matrix Qc_mat   = Matrix(Qc, true);
    Qc_mat.set_struct(struct_flag()); 

    if (Qc.rows() == Qc.cols())
    {
        struct_flag sf_u;
        sf_u.set_user(unitary_flag());
        Qc_mat.add_struct(sf_u); 
    };
    ret = Qc_mat;
    return;
};

template<class Val>
void householder_helper<Val>::eval_mult_right(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                            const umatrix& data)
{
    return householder_mult_right<Val>::make<const Matrix&>(mat, ret, t_unitary, data);
};

template<class Val>
void householder_helper<Val>::eval_mult_left(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                            const umatrix& data)
{
    return householder_mult_left<Val>::make<const Matrix&>(mat, ret, t_unitary, data);
};

template<class Val>
void householder_band_helper<Val>::eval_mult_right(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                            const umatrix& data)
{
    return householder_band_mult_right<Val>::make<const Matrix&>(mat, ret, t_unitary, data);
};

template<class Val>
void householder_band_helper<Val>::eval_mult_left(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                            const umatrix& data)
{
    return householder_band_mult_left<Val>::make<const Matrix&>(mat, ret, t_unitary, data);
};

//------------------------------------------------------------------------
//                          householder_q
//------------------------------------------------------------------------
template<class Val>
householder_q<Val>::householder_q()
    :m_mat_cols(0), m_reflectors(ti::ti_type<Val>()), m_tau_vec(ti::ti_type<Val>())
    ,m_ldiags(0), m_offset(0)
{};
template<class Val>
householder_q<Val>::householder_q(Integer N, const Mat& Qc, const Mat& tau, Integer ldiags, Integer offset)
    : m_mat_cols(N), m_reflectors(Qc), m_tau_vec(tau), m_offset(offset)
{
    Integer M   = m_reflectors.rows();
    m_ldiags    = std::max(std::min(ldiags,M-1),1);
};

template<class Val>
householder_q<Val>::~householder_q()
{};

template<class Val>
serialization_helper<unitary_matrix_data>*
householder_q<Val>::get_serialization_helper() const
{
    std::ostringstream msg;
    msg << "householder_q" << "<" << get_type_name<Val>::eval() << ">";

    return serialization_helper<unitary_matrix_data>
            ::get<householder_q<Val>>(msg.str());
};

template<class Val>
void householder_q<Val>::mult_right(matcl::Matrix& ret, const matcl::Matrix& X, 
                                     trans_type t_unit) const
{
    if (X.is_scalar())
    {
        this->to_matrix(ret);
        ret = mmul(ret,X,t_unit);
        return;
    };

    error::check_mul(rows(), cols(), X.rows(), X.cols(), t_unit, trans_type::no_trans);

    Integer M, N;

    if (t_unit == trans_type::no_trans)
        M   = this->rows();
    else
        M   = this->cols();

    N   = X.cols();

    if (X.structural_nnz() == 0)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(X.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        ret = spzeros(M, N, 0, vt);
        return;
    };

    return householder_helper<Val>::eval_mult_right(ret, X, t_unit, *this);
};
template<class Val>
void householder_q<Val>::mult_left(matcl::Matrix& ret, const matcl::Matrix& X, 
                                     trans_type t_unit) const
{
    if (X.is_scalar())
    {
        this->to_matrix(ret);
        ret = mmul(X,ret,trans_type::no_trans,t_unit);
        return;
    };

    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t_unit);

    Integer M = X.rows();
    Integer N;

    if (t_unit == trans_type::no_trans)
        N   = this->cols();
    else
        N   = this->rows();

    if (X.structural_nnz() == 0)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(X.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        ret = spzeros(M, N, 0, vt);
        return;
    };

    return householder_helper<Val>::eval_mult_left(ret, X, t_unit, *this);
};

template<class Val>
void householder_q<Val>::to_matrix(matcl::Matrix& ret) const
{
    return householder_helper<Val>::to_matrix(ret, *this);
};

template<class Val>
typename householder_q<Val>::data_ptr 
householder_q<Val>::convert(value_code new_val_code) const
{
    value_code vc   = matrix_traits::unify_value_types(new_val_code, value_code::v_float);
    Matrix v        = zeros(0,0,vc);

    data_ptr ret;
    householder_q_convert_vis<Val>::make<const Matrix&>(v,*this, ret);
    return ret;
};

template<class Val>
template<class T>
typename householder_q<Val>::data_ptr 
householder_q<Val>::convert_impl() const
{
    using TR    = typename details::unify_types<T, Float>::type;
    using Mat_C = raw::Matrix<TR,struct_dense>;

    Mat_C ret_refl  = raw::converter<Mat_C, Mat>::eval(m_reflectors);
    Mat_C ret_tau   = raw::converter<Mat_C, Mat>::eval(m_tau_vec);

    return data_ptr(new householder_q<TR>(m_mat_cols, ret_refl, ret_tau, m_ldiags, m_offset));
};

template<class Val>
void householder_q<Val>::save(oarchive& os) const
{
    //general parameters
    os << this->m_mat_cols;
    os << this->m_ldiags;
    os << this->m_offset;    

    os.get() << m_reflectors;
    os.get() << m_tau_vec;
};

template<class Val>
unitary_matrix_data* householder_q<Val>::load(std::istream& is)
{
    using ptr_type = std::unique_ptr<householder_q<Val>>;
    ptr_type ret = ptr_type(new householder_q<Val>());

    //general parameters
    is >> ret->m_mat_cols;
    is >> ret->m_ldiags;
    is >> ret->m_offset;   

    Matrix Q, tau;
    is >> Q;
    is >> tau;

    ret->m_reflectors.assign_to_fresh(Q.impl<Mat>());
    ret->m_tau_vec.assign_to_fresh(tau.impl<Mat>());

    return ret.release();
};

template<class Val>
void householder_q<Val>::save(std::ostream& os) const
{
    os << " ";

    //general parameters
    os << this->m_mat_cols  << " ";
    os << this->m_ldiags << " ";
    os << this->m_offset << " ";   

    os << Matrix(m_reflectors,false);
    os << Matrix(m_tau_vec,false);
};

template<class Val>
unitary_matrix_data* householder_q<Val>::load(iarchive& ar)
{
    using ptr_type  = std::unique_ptr<householder_q<Val>>;
    ptr_type ret    = ptr_type(new householder_q<Val>());

    //general parameters
    ar >> ret->m_mat_cols;
    ar >> ret->m_ldiags;
    ar >> ret->m_offset;

    ar.get() >> ret->m_reflectors;
    ar.get() >> ret->m_tau_vec;

    return ret.release();
};

//------------------------------------------------------------------------
//                          householder_band_q
//------------------------------------------------------------------------
template<class Val>
householder_band_q<Val>::householder_band_q()
    :m_mat_cols(0), m_reflectors(ti::ti_type<Val>()), m_tau_vec(ti::ti_type<Val>())
{};

template<class Val>
householder_band_q<Val>::householder_band_q(Integer N, const Mat_B& Qc, const Mat_D& tau)
    :m_mat_cols(N), m_reflectors(Qc), m_tau_vec(tau)
{
    if (m_reflectors.has_diag(0) == false)
        throw error::band_matrix_with_main_diag_required(m_reflectors.first_diag(), 
                                                         m_reflectors.last_diag());
};


template<class Val>
householder_band_q<Val>::~householder_band_q()
{};

template<class Val>
serialization_helper<unitary_matrix_data>*
householder_band_q<Val>::get_serialization_helper() const
{
    std::ostringstream msg;
    msg << "householder_band_q" << "<" << get_type_name<Val>::eval() << ">";

    return serialization_helper<unitary_matrix_data>
            ::get<householder_band_q<Val>>(msg.str());
};

template<class Val>
void householder_band_q<Val>::mult_right(matcl::Matrix& ret, const matcl::Matrix& X, 
                                     trans_type t_unit) const
{
    if (X.is_scalar())
    {
        this->to_matrix(ret);
        ret = mmul(ret,X,t_unit);
        return;
    };

    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t_unit, trans_type::no_trans);

    Integer M, N;

    if (t_unit == trans_type::no_trans)
        M   = this->rows();
    else
        M   = this->cols();

    N   = X.cols();

    if (X.structural_nnz() == 0)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(X.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        ret = spzeros(M, N, 0, vt);
        return;
    };

    return householder_band_helper<Val>::eval_mult_right(ret, X, t_unit, *this);
};
template<class Val>
void householder_band_q<Val>::mult_left(matcl::Matrix& ret, const matcl::Matrix& X, 
                                     trans_type t_unit) const
{
    if (X.is_scalar())
    {
        this->to_matrix(ret);
        ret = mmul(X,ret,trans_type::no_trans,t_unit);
        return;
    };

    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t_unit);

    Integer M = X.rows();
    Integer N;

    if (t_unit == trans_type::no_trans)
        N   = this->cols();
    else
        N   = this->rows();

    if (X.structural_nnz() == 0)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(X.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        ret = spzeros(M, N, 0, vt);
        return;
    };

    return householder_band_helper<Val>::eval_mult_left(ret, X, t_unit, *this);
};

template<class Val>
void householder_band_q<Val>::to_matrix(matcl::Matrix& ret) const
{
    return householder_band_helper<Val>::to_matrix(ret, *this);
};

template<class Val>
typename householder_band_q<Val>::data_ptr 
householder_band_q<Val>::convert(value_code new_val_code) const
{
    value_code vc   = matrix_traits::unify_value_types(new_val_code, value_code::v_float);
    Matrix v        = zeros(0,0,vc);

    data_ptr ret;
    householder_band_q_convert_vis<Val>::make<const Matrix&>(v,*this, ret);
    return ret;
};

template<class Val>
template<class T>
typename householder_band_q<Val>::data_ptr 
householder_band_q<Val>::convert_impl() const
{
    using TR        = typename details::unify_types<T, Float>::type;
    using Mat_BC    = raw::Matrix<TR,struct_banded>;
    using Mat_DC    = raw::Matrix<TR,struct_dense>;

    Mat_BC ret_refl = raw::converter<Mat_BC, Mat_B>::eval(m_reflectors);
    Mat_DC ret_tau  = raw::converter<Mat_DC, Mat_D>::eval(m_tau_vec);

    return data_ptr(new householder_band_q<TR>(m_mat_cols, ret_refl, ret_tau));
};

template<class Val>
void householder_band_q<Val>::save(oarchive& os) const
{
    //general parameters
    os << this->m_mat_cols;

    os.get() << m_reflectors;
    os.get() << m_tau_vec;
};

template<class Val>
unitary_matrix_data* householder_band_q<Val>::load(std::istream& is)
{
    using ptr_type = std::unique_ptr<householder_band_q<Val>>;
    ptr_type ret = ptr_type(new householder_band_q<Val>());

    //general parameters
    is >> ret->m_mat_cols;

    Matrix Q, tau;
    is >> Q;
    is >> tau;

    ret->m_reflectors.assign_to_fresh(Q.impl<Mat_B>());
    ret->m_tau_vec.assign_to_fresh(tau.impl<Mat_D>());

    return ret.release();
};

template<class Val>
void householder_band_q<Val>::save(std::ostream& os) const
{
    os << " ";

    //general parameters
    os << this->m_mat_cols << " ";

    os << Matrix(m_reflectors,false);
    os << Matrix(m_tau_vec,false);
};

template<class Val>
unitary_matrix_data* householder_band_q<Val>::load(iarchive& ar)
{
    using ptr_type  = std::unique_ptr<householder_band_q<Val>>;
    ptr_type ret    = ptr_type(new householder_band_q<Val>());

    //general parameters
    ar >> ret->m_mat_cols;

    ar.get() >> ret->m_reflectors;
    ar.get() >> ret->m_tau_vec;

    return ret.release();
};

//------------------------------------------------------------------------
//                          rand_unitary
//------------------------------------------------------------------------
template<class Val>
struct rand_unitary_dense
{
    using Mat = raw::Matrix<Val,struct_dense>;
    static void eval(unitary_matrix& ret, const Mat& mat)
    {
        Mat Ac(mat.get_type());

        if (mat.is_unique())
            Ac.assign_to_fresh(mat);
        else
            Ac.assign_to_fresh(mat.copy());

        Integer M   = Ac.rows();
        Integer N   = Ac.cols();
        Integer K   = std::min(M, N);

        Mat tau(mat.get_type(), K, 1);
        householder_helper<Val>::create_reflectors(Ac, tau);

        using impl_ptr = std::shared_ptr<householder_q<Val>>;

        bool isv = Ac.all_finite() && tau.all_finite();
        if (isv == false)
        {
            ret = unitary_matrix::from_nan(Ac.rows(), K, matrix_traits::value_code<Val>::value);
            return;
        };

        Integer offset = 0;
        impl_ptr um(new householder_q<Val>(K,Ac,tau,M-1,offset));

        ret = unitary_matrix(um);
        return;
    };
};


template<class Val, class Struct>
struct rand_unitary_mat_impl
{
    using Mat = raw::Matrix<Val,Struct>;
    static void eval(unitary_matrix& ret, const Mat& mat)
    {
        using VT    = typename md::unify_types<Val,Float>::type;
        using DM    = raw::Matrix<VT,struct_dense>;

        DM mat_d    = raw::converter<DM,Mat>::eval(mat);
        return rand_unitary_mat_impl<VT,struct_dense>::eval(ret,mat_d);
    };
};
template<class Val>
struct rand_unitary_mat_impl<Val,struct_dense>
{
    using Mat = raw::Matrix<Val,struct_dense>;
    static void eval(unitary_matrix& ret, const Mat& mat)
    {
        return rand_unitary_dense<Val>::eval(ret,mat);
    };
};
template<>
struct rand_unitary_mat_impl<Integer,struct_dense>
{
    using Mat = raw::Matrix<Integer,struct_dense>;
    static void eval(unitary_matrix& ret, const Mat& mat)
    {
        using VT    = Real;
        using DM    = raw::Matrix<VT,struct_dense>;

        DM mat_d    = raw::converter<DM,Mat>::eval(mat);
        return rand_unitary_mat_impl<VT,struct_dense>::eval(ret,mat_d);
    };
};

struct rand_unitary_vis : public extract_type_switch<void, rand_unitary_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, unitary_matrix& ret)
    {
        using VT = typename T::value_type;
        using ST = typename T::struct_type;
        return rand_unitary_mat_impl<VT,ST>::eval(ret,mat);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, unitary_matrix& ret)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, unitary_matrix&)
    {
        throw error::object_value_type_not_allowed("rand_unitary");
    }

    static void eval_scalar(const Matrix&, const Object&, unitary_matrix&)
    {
        throw error::object_value_type_not_allowed("rand_unitary");
    }
};

void details::rand_unitary_impl(unitary_matrix& ret, const matcl::Matrix& X)
{
    return rand_unitary_vis::make<const Matrix&>(X,ret);
};

template class householder_q<Real>;
template class householder_q<Float>;
template class householder_q<Complex>;
template class householder_q<Float_complex>;

template class householder_band_q<Real>;
template class householder_band_q<Float>;
template class householder_band_q<Complex>;
template class householder_band_q<Float_complex>;

}};

#pragma warning(pop)