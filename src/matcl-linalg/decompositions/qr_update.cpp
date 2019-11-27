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

#include "matcl-linalg/decompositions/qr.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/utils/optim_params.h"
#include "matcl-internals/func/inplace.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-linalg/decompositions/qr_band.h"
#include "matcl-linalg/decompositions/quern/matcl_quern_solver.h"
#include "matcl-linalg/decompositions/householder_q.h"
#include "matcl-linalg/decompositions/givens.h"

namespace matcl
{
    enum class qr_type
    {
        qr, rq
    };
};

namespace matcl { namespace details
{

template<class V, class S>
struct qr_update_str
{};

//--------------------------------------------------------------------------------
//                                  STRUCT TYPE
//--------------------------------------------------------------------------------
template<class V>
struct qr_update_str<V,struct_dense>
{
    using Mat = raw::Matrix<V,struct_dense>;
    static void eval_update(qr2_return& ret, const Mat& R0, const Matrix& mat_u, const Matrix& mat_w,
                            const Matrix& sigma, qr_type qt)
    {
        using VR        = typename md::real_type<V>::type;
        using Mat_R     = raw::Matrix<VR,struct_dense>;
        using Mat_I     = raw::Matrix<Integer,struct_dense>;

        Mat R           = R0.make_unique();
        Integer M       = R.rows();
        Integer N       = R.cols();
        Integer K       = mat_u.cols();

        const Mat& u    = mat_u.impl<Mat>();
        const Mat& w    = mat_w.impl<Mat>();
        const Mat& sig  = sigma.impl<Mat>();        

        // K >= 1, this should be checked before

        V* ptr_R            = R.ptr();
        Integer LDR         = R.ld();
        const V* ptr_sig    = sig.ptr();
        Integer inc_s       = sigma.length() == 1? 0 : 1;
        Integer pos_sig     = 0;

        const V* ptr_U  = u.ptr();
        const V* ptr_W  = w.ptr();
        Integer LDU     = u.ld();
        Integer LDW     = w.ld();

        Integer SLEN    = M * K + K * std::min(M,N);
        Integer SLEN_ALL= 0;
        Mat_R C         = Mat_R(ti::ti_type<VR>(), SLEN, 1);
        Mat S           = Mat(ti::ti_type<VR>(), SLEN, 1);
        Mat_I Ind       = Mat_I(ti::ti_type<VR>(), SLEN, 2);
        
        VR* ptr_C       = C.ptr();
        V* ptr_S        = S.ptr();
        Integer *ptr_I1 = Ind.ptr();
        Integer *ptr_I2 = Ind.ptr() + SLEN;

        Integer INFO;
        bool from_left  = false;
        Integer NU      = 0;

        if (qt == qr_type::qr)
        {
            //reduce rank-1 update to hessenberg
            lapack::qruphl(M, N, 0, lap(ptr_R), LDR, lap(ptr_sig[pos_sig]), lap(ptr_U), 1, lap(ptr_W), 1, 
                           lap(ptr_C), lap(ptr_S), ptr_I1, ptr_I2, SLEN, INFO);

            if (INFO != 0)
                throw error::error_general("invalid argument passed to qruphl");

            SLEN_ALL            += SLEN;

            if (K > 1)
            {
                Mat u2          = u.copy();
                LDU             = u2.ld();
                V* ptr_u2       = u2.ptr() + LDU;

                for (Integer j = 1; j < K; ++j)
                {
                    ptr_W       += LDW;

                    //apply rotations to u
                    Integer LD  = M;
                    Integer UD  = K;
                    lapack::rotseq("left", "left", "No", SLEN, lap(ptr_C), lap(ptr_S), ptr_I1, ptr_I2, M, K-j, LD, UD,
                       lap(ptr_u2), LDU, INFO);

                    if (INFO != 0)
                        throw error::error_general("invalid argument passed to rotseq");

                    pos_sig     += inc_s;
                    ptr_C       += SLEN;
                    ptr_S       += SLEN;
                    ptr_I1      += SLEN;
                    ptr_I2      += SLEN;

                    //reduce rank-1 update to hessenberg
                    lapack::qruphl(M, N, j, lap(ptr_R), LDR, lap(ptr_sig[pos_sig]), lap(ptr_u2), 1, lap(ptr_W), 1, 
                           lap(ptr_C), lap(ptr_S), ptr_I1, ptr_I2, SLEN, INFO);

                    if (INFO != 0)
                        throw error::error_general("invalid argument passed to qruphl");

                    SLEN_ALL    += SLEN;
                    ptr_u2      += LDU;
                };
            };

            // reduce hessenberg to upper triangular
            ptr_C       += SLEN;
            ptr_S       += SLEN;
            ptr_I1      += SLEN;
            ptr_I2      += SLEN;

            lapack::huundl(M, N, std::min(K, M-1), lap(ptr_R), LDR, lap(ptr_C), lap(ptr_S), ptr_I1, 
                        ptr_I2, SLEN, INFO);

            if (INFO != 0)
                throw error::error_general("invalid argument passed to huundl");

            SLEN_ALL    += SLEN;
            from_left   = true;
            NU          = M;
        }
        else if (qt == qr_type::rq)
        {
            //it should be already checked that M = N

            //reduce rank-1 update to hessenberg
            lapack::qruphr(N, 0, lap(ptr_R), LDR, lap(ptr_sig[pos_sig]), lap(ptr_U), 1, lap(ptr_W), 1, 
                           lap(ptr_C), lap(ptr_S), ptr_I1, ptr_I2, SLEN, INFO);

            if (INFO != 0)
                throw error::error_general("invalid argument passed to qruphr");

            SLEN_ALL            += SLEN;

            if (K > 1)
            {
                Mat w2          = w.copy();
                LDW             = w2.ld();
                V* ptr_w2       = w2.ptr() + LDW;

                for (Integer j = 1; j < K; ++j)
                {
                    ptr_U       += LDU;

                    //apply rotations to w
                    Integer LD  = N;
                    Integer UD  = K;
                    lapack::rotseq("right", "left", "conj tr", SLEN, lap(ptr_C), lap(ptr_S), ptr_I1, ptr_I2, N, 
                                    K-j, LD, UD, lap(ptr_w2), LDW, INFO);

                    if (INFO != 0)
                        throw error::error_general("invalid argument passed to rotseq");

                    pos_sig     += inc_s;
                    ptr_C       += SLEN;
                    ptr_S       += SLEN;
                    ptr_I1      += SLEN;
                    ptr_I2      += SLEN;

                    //reduce rank-1 update to hessenberg
                    lapack::qruphr(N, j, lap(ptr_R), LDR, lap(ptr_sig[pos_sig]), lap(ptr_U), 1, lap(ptr_w2), 1, 
                           lap(ptr_C), lap(ptr_S), ptr_I1, ptr_I2, SLEN, INFO);

                    if (INFO != 0)
                        throw error::error_general("invalid argument passed to qruphl");

                    SLEN_ALL    += SLEN;
                    ptr_w2      += LDW;
                };
            };

            // reduce hessenberg to upper triangular
            ptr_C       += SLEN;
            ptr_S       += SLEN;
            ptr_I1      += SLEN;
            ptr_I2      += SLEN;

            lapack::huundr(N, std::min(K, N-1), 0, lap(ptr_R), LDR, lap(ptr_C), lap(ptr_S), ptr_I1, 
                        ptr_I2, SLEN, INFO);

            if (INFO != 0)
                throw error::error_general("invalid argument passed to huundr");

            SLEN_ALL    += SLEN;
            from_left   = false;
            NU          = N;
        };

        C.assign_to_fresh(C.resize(SLEN_ALL,1));
        S.assign_to_fresh(S.resize(SLEN_ALL,1));
        Ind.assign_to_fresh(Ind.resize(SLEN_ALL,2));

        unitary_matrix Q = givens_to_unitary(NU, Matrix(C,false), Matrix(S,false), 
                                             Matrix(Ind,false), from_left);
        Q                = ctrans(Q);

        R.get_struct().reset();
        R.get_struct().add(predefined_struct_type::triu);

        Matrix RR = Matrix(R,true);
        ret = qr2_return(Q,RR);

        return;
    };
};

template<class V>
struct qr_update_str<V,struct_sparse>
{
    using Mat   = raw::Matrix<V,struct_sparse>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval_update(qr2_return& ret, const Mat& R, const Matrix& u, const Matrix& w, const Matrix& sigma,
                            qr_type qt)
    {
        Mat_D Ac    = raw::converter<Mat_D, Mat>::eval(R);
        return qr_update_str<V,struct_dense>::eval_update(ret, Ac, u, w, sigma, qt);
    };
};

template<class V>
struct qr_update_str<V,struct_banded>
{
    using Mat   = raw::Matrix<V,struct_banded>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval_update(qr2_return& ret, const Mat& R, const Matrix& u, const Matrix& w, const Matrix& sigma,
                            qr_type qt)
    {
        Mat_D Ac    = raw::converter<Mat_D, Mat>::eval(R);
        return qr_update_str<V,struct_dense>::eval_update(ret, Ac, u, w, sigma, qt);
    };
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct qr_update
{
    using M = raw::Matrix<V,S>;

    static void eval_update(matcl::qr2_return& ret, const M& R, const Matrix& u, const Matrix& w, 
                            const Matrix& sigma, qr_type qt)
    {
        return qr_update_str<V,S>::eval_update(ret, R, u, w, sigma, qt);
    };
};

template<class S>
struct qr_update<Integer,S>
{
    using M = raw::Matrix<Integer,S>;

    static void eval_update(matcl::qr2_return& ret, const M& A, const Matrix& u, const Matrix& w, 
                            const Matrix& sigma, qr_type qt)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        return qr_update<Real,S>::eval_update(ret, AC, u, w, sigma, qt);
    };
};

struct visitor_qr_update : public extract_type_switch<void,visitor_qr_update,true>
{
    template<class T>
    static void eval(const Matrix& h, const T& mat, qr2_return& ret, const Matrix& u, const Matrix& v, 
                     const Matrix& sigma, qr_type qt)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        value_code vt1  = matrix_traits::value_code<V>::value;
        value_code vt2  = u.get_value_code();
        value_code vt3  = v.get_value_code();
        value_code vt4  = sigma.get_value_code();
        value_code vt0  = matrix_traits::unify_value_types(vt1, vt2);
        vt0             = matrix_traits::unify_value_types(vt0, vt3);
        vt0             = matrix_traits::unify_value_types(vt0, vt4);
        value_code vt   = matrix_traits::unify_value_types(vt0, value_code::v_float);

        if (vt != vt1)
        {
            //conversion is required
            struct_code st  = h.get_struct_code();
            mat_code mc     = matrix_traits::get_matrix_type(vt, st);
            Matrix Ac       = convert(h, mc);
            return visitor_qr_update::make<const Matrix&>(Ac,ret,u,v,sigma,qt);
        };

        return details::qr_update<V,S>::eval_update(ret, mat, u, v, sigma, qt);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix&, const T& A, qr2_return& ret, const Matrix& u, const Matrix& v, 
                            const Matrix& Sigma, qr_type)
    {         
        Matrix RR   = A + Sigma * u * ctrans(v);
        Matrix Q    = speye(1, 1, RR.get_value_code());
        ret         = qr2_return(unitary_matrix(Q,false),RR);
        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, qr2_return&, 
                     const Matrix&, const Matrix&, const Matrix&, qr_type )
    {
        throw error::object_value_type_not_allowed("qr_update");
    };
    static void eval_scalar(const Matrix&, const Object&, qr2_return&, const Matrix&, const Matrix&, 
                            const Matrix&, qr_type )
    {
        throw error::object_value_type_not_allowed("qr_update");
    };
};

static void qr_update_impl(qr2_return& ret, const Matrix& R, const Matrix& u, const Matrix& v, 
                             const Matrix& sigma, qr_type qt)
{
    if (is_triu(R) == false)
        throw error::error_invalid_qr_factor();

    if (qt == qr_type::rq && R.rows() != R.cols())
        throw error::error_invalid_rq_factor(R.rows(), R.cols());

    if (u.rows() != R.rows() || v.rows() != R.cols())
        throw error::error_invalid_qr_update(R.rows(), R.cols(), u.rows(), u.cols(), v.rows(), v.cols());

    Integer K_vec   = u.cols();
    Integer K_r     = sigma.rows();
    Integer K_c     = sigma.cols();
    bool sigma_ok   = (K_r == 1 && (K_c == 1 || K_c == K_vec) )
                    || (K_c == 1 && (K_r == 1 || K_r == K_vec) );

    if (sigma_ok == false)
        throw error::error_invalid_update_sigma(K_vec, K_r, K_c);

    if (sigma.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("qr_update");

    if (K_vec == 0)
    {
        Integer NU  = qt == qr_type::qr ? R.rows() : R.cols();
        Matrix Q    = speye(NU, NU, R.get_value_code()); 
        ret         = qr2_return(unitary_matrix(Q,false), R);
        return;
    };

    return details::visitor_qr_update::make<const Matrix&>(R,ret,u,v,sigma,qt);
};

}};

namespace matcl
{

qr2_return matcl::qr_update(const Matrix& R0, const Matrix& U, const Matrix& W, 
                   const Matrix& sigma)
{
    Matrix R(R0);

    qr2_return ret;
    details::qr_update_impl(ret, R, U, W, sigma, qr_type::qr);
    return ret;
};
qr2_return matcl::qr_update(Matrix&& R0, const Matrix& U, const Matrix& W, 
                   const Matrix& sigma)
{
    Matrix R(std::move(R0));

    qr2_return ret;
    details::qr_update_impl(ret, R, U, W, sigma, qr_type::rq);
    return ret;
};

qr2_return matcl::rq_update(const Matrix& R0, const Matrix& U, const Matrix& W, 
                   const Matrix& sigma)
{
    Matrix R(R0);

    qr2_return ret;
    details::qr_update_impl(ret, R, U, W, sigma, qr_type::rq);
    return ret;
};
qr2_return matcl::rq_update(Matrix&& R0, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma)
{
    Matrix R(std::move(R0));

    qr2_return ret;
    details::qr_update_impl(ret, R, U, W, sigma, qr_type::rq);
    return ret;
};

qr2_return matcl::qr_update(const unitary_matrix& Q, const Matrix& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma)
{
    //Q * R + U * Sigma * W = Q * (R + Q' * U * Sigma * W)
    Matrix U2   = mmul(Q, U, trans_type::conj_trans);

    Matrix R2;
    unitary_matrix Q2;
    tie(Q2, R2)         = qr_update(R, U2, W, sigma);
    unitary_matrix Qret = Q * Q2;

    return qr2_return(Qret, R2);
};
qr2_return matcl::qr_update(const unitary_matrix& Q, Matrix&& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma)
{
    if (Q.is_square() == false)
        throw error::invalid_size(Q.rows(), Q.cols());

    //Q * R + U * Sigma * W = Q * (R + Q' * U * Sigma * W)
    Matrix U2   = mmul(Q, U, trans_type::conj_trans);

    Matrix R2;
    unitary_matrix Q2;
    tie(Q2, R2)         = qr_update(std::move(R), U2, W, sigma);
    unitary_matrix Qret = Q * Q2;

    return qr2_return(Qret, R2);
};
qr2_return matcl::qr_update(const Matrix& Q, const Matrix& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma)
{
    if (Q.is_square() == false)
        throw error::invalid_size(Q.rows(), Q.cols());

    //Q * R + U * Sigma * W = Q * (R + Q' * U * Sigma * W)
    Matrix U2   = mmul(Q, U, trans_type::conj_trans);

    Matrix R2;
    unitary_matrix Q2;
    tie(Q2, R2)         = qr_update(R, U2, W, sigma);
    unitary_matrix Qret = unitary_matrix(Q * Q2, true);

    return qr2_return(Qret, R2);
};
qr2_return matcl::qr_update(const Matrix& Q, Matrix&& R, const Matrix& U, const Matrix& W, 
                                       const Matrix& sigma)
{
    if (Q.is_square() == false)
        throw error::invalid_size(Q.rows(), Q.cols());

    //Q * R + U * Sigma * W = Q * (R + Q' * U * Sigma * W)
    Matrix U2   = mmul(Q, U, trans_type::conj_trans);

    Matrix R2;
    unitary_matrix Q2;
    tie(Q2, R2)         = qr_update(std::move(R), U2, W, sigma);
    unitary_matrix Qret = unitary_matrix(Q * Q2, true);

    return qr2_return(Qret, R2);
};

qr2_return matcl::rq_update(const Matrix& Q, const Matrix& R, const Matrix& U, const Matrix& W, 
                   const Matrix& sigma)
{
    //R*Q + U * Sigma * W' = (R + U * Sigma * (Q*W)')*Q

    if (Q.is_square() == false)
        throw error::invalid_size(Q.rows(), Q.cols());

    Matrix W2   = mmul(Q, W);

    Matrix R2;
    unitary_matrix Q2;
    tie(Q2, R2)         = rq_update(R, U, W2, sigma);
    unitary_matrix Qret = unitary_matrix(Q2*Q, true);

    return qr2_return(Qret, R2);
};
qr2_return matcl::rq_update(const Matrix& Q, Matrix&& R, const Matrix& U, const Matrix& W, 
                   const Matrix& sigma)
{
    //R*Q + U * Sigma * W' = (R + U * Sigma * (Q*W)')*Q

    if (Q.is_square() == false)
        throw error::invalid_size(Q.rows(), Q.cols());

    Matrix W2   = mmul(Q, W);

    Matrix R2;
    unitary_matrix Q2;
    tie(Q2, R2)         = rq_update(std::move(R), U, W2, sigma);
    unitary_matrix Qret = unitary_matrix(Q2*Q,true);

    return qr2_return(Qret, R2);
};
qr2_return matcl::rq_update(const unitary_matrix& Q, const Matrix& R, const Matrix& U, const Matrix& W, 
                   const Matrix& sigma)
{
    //R*Q + U * Sigma * W' = (R + U * Sigma * (Q*W)')*Q

    if (Q.is_square() == false)
        throw error::invalid_size(Q.rows(), Q.cols());

    Matrix W2   = mmul(Q, W);

    Matrix R2;
    unitary_matrix Q2;
    tie(Q2, R2)         = rq_update(R, U, W2, sigma);
    unitary_matrix Qret = unitary_matrix(Q2*Q);

    return qr2_return(Qret, R2);
};
qr2_return matcl::rq_update(const unitary_matrix& Q, Matrix&& R, const Matrix& U, const Matrix& W, 
                   const Matrix& sigma)
{
    //R*Q + U * Sigma * W' = (R + U * Sigma * (Q*W)')*Q

    if (Q.is_square() == false)
        throw error::invalid_size(Q.rows(), Q.cols());

    Matrix W2   = mmul(Q, W);

    Matrix R2;
    unitary_matrix Q2;
    tie(Q2, R2)         = rq_update(std::move(R), U, W2, sigma);
    unitary_matrix Qret = unitary_matrix(Q2*Q);

    return qr2_return(Qret, R2);
};

};
