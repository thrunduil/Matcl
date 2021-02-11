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

#include "arnoldi_process.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-ext/lapack_ext/blas_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-core/details/complex_details.h"
#include "matcl-blas-lapack/level1/level1.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-linalg/utils/optim_params.h"
#include <algorithm>

namespace matcl
{

namespace md = matcl::details;

//---------------------------------------------------------------------------
//                      UTILS
//---------------------------------------------------------------------------
inline Real real(Real a)
{
    return a;
};
inline Float real(Float a)
{
    return a;
};

inline Real abs(const Real& v)
{
    return std::abs(v);
};
inline Float abs(const Float& v)
{
    return std::abs(v);
};

inline Real abs(const Complex& v)
{
    return std::abs(v.value);
}
inline Float abs(const Float_complex& v)
{
    return std::abs(v.value);
}

inline Complex operator+(const Complex& a, const Complex& b)
{
    return md::plus_c(a,b);
};
inline Float_complex operator+(const Float_complex& a, const Float_complex& b)
{
    return md::plus_c(a,b);
};
inline Complex operator-(const Complex& a, const Complex& b)
{
    return md::minus_c(a,b);
};
inline Float_complex operator-(const Float_complex& a, const Float_complex& b)
{
    return md::minus_c(a,b);
};

inline Complex operator-(const Complex& a)
{
    return md::uminus_c(a);
};
inline Float_complex operator-(const Float_complex& a)
{
    return md::uminus_c(a);
};

};

namespace matcl { namespace raw
{

//---------------------------------------------------------------------------
//                          aitr
//---------------------------------------------------------------------------
template<class V>
struct callback_gsg
{    
    using VL            = typename md::lusol_value_type<V>::type;

    callback_aitr<V>*   m_callback;

    callback_gsg(callback_aitr<V>* callback)
        :m_callback(callback)
    {};

    void eval(Integer N, VL* X, VL* WORK)
    {
        m_callback->eval(2, N, (V*)X, (V*)WORK, (V*)nullptr);
    };

    static void eval(void* CTX, Integer N, VL* X, VL* WORK)
    {
        callback_gsg* call = reinterpret_cast<callback_gsg*>(CTX);
        call->eval(N, X, WORK);
    };
};

template<class V>
void raw::aitr(bool SYM, bool GEN, Integer N, Integer K, Integer NP, 
        typename md::real_type<V>::type TOL, typename md::real_type<V>::type& RNORM, 
        V* U, Integer LDU, V* H_NS, typename md::real_type<V>::type* H_S, Integer LDH, 
        typename md::real_type<V>::type& SCALE, typename md::real_type<V>::type& SUMSQ,
        V* WORK, callback_aitr<V>* callback, Integer& INFO)
{
    using VR = typename md::real_type<V>::type;

    using md::lap;

    //check arguments
    INFO                = 0;
    if (N < 0)
        INFO            = -3;
    else if (K < 0 || K > N)
        INFO            = -4;
    else if (NP < 0 || K + NP > N)
        INFO            = -5;
    else if (LDU < N)
        INFO            = -9;
    else if (LDH < K + NP)
        INFO            = -12;

    if (INFO < 0)
        return;

    //fast exit
    if (N == 0)
        return;

    // set machine constants
    VR SAFEMIN          = lapack::lamch<VR>("safmin");
    VR ULP              = lapack::lamch<VR>("PRECISION");
    TOL                 = std::max(TOL, ULP);

    INFO                = 0;
    Integer POS_RES     = 0;                //store resid here
    Integer POS_B       = POS_RES + N;      //store B*r here
    Integer POS_2       = POS_B + N;        // work
    Integer POS_H       = POS_2 + N;        //work for H_S

    callback_gsg<V> call_gsg(callback);

    //compute norm of the initial vector
    if (GEN)
    {
        callback->eval(2, N, WORK + POS_RES, WORK + POS_B, nullptr);

        RNORM           = abs(lapack::dotc(N, lap(WORK + POS_RES), 1, lap(WORK + POS_B), 1));
        RNORM           = sqrt(RNORM);
    }
    else
    {
        RNORM           = lapack::nrm2(N, lap(WORK + POS_RES), 1);
    };

    if (K == 0)
    {
        SCALE    = VR(0.0);
        SUMSQ    = VR(1.0);
    };

    //-------------------------------------------------------------------
    //       A R N O L D I     I T E R A T I O N     L O O P
    //-------------------------------------------------------------------
    for (Integer J = K + 1; J <= K + NP; ++J)
    {
        VR BETAJ        = RNORM;
        VR NORM_F       = SCALE * sqrt(SUMSQ); //lower bound on the second norm
        NORM_F          = NORM_F / sqrt(VR(J));

        if (RNORM == VR(0.0))
        {
            // Invariant subspace found
            INFO        = J - 1;
            return;
        };

        //deflation check
        if (RNORM < TOL * NORM_F)
        {
            //torelance reached
            INFO        = J - 1;
            return;
        };
 
        // STEP 2:  v_{J} = r_{J-1}/RNORM and p_{J} = p_{J}/RNORM
        // Note that p_{J} = B*r_{J-1}. In order to avoid overflow when reciprocating 
        // a small RNORM, test against lower machine bound.

        if ( RNORM >= SAFEMIN)
        {
            V TEMP1     = VR(1.0) / RNORM;

            //put scaled resid in WORK(POS_2)
            lapack::yax(N, lap(TEMP1), lap(WORK + POS_RES), 1, lap(WORK + POS_2), 1);

            if (GEN)
                lapack::scal(N, lap(TEMP1), lap(WORK + POS_B), 1);
        }
        else
        {
            Integer DUM = 0;
            Integer INFOL;

            lapack::copy(N, lap(WORK + POS_RES), 1, lap(WORK + POS_2), 1);

            // To scale both v_{J} and p_{J} carefully use LAPACK routine AR_SLASCL
            lapack::lascl("General", DUM, DUM, RNORM, VR(1.0), N, 1, lap(WORK + POS_2), N, INFOL);
            
            if (GEN)
                lapack::lascl("General", DUM, DUM, RNORM, VR(1.0), N, 1, lap(WORK + POS_B), N, INFOL);
        };

        lapack::copy(N, lap(WORK + POS_2), 1, lap(U + (J-1)*LDU), 1);

        // STEP 3:  r_{J} = OP*v_{J}; Note that p_{J} = B*v_{J}
        // Note that this is not quite yet r_{J}. See STEP 4 

        callback->eval(1, N, WORK + POS_2, WORK + POS_RES, WORK + POS_B);

        // STEP 4:  Finish extending the Arnoldi to length J.        
        // orthogonalize OP * r against U and update H

        if (SYM == false)
        {            
            Integer INFO2;
            if (GEN)
            {
                lapack::gsg_ort(callback_gsg<V>::eval, &call_gsg, "N", N, J, lap(U), LDU, 
                    lap(WORK + POS_RES), 1, lap(H_NS + (J-1)*LDH), 1, RNORM, lap(WORK + POS_B), 
                    lap(WORK + POS_2), INFO2);
            }
            else
            {
                lapack::gs_ort("N", N, J, lap(U), LDU, lap(WORK + POS_RES), 1, 
                            lap(H_NS + (J-1)*LDH), 1, RNORM, lap(WORK + POS_2), INFO2);
            };

            //update the Frobenius norm
            lapack::lassq(J, lap(H_NS + (J-1)*LDH), 1, SCALE, SUMSQ);

            if (J > 1) 
            {
                H_NS[J-1 + (J-1-1)*LDH] = BETAJ;
                lapack::lassq(1, lap(&BETAJ), 1, SCALE, SUMSQ);
            }
        }
        else
        {
            Integer INFO2;
            if (GEN)
            {
                lapack::gsg_ort(callback_gsg<V>::eval, &call_gsg, "N", N, J, lap(U), LDU, 
                    lap(WORK + POS_RES), 1, lap(WORK + POS_H), 1, RNORM, lap(WORK + POS_B), 
                    lap(WORK + POS_2), INFO2);
            }
            else
            {
                lapack::gs_ort("N", N, J, lap(U), LDU, lap(WORK + POS_RES), 1, 
                            lap(WORK + POS_H), 1, RNORM, lap(WORK + POS_2), INFO2);
            };

            H_S[J-1 + LDH]  = real(WORK[POS_H + J - 1]);

            //update the Frobenius norm
            lapack::lassq(1, lap(&H_S[J-1 + LDH]), 1, SCALE, SUMSQ);

            if (J == 1)
                H_S[0]      = VR(0.0);
            else
            {
                H_S[J-1]    = BETAJ;
                lapack::lassq(1, lap(&BETAJ), 1, SCALE, SUMSQ);
                lapack::lassq(1, lap(&BETAJ), 1, SCALE, SUMSQ);
            }
        };

        if (SYM == true)
        {
            // Make sure the last off-diagonal element is non negative
            // If not perform a similarity transformation on H(1:j,1:j) 
            // and scale U(:,j) by -1.                                  |
            if (H_S[J-1] < VR(0.0))
            {
                H_S[J-1]    = -H_S[J-1];

                if (J < K + NP)
                    lapack::scal(N, lap(-V(1.0)), lap(U + J*LDU), 1);
                else
                    lapack::scal(N, lap(-V(1.0)), lap(WORK + POS_RES), 1);
            };
        };
    };

    return;
};

template
void aitr<Real>(bool SYM, bool BMAT, Integer N, Integer K, Integer NP, 
        Real TOL, Real& RNORM, Real* U, Integer LDU, Real* H_NS, Real* H_S, 
        Integer LDH, Real& SCALE, Real& SUMSQ, Real* WORK, callback_aitr<Real>* callback, 
        Integer& INFO);
template
void aitr<Float>(bool SYM, bool BMAT, Integer N, Integer K, Integer NP, 
        Float TOL, Float& RNORM, Float* U, Integer LDU, Float* H_NS, Float* H_S, 
        Integer LDH, Float& SCALE, Float& SUMSQ, Float* WORK, callback_aitr<Float>* callback, 
        Integer& INFO);

template
void aitr<Complex>(bool SYM, bool BMAT, Integer N, Integer K, Integer NP, 
        Real TOL, Real& RNORM, Complex* U, Integer LDU, Complex* H_NS, Real* H_S, 
        Integer LDH, Real& SCALE, Real& SUMSQ, Complex* WORK, callback_aitr<Complex>* callback,
        Integer& INFO);

template
void aitr<Float_complex>(bool SYM, bool BMAT, Integer N, Integer K, Integer NP, 
        Float TOL, Float& RNORM, Float_complex* U, Integer LDU, Float_complex* H_NS, Float* H_S, 
        Integer LDH, Float& SCALE, Float& SUMSQ, Float_complex* WORK, 
        callback_aitr<Float_complex>* callback, Integer& INFO);

//---------------------------------------------------------------------------
//                          baitr
//---------------------------------------------------------------------------
template<class V>
void orm2r(const char* TRANS, Integer M, Integer N, Integer K, V* A, Integer LDA, const V* TAU, V* C, 
           Integer LDC, V* WORK)
{
    using md::lap;

    bool NOTRAN = TRANS[0] == 'N' || TRANS[0] == 'n';

    // Quick return if possible
    if (M == 0 || N == 0 || K == 0 )
        return;

    if (NOTRAN == false)
    {
        for (Integer I = 0; I < K; I += 1)
        {
            // Apply H(i)
            V AII           = A[I + I*LDA];
            A[I + I*LDA]    = V(1.0);
        
            lapack::larf("L", M - I, N, lap(A + I + I*LDA), 1, lap(TAU + I), lap(C + I), LDC, lap(WORK) );
            A[I + I*LDA]    = AII;
        };
    }
    else
    {
        for (Integer I = K-1; I >= 0; I -= 1)
        {
            // Apply H(i)
            V AII           = A[I + I*LDA];
            A[I + I*LDA]    = V(1.0);
        
            lapack::larf("L", M - I, N, lap(A + I + I*LDA), 1, lap(TAU + I), lap(C + I), LDC, lap(WORK) );
            A[I + I*LDA]    = AII;
        };
    };

    return;
};

template<class V>
void raw::apply_U(const char* TRANS, Integer M, Integer N, Integer K, V* U, Integer LDU, const V* TAU,
        V* X, Integer LDX, V* T, Integer LDT, Integer& NT, Integer NTB, V* WORK, Integer LDWORK0)
{
    using VR    = typename md::real_type<V>::type;
    using md::lap;

    bool NOTRAN     = TRANS[0] == 'N' || TRANS[0] == 'n';
    static const Integer NBMAX   = 64;

    // Determine the block size.
    Integer NB      = std::min( NBMAX, NTB);

    // Quick return if possible
    if ( M == 0 || N == 0 || K == 0 )
        return;

    Integer NBMIN   = 2;    
    Integer NBX     = LDWORK0/NB;    

    if ( NB < NBMIN || NB > K || NBX < NBMIN)
    {
        // Use unblocked code
        orm2r<V>(TRANS, M, N, K, U, LDU, TAU, X, LDX, WORK);
    }
    else
    {
        //accumulate next triangular factor
        Integer NT2 = K/NB;

        for (Integer J = NT; J < NT2; ++J)
        {      
            Integer I   = J * NB;

            // Form the triangular factor of the block reflector
            // H = H(i) H(i+1) . . . H(i+ib-1)

            lapack::larft("F", "C", M-I, NB, lap(U + I + I*LDU), LDU, lap(TAU + I), lap(T + J*LDT), NB );
            NT  = NT + 1;
        };

        Integer K1  = K / NB;
        Integer K2  = K1 * NB;        

        // Use blocked code
        if (NOTRAN == false)
        {
            for (Integer P = 0; P < N; P += NBX)
            {
                Integer JNBX    = std::min(NBX, N - P);
                Integer LDWORKB = JNBX;

                for (Integer I = 0, J = 0; I < K2; I += NB, J += 1)
                {
                    // Apply H or H'
                    lapack::larfb("L", TRANS, "F", "C", M - I, JNBX, NB, lap(U + I + I*LDU), LDU, lap(T + J*LDT), NB, 
                                  lap(X + I + P*LDX), LDX, lap(WORK), LDWORKB );
                };
            };

            if (K > K2)
            {
                //last block
                orm2r<V>(TRANS, M - K2, N, K - K2, U + K2 + K2 * LDU, LDU, TAU + K2, X + K2, LDX, WORK);
            };
        }
        else
        {            
            //first block
            if (K > K2)
                orm2r<V>(TRANS, M - K2, N, K - K2, U + K2 + K2 * LDU, LDU, TAU + K2, X + K2, LDX, WORK);

            for (Integer P = 0; P < N; P += NBX)
            {
                Integer JNBX    = std::min(NBX, N - P);
                Integer LDWORKB = JNBX;

                for (Integer I = K2 - NB, J = K1 - 1; I >= 0; I -= NB, J -= 1)
                {
                    // Apply H or H'
                    lapack::larfb("L", TRANS, "F", "C", M - I, JNBX, NB, lap(U + I + I*LDU), LDU, lap(T + J*LDT), NB, 
                                  lap(X + I + P*LDX), LDX, lap(WORK), LDWORKB);
                };
            };
        };
    };

    return;
};

template<class V>
Integer baitr_block_size()
{
    return md::linalg_optim_params::block_size_ORMQR();
};

template<class V>
void raw::baitr(bool SYM, Integer N, Integer K, Integer NP, Integer& KB, typename md::real_type<V>::type TOL, 
    V* X, Integer LDX, bool INIT, V* RESID, Integer LDR, typename md::real_type<V>::type& NORMR, V* U, Integer LDU,
    V* TAU, V* H_NS, V* H_S, Integer LDH, typename md::real_type<V>::type& SCALE, typename md::real_type<V>::type& SUMSQ,
    V* T, Integer LDT, Integer& NT, V* WORK, Integer LWORK, callback_baitr<V>* callback, Integer& INFO)
{
    using VR    = typename md::real_type<V>::type;
    using md::lap;

    Integer NTB         = baitr_block_size<V>();

    //check arguments
    INFO                = 0;
    if (N < 0)
        INFO            = -2;
    else if (K < 0 || K > N)
        INFO            = -3;
    else if (NP < 0 || K + NP > N)
        INFO            = -4;
    else if (KB < 1)
        INFO            = -5;
    else if (LDX < N)
        INFO            = -8;
    else if (LDR < N)
        INFO            = -11;
    else if (LDU < N)
        INFO            = -14;
    else if (LDH < K + NP && SYM == false)
        INFO            = -18;
    else if (LDH < std::min(2*KB, K + NP) && SYM == true)
        INFO            = -18;
    else if (LDT < NTB * NTB)
        INFO            = -22;
    else if (NT < 0)
        INFO            = -23;

    if (INFO < 0)
        return;

    Integer ierr;

    if (K == 0 || INIT == true)
        KB              = std::min(KB, NP);

    Integer LWRK_PIV    = KB * 2;

    //query LWKR_QR
    lapack::geqp3(N, KB, lap(X), LDX, (Integer*)nullptr, lap(TAU), lap(WORK), -1, &ierr);
    Integer LWKR_QR     = (Integer)real(WORK[0]);
    
    Integer KBO         = std::min(KB, 2 * NTB);
    Integer LWRK_ORMQR  = NTB * KBO;

    Integer LWORK_QR    = std::max(LWKR_QR, LWRK_ORMQR);
    Integer LWRK        = LWRK_PIV + LWORK_QR;

    if (LWORK == -1)
    {
        WORK[0]         = V(VR(LWRK));
        return;
    };

    if (LWORK < LWRK)
    {
        INFO            = -25;
        return;
    };

    NORMR               = VR(0.0);

    //fast exit
    if (N == 0)
    {
        WORK[0]         = V(VR(LWRK));
        return;
    };       

    Integer POS_PIV     = 0;           //size 2 * KB;
    Integer POS_QR      = POS_PIV + 2*KB;   //size LWORK_QR

    Integer* JPVT       = reinterpret_cast<Integer*>(WORK + POS_PIV);
    Integer* WORK_JPVT  = JPVT + KB;
    V* WORK_QR          = WORK + POS_QR;    

    bool is_complex         = md::is_complex<V>::value;
    const char* trans_char  = is_complex ? "C" : "T";

    VR ULP              = lapack::lamch<VR>("PRECISION");    

    if (K > 0 && INIT == true)
    {
        // form s = V' * X
        apply_U<V>(trans_char, N, KB, K, U, LDU, TAU, X, LDX, T, LDT, NT, NTB, WORK_QR, LWORK_QR);
    };

    if (K == 0)
    {
        SCALE    = VR(0.0);
        SUMSQ    = VR(1.0);
    };

    // begin of main iteration loop for the Augmented Block Householder Arnoldi
    // decomposition.
    Integer JKB         = KB;
    Integer J           = K;

    for (; J < K + NP; )
    {
        VR NORM_F       = SCALE * sqrt(SUMSQ); //lower bound on the second norm
        NORM_F          = NORM_F / sqrt(VR(J+1));
        Real JTOL       = ULP * NORM_F;

        //QR factorization with pivoting
        level1::set_val<Integer,0>::eval(JPVT, JKB, 0);
        lapack::geqp3(N-J, JKB, lap(X + J), LDX, JPVT, lap(TAU + J), lap(WORK_QR), LWORK_QR, &ierr);

        //permutations to col interchanges
        lapack::perm2int(JKB, JPVT, WORK_JPVT, false);

        //deflate 
        Integer NDEFL   = 0;
        for (Integer k = std::min(N-J,JKB) - 1; k >= 0; --k)
        {
            V diag      = X[J + k + k * LDX];

            if (abs(diag) < JTOL)
                ++NDEFL;
            else
                break;
        };

        Integer JKB_old = JKB;
        JKB             = JKB - NDEFL;        
        JKB             = std::min(JKB, K + NP - J);

        if (JKB == 0)
            break;

        //update H
        if (J > 0)
        {
            if (SYM == false)
            {
                //H(J+1:J+JKB,J+1:J+JKB)    = triu(X(1:J, 1:J))
                lapack::lacpy("U", JKB, JKB_old, lap(X + J), LDX, lap(H_NS + J + (J-JKB_old)*LDH), LDH);
                lapack::laset("L", JKB - 1, JKB_old - 1, lap(V(0.0)), lap(V(0.0)), 
                              lap(H_NS + J + 1 + (J-JKB_old)*LDH), LDH);

                //update the Frobenius norm
                for (Integer k = 0; k < JKB_old; ++k)
                    lapack::lassq(JKB, lap(H_NS + J + (J-JKB_old+k)*LDH), 1, SCALE, SUMSQ);

                //apply permutations                
                lapack::laswpc(JKB, lap(H_NS + J + (J-JKB_old)*LDH), LDH, 1, JKB_old, JPVT, -1);
            }
            else
            {
                //H(J+1:J+JKB,J+1:J+JKB)    = triu(X(1:J, 1:J))
                lapack::lacpy("U", JKB, JKB_old, lap(X + J), LDX, lap(H_S + JKB_old + (J-JKB_old)*LDH), LDH);
                lapack::laset("L", JKB - 1, JKB_old - 1, lap(V(0.0)), lap(V(0.0)), 
                              lap(H_S + JKB_old + 1 + (J-JKB_old)*LDH), LDH);

                //update the Frobenius norm
                for (Integer k = 0; k < JKB_old; ++k)
                {
                    //two times the same vector
                    lapack::lassq(JKB, lap(H_S + JKB_old + (J-JKB_old+k)*LDH), 1, SCALE, SUMSQ);
                    lapack::lassq(JKB, lap(H_S + JKB_old + (J-JKB_old+k)*LDH), 1, SCALE, SUMSQ);
                }

                //apply permutations
                lapack::laswpc(JKB, lap(H_S + JKB_old + (J-JKB_old)*LDH), LDH, 1, JKB_old, JPVT, -1);

                //extract diagonals
                for (Integer i = 1; i < JKB_old; ++i)
                {
                    level1::copy<V,V,0>::eval(H_S + JKB_old - i + (J-JKB_old+i)*LDH, 1,
                                               H_S + JKB_old + (J-JKB_old+i)*LDH, 1, JKB);
                    level1::set_val<V,0>::eval(H_S + JKB_old + JKB - i + (J-JKB_old+i)*LDH, 1, i, V(0.0));
                };
            };
        };

        // copy residual vector to U
        level1::copy_mat<true,V,V,0,0,0>::eval(U + J + J * LDU, LDU, X + J, LDX, N-J, JKB);

        // prepare ( I + V_m * T_m * V_m') * E_m -> orthogonal residual vector.
        level1::set_val_mat<true,V,0,0,0>::eval(RESID, LDR, N, JKB, V(0.0));

        for (Integer i = 0; i < JKB; ++i)
            RESID[J + i + i *LDR]   = V(1.0);

        apply_U("N", N, JKB, J + JKB, U, LDU, TAU, RESID, LDR, T, LDT, NT, NTB, WORK_QR, LWORK_QR);

        // form r2 = OP * r
        callback->eval(N, JKB, RESID, LDR, X, LDX);

        // form s = V' * r2
        apply_U(trans_char, N, JKB, J + JKB, U, LDU, TAU, X, LDX, T, LDT, NT, NTB, WORK_QR, LWORK_QR);
         
        // Compute the columns of the upper block hessenberg matrix H.
        if (SYM == false)
        {            
            //H(1:J+JKB, J+1:J+JKB) = X(1:J+JKB,J+1:J+JKB);
            level1::copy_mat<true,V,V,0,0,0>::eval(H_NS + J*LDH, LDH, X, LDX, J+JKB, JKB);

            //update the Frobenius norm
            for (Integer k = 0; k < JKB; ++k)
                lapack::lassq(J+JKB, lap(H_NS + (J+k)*LDH), 1, SCALE, SUMSQ);
        }
        else
        {
            for (Integer i = 0; i < JKB; ++i)
            {
                level1::copy<V,V,0>::eval(H_S + (J+i)*LDH, 1, X + J + i + i * LDX, 1, JKB - i);
                level1::set_val<V,0>::eval(H_S + JKB - i + (J+i)*LDH, 1, i, V(0.0));

                //update the Frobenius norm
                lapack::lassq(JKB, lap(X + J + i * LDX), 1, SCALE, SUMSQ);
            };
        };

        //evaluate the Frobenius norm of residuals
        VR R_SCALE      = VR(0.0);
        VR R_SUMSQ      = VR(1.0);

        for (Integer k = 0; k < JKB; ++k)
            lapack::lassq(N - J - JKB, lap(X + J + JKB + k*LDX), 1, R_SCALE, R_SUMSQ);

        J               = J + JKB;

        NORMR           = R_SCALE * sqrt(R_SUMSQ);

        if (NORMR < TOL * NORM_F)
            break;

        if (J + JKB > K + NP && K + NP < N)
        {
            //not enough place to store next Arnoldi vectors
            break;
        };
    };

    // prepare residual
    if (J < N)
    {
        // X stores U'*RESID,
        level1::set_val_mat<true, V, 0,0,0>::eval(RESID, LDR, J, JKB, V(0.0));
        level1::copy_mat<true,V,V,0,0,0>::eval(RESID + J, LDR, X + J, LDX, N - J, JKB);

        //recreate RESID
        apply_U("N", N, JKB, J, U, LDU, TAU, RESID, LDR, T, LDT, NT, NTB, WORK_QR, LWORK_QR); 

        KB = JKB;
    }
    else
    {
        KB      = 0;
        NORMR   = VR(0.0);
    }

    if (J < K + NP)
        INFO    = J;

    WORK[0]         = V(VR(LWRK));
};

template
void baitr<Real>(bool SYM, Integer N, Integer K, Integer NP, Integer& KB,
        Real TOL, Real* X, Integer LDX, bool INIT, Real* RESID, Integer LDR, Real& NORMR, Real* U, Integer LDU, 
        Real* TAU, Real* H_NS, Real* H_S, Integer LDH, Real& SCALE, Real& SUMSQ, Real* T, Integer LDT, Integer& NT, 
        Real* WORK, Integer LWORK, callback_baitr<Real>* callback, Integer& INFO);

template
void baitr<Float>(bool SYM, Integer N, Integer K, Integer NP, Integer& KB,
        Float TOL, Float* X, Integer LDX, bool INIT, Float* RESID, Integer LDR, Float& NORMR, Float* U, Integer LDU, 
        Float* TAU, Float* H_NS, Float* H_S, Integer LDH, Float& SCALE, Float& SUMSQ, Float* T, Integer LDT,
        Integer& NT, Float* WORK, Integer LWORK, callback_baitr<Float>* callback, Integer& INFO);

template
void baitr<Complex>(bool SYM, Integer N, Integer K, Integer NP, Integer& KB,
        Real TOL, Complex* X, Integer LDX, bool INIT, Complex* RESID, Integer LDR, Real& NORMR, Complex* U, Integer LDU, 
        Complex* TAU, Complex* H_NS, Complex* H_S, Integer LDH, Real& SCALE, Real& SUMSQ, Complex* T, Integer LDT,
        Integer& NT, Complex* WORK, Integer LWORK, callback_baitr<Complex>* callback, Integer& INFO);

template
void baitr<Float_complex>(bool SYM, Integer N, Integer K, Integer NP, Integer& KB,
        Float TOL, Float_complex* X, Integer LDX, bool INIT, Float_complex* RESID, Integer LDR, Float& NORMR,
        Float_complex* U, Integer LDU, Float_complex* TAU,  Float_complex* H_NS,  Float_complex* H_S, Integer LDH, 
        Float& SCALE, Float& SUMSQ, Float_complex* T, Integer LDT, Integer& NT,
        Float_complex* WORK, Integer LWORK, callback_baitr<Float_complex>* callback, Integer& INFO);

}}