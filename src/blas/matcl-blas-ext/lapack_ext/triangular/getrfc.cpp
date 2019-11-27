/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"


namespace matcl { namespace lapack
{

template<class T>
void lapack::getrfc(i_type M, i_type N, T *A, i_type LDA, i_type *IPIV, i_type *IQIV, 
                    typename lapack::details::real_type<T>::type TOL, i_type *INFO)
{
//
// =====================================================================
    --A;
    --IPIV;
    --IQIV;

    using TR = typename lapack::details::real_type<T>::type;

    //  Test the input parameters.
    *INFO = 0;

    if ( M < 0 )
        *INFO = -1;
    else if ( N < 0 )
        * INFO = -2;
    else if( LDA < lapack::maximum((i_type)1,M) )
        *INFO = -4;

    if( *INFO != 0 )
        return;

    // Quick return if possible
    if( M == 0 || N == 0)
    {
        *INFO = 0;
        return;   
    };

    i_type MN   = lapack::minimum(M, N);

    TR TOLV     = TOL;
    TR NORM_F   = TR(1.0);

    if (TOLV < 0)
    {
        NORM_F  = lapack::lange<T>("F", M, N, A+1, LDA, nullptr);
        NORM_F  = NORM_F / sqrt(TR(MN));

        TR EPS  = lapack::lamch<TR>("Eps");
        TOL     = TR(10.0) * EPS * NORM_F;
    };    

    i_type  I, J, K, L, LAST, IMAX, JMAX, JLAST, JNEW, IDA1, IDA2, RANK;
    T       V;
    TR      AIJMAX, AJMAX;

    bool is_sing;
    LAST = N;
    RANK = 0;

    if (NORM_F == TR(0.0))
    {
        LAST    = 0;
        goto exit_label;
    };

    //-----------------------------------------------------------------
    //    Start of elimination loop.
    //-----------------------------------------------------------------
    for(K = 1; K <= N; K++) 
    {
        // Find the biggest aij in row imax and column jmax.
        AIJMAX  = 0;
        IMAX    = K;
        JMAX    = K;
        JLAST   = LAST;
        is_sing = true;

        for(J = K; J <= JLAST; J++) 
        {
            L       = amax(M-K+1,A+K+(J-1)*LDA,1)+K-1;
            AJMAX   = abs(A[L+(J-1)*LDA]);

            if(AJMAX <= TOLV) 
            {
                //========================================================
                //    Do column interchange, changing old column to zero.
                //    Reduce  "last"  and try again with same j.
                //========================================================
                JNEW        = IQIV[LAST];
                IQIV[LAST]  = IQIV[J];
                IQIV[J]     = JNEW;

                if (J != LAST)
                {
                    for(I = 1; I <= K-1; I++) 
                    {
                        IDA1    = I+(LAST-1)*LDA;
                        IDA2    = I+(J-1)*LDA;
                        V       = A[IDA1];
                        A[IDA1] = A[IDA2];
                        A[IDA2] = V;
                    }
                    for(I = K; I <= M; I++) 
                    {
                        IDA1            = I+(LAST-1)*LDA;
                        V               = A[IDA1];
                        A[IDA1]         = 0;
                        A[I+(J-1)*LDA]  = V;
                    }
                }
                else
                {
                    for(I = K; I <= M; I++) 
                    {
                        IDA1            = I+(LAST-1)*LDA;
                        A[IDA1]         = 0;
                    }
                };
                --LAST;

                if(J<=LAST)
                {
                    --J;
                    continue;
                };

                break;
            };

            is_sing = false;
            //  Check if this column has biggest aij so far.
            if(AIJMAX < AJMAX) 
            {
                AIJMAX  = AJMAX;
                IMAX    = L;
                JMAX    = J;
            }
            if(J >= LAST)
                break;
        };

        if (is_sing == true)
            goto exit_label;

        ++RANK;
        IPIV[K] = IMAX;

        if(JMAX != K) 
        {
            //==========================================================
            //    Do column interchange (k and jmax).
            //==========================================================
            JNEW        = IQIV[JMAX];
            IQIV[JMAX]  = IQIV[K];
            IQIV[K]     = JNEW;

            IDA1        = 1+(JMAX-1)*LDA;
            IDA2        = 1+(K-1)*LDA;
            for(I = 1; I <= M; ++I,++IDA1,++IDA2) 
            {
                V       = A[IDA1];
                A[IDA1] = A[IDA2];
                A[IDA2] = V;
            }
        }
        if(M>K) 
        {
            //===========================================================
            //       Do row interchange if necessary.
            //===========================================================
            if(IMAX!=K) 
            {
                IDA1        = IMAX;
                IDA2        = K;
                for (I = 1; I <= K; ++I,IDA1+=LDA,IDA2+=LDA)
                {
                    V       = A[IDA1];
                    A[IDA1] = A[IDA2];
                    A[IDA2] = V;
                };
            }
            //===========================================================
            //    Compute multipliers.
            //    Do row elimination with column indexing.
            //===========================================================
            V = TR(1.)/A[K+(K-1)*LDA];
            scal(M-K,V,A+K+1+(K-1)*LDA,1);

            for(J = K+1; J <= LAST; J++) 
            {
                IDA1    = IMAX+(J-1)*LDA;
                V       = A[IDA1];
                if(IMAX!=K) 
                {
                    IDA2    = K+(J-1)*LDA;
                    A[IDA1] = A[IDA2];
                    A[IDA2] = V;
                }
                axpy(M-K,-V,A+K+1+(K-1)*LDA,1, A+K+1+(J-1)*LDA,1);
            }
        }
        else
            break;

        if(K>=LAST)
            break;
    }

  exit_label:

    *INFO = RANK;
    // Set ipvt(*) for singular rows.
    for(K = LAST+1; K <= M; K++)
        IPIV[K] = K;

    return;
};

template void BLAS_EXT_EXPORT
getrfc<s_type>(i_type M, i_type N, s_type *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
                       s_type TOLV, i_type *INFO );

template void BLAS_EXT_EXPORT
getrfc<d_type>(i_type M, i_type N, d_type *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
                       d_type TOLV, i_type *INFO );

template void BLAS_EXT_EXPORT
getrfc<z_type>(i_type M, i_type N, z_type *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
                       d_type TOLV, i_type *INFO );

template void BLAS_EXT_EXPORT
getrfc<c_type>(i_type M, i_type N, c_type *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
                       s_type TOLV, i_type *INFO );

};};