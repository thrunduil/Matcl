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
#include <algorithm>

namespace matcl { namespace lapack
{

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
lapack::trsen_cond(const char* JOB, i_type N, i_type M, const V* T, i_type LDT, 
           typename details::real_type<V>::type& S, typename details::real_type<V>::type& SEP, 
           V* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO )
{
    using VR = typename details::real_type<V>::type;
    bool is_V_real  = details::is_complex<V>::value == false;

    i_type N1, N2, NN, LWMIN, LIWMIN;

    // Decode and test the input parameters
    bool WANTBH     = JOB[0] == 'B' || JOB[0] == 'b';
    bool WANTS      = JOB[0] == 'E' || JOB[0] == 'e' ||  WANTBH;
    bool WANTSP     = JOB[0] == 'V' || JOB[0] == 'v' || WANTBH;

    const char* trans_conj_char = (is_V_real)? "T" : "C";

    INFO            = 0;
    bool LQUERY     = ( LWORK == -1 );

    if (WANTS == false &&  WANTSP == false)
        INFO = -1;
    else if ( N < 0 )
        INFO = -2;
    else if ( M < 0 || M > N)
        INFO = -3;
    else if ( LDT < std::max(1,N ) )
        INFO = -5;
    else
    {
        // Set M to the dimension of the specified invariant subspace,
        // and test LWORK and LIWORK.
        if (M < N)
        {
            if (T[M+1-1 + (M-1)*LDT] != V(0))
                M = M + 1;
        };

        N1      = M;
        N2      = N - M;
        NN      = N1*N2;

        if (WANTSP)
        {
            LWMIN   = std::max( 1, 2*NN );
            LIWMIN  = std::max( 1, NN );
        }
        else if (WANTS)
        {
            LWMIN   = std::max( 1, NN );
            LIWMIN  = 1;
        };

        if( LWORK < LWMIN && !LQUERY )
            INFO = -9;
        else if( LIWORK < LIWMIN && !LQUERY )
            INFO = -11;

    };

    if( INFO == 0 )
    {
        WORK[1-1]   = V(VR(LWMIN));
        IWORK[1-1]  = (i_type)LIWMIN;
    };

    if (INFO != 0 )
         return;
    else if( LQUERY )
         return;

    //Quick return if possible.
    if ( M == N || M == 0 )
    {
        if ( WANTS )
            S   = 1.0;

        if ( WANTSP )
        {
            VR RWORK;
            SEP     = lapack::lange("1", N, N, T, LDT, &RWORK );
        }

        goto lab_40;
    };    

    if ( WANTS )
    {
        // Solve Sylvester equation for R:
        //
        //  T11*R - R*T22 = scale*T12

        VR SCALE;
        i_type IERR;

        lapack::lacpy<V>( "F", N1, N2, &T[1-1+(N1+1-1)*LDT], LDT, WORK, N1 );
        lapack::trsyl<V>( "N", "N", -1, N1, N2, T, LDT, &T[N1+1-1 + (N1+1-1)*LDT],
                        LDT, WORK, N1, &SCALE, &IERR );

        // Estimate the reciprocal of the condition number of the cluster
        // of eigenvalues.
        VR RWORK;
        VR RNORM = lapack::lange( "F", N1, N2, WORK, N1, &RWORK );
        if ( RNORM == V(0.0) )
            S = 1.0;
        else
            S = SCALE / ( sqrt( SCALE*SCALE / RNORM+RNORM )* sqrt( RNORM ) );
    };

    if ( WANTSP )
    {
        // Estimate sep(T11,T22).
        VR EST      = 0.0;
        i_type KASE = 0;
        i_type IERR;
        i_type ISAVE[3] = {0};
        VR SCALE;

      lab_30:
   
         lapack::lacn2(NN, &WORK[NN+1-1], WORK, IWORK, &EST, &KASE, ISAVE );

         if ( KASE != 0 )
         {
            if ( KASE == 1 )
            {
                // Solve  T11*R - R*T22 = scale*X.
                lapack::trsyl( "N", "N", -1, N1, N2, T, LDT, &T[N1+1-1 +(N1+1-1)*LDT], 
                              LDT, WORK, N1, &SCALE, &IERR );
            }
            else
            {
                // Solve T11**T*R - R*T22**T = scale*X.
                lapack::trsyl(trans_conj_char, trans_conj_char, -1, N1, N2, T, LDT, &T[N1+1-1 +(N1+1-1)*LDT], LDT, 
                              WORK, N1, &SCALE, &IERR );
            };

            goto lab_30;
         };

         SEP = SCALE / EST;
    };

  lab_40:

    WORK[1-1]   = V(VR(LWMIN));
    IWORK[1-1]  = LIWMIN;
};

template BLAS_EXT_EXPORT
void trsen_cond<s_type>(const char* JOB, i_type N, i_type M, const s_type* T, i_type LDT, s_type& S, s_type& SEP, 
           s_type* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO );

template BLAS_EXT_EXPORT
void trsen_cond<d_type>(const char* JOB, i_type N, i_type M, const d_type* T, i_type LDT, d_type& S, d_type& SEP, 
           d_type* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO );

template BLAS_EXT_EXPORT
void trsen_cond<c_type>(const char* JOB, i_type N, i_type M, const c_type* T, i_type LDT, s_type& S, s_type& SEP, 
           c_type* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO );

template BLAS_EXT_EXPORT
void trsen_cond<z_type>(const char* JOB, i_type N, i_type M, const z_type* T, i_type LDT, d_type& S, d_type& SEP, 
           z_type* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO );

}};
