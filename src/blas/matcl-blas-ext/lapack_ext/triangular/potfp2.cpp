/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include <iostream>

namespace matcl { namespace lapack
{

/*
* C++ version of:
*
* function LEV2PCHOL
* Craig Lucas, University of Manchester. January, 2004
*/
template<class T>
void lapack::potfp2(const char * uplo, i_type N, T *A, i_type LDA, i_type* PIV, 
                    i_type& RANK, typename details::real_type<T>::type TOL, 
                    typename details::real_type<T>::type* WORK, i_type& INFO)
{
    using real_type = typename details::real_type<T>::type;
    const real_type ZERO = real_type();
    const real_type ONE  = 1.;

    // Test the input parameters
    INFO        = 0;
    bool UPPER  = (uplo[0] == 'U' || uplo[0] == 'u');
    bool LOWER  = (uplo[0] == 'L' || uplo[0] == 'l');

    if( UPPER == false && LOWER == false)
        INFO    = -1;
    else if( N < 0 )
        INFO    = -2;
    else if( LDA < lapack::maximum( (i_type)1, N ) )
        INFO    = -4;
    
    if( INFO != 0 )
        return;

    // Quick return if possible
    if( N == 0 )
    {
        RANK = 0;
        return;
    };

    // Initialize PIV
    for (i_type i = 0; i < N; ++i)
        PIV[i] = i + 1;

    // Get unit roundoff
    real_type U     = lamch<real_type>( "E" );

    real_type TOLV      = TOL;
    real_type NORM_F    = real_type(1.0);

    if (TOLV < real_type(0.0))
    {
        NORM_F          = lapack::lange<T>("F", N, N, A, LDA, nullptr);
        NORM_F          = NORM_F / sqrt(real_type(N));
        TOLV            = real_type(10.0) * U * NORM_F;
    };

    // Compute stopping value if not supplied
    real_type DSTOP     = TOLV;

    if ( NORM_F == ZERO )
    {
        RANK = 0;
        return;
    };

    i_type PVT          = lapack::amax(N, A, LDA+1);
    PVT                 = PVT - 1;
    T AJJ               = A[PVT + PVT * LDA];
    real_type AJJ_r     = lapack::real(AJJ);
    real_type AJJ_a     = lapack::abs(AJJ);

    if( AJJ_r <= DSTOP)
    {
        RANK = 0;
        return;
    };

    // Set first half of WORK to zero, holds dot products
    for (i_type P = 0; P < N; ++P)
        WORK[P] = ZERO;

    i_type J;
    real_type DTEMP;

    if ( UPPER )
    {
        // Compute the Cholesky factorization P' * A * P = U' * U
        for(J = 0; J <  N; ++J)
        {
            // Find pivot, test for exit, else swap rows and columns
            // Update dot products, compute possible pivots which are
            // stored in the second half of WORK
            for (i_type P = J; P < N; ++P)
            {
                if ( J > 0 )
                {
                    T tmp = A[J-1 + P*LDA];
                    WORK[P] = WORK[P] + lapack::abs2(tmp);
                };
                WORK[N+P] = lapack::real(A[P + P*LDA]) - WORK[P];
            };

            if ( J > 0 )
            {
                i_type ITEMP = lapack::amax( N-J, WORK + N+J, 1) - 1;
                PVT         = ITEMP + J;
                AJJ         = WORK[N+PVT];

                AJJ_r       = lapack::real(AJJ);
                AJJ_a       = lapack::abs(AJJ);

                if( AJJ_r <= DSTOP)
                {
                    A[J + J*LDA] = AJJ_r;
                    goto exit_label;
                };
            };

            if( J != PVT )
            {
                // Pivot OK, so can now swap pivot rows and columns
                A[PVT + PVT*LDA] = A[J + J*LDA];
                lapack::swap(J, A + J*LDA, 1, A + PVT*LDA, 1 );
                lapack::swap(N-PVT-1, A + J + (PVT+1)*LDA, LDA, A + PVT + (PVT+1)*LDA, LDA );
                lapack::swap(PVT-J-1, A + J + (J+1)*LDA, LDA, A + J+1 + PVT*LDA, 1 );
                lapack::lacgv<T>(PVT-J, A + J + (J+1)*LDA, LDA);
                lapack::lacgv<T>(PVT-J, A + J+1 + PVT*LDA, 1);

                // Swap dot products and PIV
                DTEMP           = WORK[J];
                WORK[J]         = WORK[PVT];
                WORK[PVT]       = DTEMP;

                i_type ITEMP    = PIV[PVT];
                PIV[PVT]        = PIV[J];
                PIV[J]          = ITEMP;
            };

            AJJ_r           = std::sqrt(AJJ_r);
            A[J + J*LDA]    = AJJ_r;

            // Compute elements J+1:N of row J
            if ( J  < N )
            {
                lapack::lacgv<T>(J, A + J*LDA, 1);
                lapack::gemv<T>( "Trans", J, N-J - 1, -ONE, A + (J+1)*LDA, LDA,
                            A + J*LDA, 1, ONE, A + J + (J+1)*LDA, LDA );
                lapack::lacgv<T>(J, A + J*LDA, 1);

                lapack::scal<T>( N-J - 1, ONE/AJJ_r, A + J + (J+1)*LDA, LDA );
            };
        };
    }
    else
    {
        // Compute the Cholesky factorization P' * A * P = L * L'

        for ( J = 0; J < N; ++J)
        {
            // Find pivot, test for exit, else swap rows and columns
            // Update dot products, compute possible pivots which are
            // stored in the second half of WORK

            for (i_type P = J; P < N; ++P)
            {
                if ( J > 0 )
                {
                    T tmp = A[P + (J-1)*LDA];
                    WORK[P] = WORK[P] + lapack::abs2(tmp);
                };

                WORK[N+P] = lapack::real(A[P + P*LDA]) - WORK[P];
            };

            if( J > 0 )
            {
                i_type ITEMP    = lapack::amax( N-J, WORK + N+J, 1) - 1;
                PVT             = ITEMP + J;
                AJJ             = WORK[N+PVT];

                AJJ_r     = lapack::real(AJJ);
                AJJ_a     = lapack::abs(AJJ);

                if( AJJ_r <= DSTOP )
                {
                    A[J + J*LDA] = AJJ_r;
                    goto exit_label;
                };

            }
            if ( J != PVT )
            {
                // Pivot OK, so can now swap pivot rows and columns
                A[PVT + PVT*LDA] = A[J + J*LDA];
                lapack::swap<T>(J, A + J, LDA, A + PVT, LDA );
                lapack::swap<T>(N - PVT - 1, A + PVT+1 + J*LDA, 1, A + PVT+1 + PVT*LDA, 1 );
                lapack::swap<T>(PVT-J-1, A + J+1 + J*LDA, 1, A + PVT + (J+1)*LDA, LDA );

                lapack::lacgv<T>(PVT-J, A + J+1 + J*LDA, 1);
                lapack::lacgv<T>(PVT-J, A + PVT + (J+1)*LDA, LDA);

                // Swap dot products and PIV
                DTEMP           = WORK[J];
                WORK[J]         = WORK[PVT];
                WORK[PVT]       = DTEMP;

                i_type ITEMP    = PIV[PVT];
                PIV[PVT]        = PIV[J];
                PIV[J]          = ITEMP;
            };

            AJJ_r           = std::sqrt(AJJ_r);
            A[J + J*LDA]    = AJJ_r;

            // Compute elements J+1:N of column J

            if ( J+1 < N )
            {
                lapack::lacgv<T>(J, A + J, LDA);
                lapack::gemv<T>("No tran", N -J - 1, J, -ONE, A + J+1, LDA,
                                A + J, LDA, ONE, A + J+1 + J*LDA, 1 );
                lapack::lacgv<T>(J, A + J, LDA);

                lapack::scal<T>( N-J-1, ONE / AJJ_r, A + J+1 + J*LDA, 1 );
            };
        };
    };

    RANK = N;
    return;

  exit_label:

    // Rank is number of steps completed
    RANK = J;
    return;
};

template void BLAS_EXT_EXPORT
potfp2<s_type>(const char * uplo, i_type N, s_type *A, i_type LDA, i_type* PIV, 
                    i_type& RANK, s_type TOL, s_type* WORK, i_type& INFO);

template void BLAS_EXT_EXPORT
potfp2<d_type>(const char * uplo, i_type N, d_type *A, i_type LDA, i_type* PIV, 
                    i_type& RANK, d_type TOL, d_type* WORK, i_type& INFO);

template void BLAS_EXT_EXPORT
potfp2<z_type>(const char * uplo, i_type N, z_type *A, i_type LDA, i_type* PIV, 
                    i_type& RANK, d_type TOL, d_type* WORK, i_type& INFO);

template void BLAS_EXT_EXPORT
potfp2<c_type>(const char * uplo, i_type N, c_type *A, i_type LDA, i_type* PIV, 
                    i_type& RANK, s_type TOL, s_type* WORK, i_type& INFO);

}}
