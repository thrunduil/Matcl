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
#include "matcl-blas-ext/lapack_ext/blas_ext.h"

#include <algorithm>

namespace matcl { namespace lapack
{

/*
*  Purpose
*  =======
*
*  DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
*  matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
*
*  where either
*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
*  conjugate eigenvalues.
*
*  Arguments
*  =========
*
*  A       (input/output) DOUBLE PRECISION
*  B       (input/output) DOUBLE PRECISION
*  C       (input/output) DOUBLE PRECISION
*  D       (input/output) DOUBLE PRECISION
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1R    (output) DOUBLE PRECISION
*  RT1I    (output) DOUBLE PRECISION
*  RT2R    (output) DOUBLE PRECISION
*  RT2I    (output) DOUBLE PRECISION
*          The real and imaginary parts of the eigenvalues. If the
*          eigenvalues are a complex conjugate pair, RT1I > 0.
*
*  CS      (output) DOUBLE PRECISION
*  SN      (output) DOUBLE PRECISION
*          Parameters of the rotation matrix.
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
lapack::lanv2(V& A, V& B, V& C, V& D, V& RT1R, V& RT1I, V& RT2R, V& RT2I, V& CS, V& SN)
{
    const V ZERO    = V(0.0);
    const V ONE     = V(1.0);
    const V HALF    = V(0.5);
    const V MULTPL  = V(4.0);
    const V EPS     = lapack::lamch<V>("P");

    if (C == ZERO )
    {
        CS      = ONE;
        SN      = ZERO;
        goto lab_10;
    }
    else if ( B == ZERO )
    {
        // Swap rows and columns
        CS      = ZERO;
        SN      = ONE;
        V TEMP  = D;
        D       = A;
        A       = TEMP;
        B       = -C;
        C       = ZERO;
        goto lab_10;
    }
    else if ( (A - D) == ZERO && copysign(ONE,B) != copysign(ONE, C) )
    {
        CS      = ONE;
        SN      = ZERO;
        goto lab_10;
    }
    else
    {
        V TEMP      = A - D;
        V P         = HALF * TEMP;
        V BCMAX     = std::max( abs(B), abs(C) );
        V BCMIS     = std::min( abs(B), abs(C) ) * copysign(ONE, B) * copysign(ONE, C );
        V SCALE     = std::max( abs(P), BCMAX );
        V Z         = ( P / SCALE ) * P + ( BCMAX / SCALE ) * BCMIS;

        // If Z is of the order of the machine accuracy, postpone the
        // decision on the nature of eigenvalues
        if ( Z >= MULTPL * EPS )
        {
            // d_type eigenvalues. Compute A and D.
            Z       = P + copysign(std::sqrt(SCALE) * std::sqrt(Z), P);
            A       = D + Z;
            D       = D - ( BCMAX / Z ) * BCMIS;

            // Compute B and the rotation matrix
            V TAU   = lapack::lapy2( C, Z );
            CS      = Z / TAU;
            SN      = C / TAU;
            B       = B - C;
            C       = ZERO;
        }
        else
        {
            // Complex eigenvalues, or real (almost) equal eigenvalues.
            // Make diagonal elements equal.
            V SIGMA = B + C;
            V TAU   = lapack::lapy2( SIGMA, TEMP );
            CS      = std::sqrt( HALF * (ONE + abs(SIGMA) / TAU ) );
            SN      = - (P / (TAU * CS) ) * copysign(ONE, SIGMA);

            // Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
            //         [ CC  DD ]   [ C  D ] [ SN  CS ]

            V AA    = A*CS + B*SN;
            V BB    = B*CS - A*SN;
            V CC    = C*CS + D*SN;
            V DD    = D*CS - C*SN;

            // Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
            //         [ C  D ]   [-SN  CS ] [ CC  DD ]
            A       = AA*CS + CC*SN;
            B       = BB*CS + DD*SN;
            C       = CC*CS - AA*SN;
            D       = DD*CS - BB*SN;

            TEMP    = HALF * ( A+D );
            A       = TEMP;
            D       = TEMP;

            if (C != ZERO )
            {
                if (B != ZERO )
                {
                    if( copysign(ONE, B) == copysign(ONE, C) )
                    {
                        // d_type eigenvalues: reduce to upper triangular form
                         V SAB  = std::sqrt( abs(B) );
                         V SAC  = std::sqrt( abs(C) );
                         P      = copysign(SAB * SAC, C);
                         TAU    = ONE / std::sqrt( abs( B+C ) );
                         A      = TEMP + P;
                         D      = TEMP - P;
                         B      = B - C;
                         C      = ZERO;
                         V CS1  = SAB*TAU;
                         V SN1  = SAC*TAU;
                         TEMP   = CS*CS1 - SN*SN1;
                         SN     = CS*SN1 + SN*CS1;
                         CS     = TEMP;
                    };
                }
                else
                {
                    B           = -C;
                    C           = ZERO;
                    TEMP        = CS;
                    CS          = -SN;
                    SN          = TEMP;
                };
            };
        };

    };

  lab_10:

    // Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
    RT1R            = A;
    RT2R            = D;

    if ( C == ZERO )
    {
        RT1I        = ZERO;
        RT2I        = ZERO;
    }
    else
    {
        RT1I        = std::sqrt( abs( B ) ) * std::sqrt( abs( C ) );
        RT2I        = -RT1I;
    };
};

template BLAS_EXT_EXPORT void 
lapack::lanv2(d_type& A, d_type& B, d_type& C, d_type& D, 
           d_type& RT1R, d_type& RT1I, d_type& RT2R, d_type& RT2I, d_type& CS, d_type& SN);

template BLAS_EXT_EXPORT void 
lapack::lanv2(s_type& A, s_type& B, s_type& C, s_type& D, 
           s_type& RT1R, s_type& RT1I, s_type& RT2R, s_type& RT2I, s_type& CS, s_type& SN);

}};
